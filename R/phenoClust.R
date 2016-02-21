#' @title  Cluster samples by maximizing graph modularity
#'
#' @description
#' This method clusters sample data using a two-step approach. In the first
#' step it converts the high dimensional data into a knn-nearest neighbor graph,
#' thus connecting samples that have a similar profile.
#' The edges in knn-nearest neighbors graph are given weights according
#' also known as the Hadamard index.
#' In the second step the Louvaine algorithm, which was developed to find
#' communities in social networks, is applied to this graph in order to find
#' samples that share similar spaces in multi-dimensional space.
#'
#' @param X an M-by-N numeric matrix (or any object that can be coerced to a
#' matrix) holding N samples in columns, where each column has M values.
#' @param method is a string that is passed to the \code{\link[amap]{Dist}}
#' function (default="spearman")
#' @param D a dist object. If provided then \code{X} and \code{method} are ignored
#' (default=NULL, and is calculated from \code{X}).
#' @param knn a numeric value indicating the number of nearest neighbors to use
#' in the initial neighborhood (default=NULL uses nclass.Sturges)
#' @param G an N-by-N adjacency matrix on which modularity will be maximized. If
#' provided \code{X}, \code{method}, \code{D}, and \code{knn} are ignored
#' (default=NULL, and is calculated from \code{D})
#' @param repeats a numeric value for how many times to run the Louvain algorithm
#' in order to find a maximal local maximum (default=50).
#' @param verbose a boolean vaue indicating whether to print progress messages
#' (default=TRUE)
#' @param C a numeric vector of length N indicating the initial cluster assignment.
#' Value of NULL (default) is equivalent to \code{C=1:N}.
#' @details
#' It should be noted that lower values of \code{knn} will result in faster runtime.
#' However, values of knn that are too low will result in inconsistent clustering.
#' Similarly, providing an intial value of \code{C} will reduce the runtime but
#' will affect the clustering.
#' @return a list holding three values:
#' \code{C} - a numeric vector of length N indicating the cluster number.
#' \code{G} - an N-by-N matrix representing the Hadamrd graph.
#' \code{Q} - the modularity provided by \code{G} and \code{C}.
#' @references Levine et al, "Data-Driven Phenotypic Dissection of AML Reveals
#' Progenitor-like Cells that Correlate with Prognosis", Cell, Volume 162,
#' Issue 1, 2 July 2015, Pages 184-197
#' @examples
#' data(iris)
#' C <- phenoClust(X=t(iris[,1:4]), verbose=FALSE, method="euclid")
#' print(table(C$C,iris[,5]))
#' @seealso hadamard, modularity, louvain
#' @export

phenoClust <- function(X=NULL, method="spearman", D=NULL, knn=NULL, G=NULL,
                       C=NULL, repeats=50, verbose=TRUE){
  vmsg <- function(x, verb=verbose){
    if (verb)
      message(x)
  }

  if (is.null(G) && is.null(D) && is.null(X)){
    stop("One of X, D, or G must be provided")
  }
  if (is.null(G) && is.null(D)){
    if (!is.numeric(as.matrix(X)))
      stop("X must be a numeric matrix or coercable to a numeric matrix")
    vmsg("Calculating distances")
    D <- amap::Dist(t(X), method = method) # this is a memory and time-consuming part
  }
  if (is.null(G)){
    if (class(D)=="dist"){
      D <- as.matrix(D)
    }
    D <- as.matrix(D)
    diag(D) <- 0 # a sample is always the nearest neighbor to itself

    # for each column identify the knn for nearest neighbors
    M <- nrow(D)
    N <- ncol(D)
    if (is.null(knn))
      knn <- ceiling(N/nclass.Sturges(1:N))

    # find the knn nearest neighbors
    vmsg(paste("Finding the", knn, "nearest samples"))
    nearest <- apply(X=D, MARGIN = 2, FUN=function(x, knn){
      x <= bottomMin(x, knn+1)
    },knn)
    # this provides a matrix whose columns indicate the nearest neighbors of each sample
    # such that nearest[i,j]==1 means that v(j) is a near neighbor of v(i)
    # the knn+1 ensures that the k-nearest-neighbors are on top of choosing the smaple itself

    vmsg("Converting the nearest neighbor network to a Hadamard distance")
    G <- hadamard(nearest)
    diag(G) <- 0
  }

  if (is.null(C)){
    initC <- 1:N
  } else {
    initC <- C
  }

  vmsg("Running the Louvain algorithm")
  Q <- -Inf
  for (t in 1:repeats){
    L <- louvain(G, C=initC, verbose = verbose)
    if (L$newQ > Q){
      Q <- L$newQ
      C <- L$clusters
    }
  }
  return(list(C=C,G=G,Q=Q))
}
