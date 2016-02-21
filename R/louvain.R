#' @title Louvain algorithm
#'
#' @description
#' Implementation of the Louvain method for local maximization of \code{\link{modularity}}.
#' @param G a symmetric N-by-N numeric matrix representing the weights of edges
#' between the N nodes
#' @param C a vector of length N with the initial numeric cluster assignments
#' @param maxreps a numeric value indicating the maximal number of iteration
#' in the Louvain algorithm (default=100)
#' @param verbose a boolean vaue indicating whether to print progress messages
#' (default=TRUE)
#' @return a list containing the a vector of length N (clusters) with cluster assignments
#' that maximizes modularity, and the new modularity (newQ)
#' @seealso modularity
#' @export

louvain <- function(G, C=1:ncol(G), maxreps=100, verbose=TRUE){
  vmsg <- function(x, verb=verbose){
    if (verb)
      message(x)
  }

  # choose the order
  O <- sample(1:ncol(G), replace = F)
  # start with an assignment of each node to its own class
  Q <- modularity(G, C)
  for (i in 1:maxreps){
    L <- louvainStep(G, C, O, Q, verbose)
    vmsg(paste(i, '  : Q =', L$newQ, '| Clusters :', paste(sort(table(L$clusters),decreasing=T), collapse=';')))
    if (L$newQ > Q){
      C <- L$clusters
      C <- match(C,unique(C))
      Q <- L$newQ

      # run the second phase, trying to combine each cluster
      nClust <- length(unique(C))
      metaG <- matrix(nrow=nClust, ncol=nClust)
      for (ci in 1:nClust){
        for (cj in ci:nClust){
          metaG[ci,cj] <- sum(G[C==ci,C==cj])
          metaG[cj,ci] <- metaG[ci,cj]
        }
      }
      metaL <- louvainStep(metaG, C = 1:nClust, O = sample(1:nClust, replace = F))
      if (metaL$newQ - Q > 1e-15){
        tempC <- C
        for (ci in 1:nClust){
          tempC[C==ci] <- metaL$clusters[ci]
        }
        tempC <- match(tempC,unique(tempC))
        C <- tempC
        L$clusters <- C
        L$newQ <- metaL$newQ
        Q <- L$newQ
        vmsg(paste(i+.5, ': Q =', L$newQ, '| Clusters :', paste(sort(table(L$clusters),decreasing=T), collapse=';')))
      }
    } else {
      break
    }
  }
  if (i == maxreps)
    warning("Maximal modularity not found after the allowed number of repetitions")
  vmsg("-----")
  return(L)
}
