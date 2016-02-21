\name{phenoClust}
\alias{phenoClust}
\title{Cluster samples by maximizing graph modularity}
\usage{
  phenoClust(X=NULL, method="spearman", D=NULL, knn=NULL, G=NULL,
             C=NULL, repeats=50,  verbose=TRUE)
}
\arguments{
  \item{X}{an M-by-N numeric matrix (or any object that can
  be coerced to a matrix) holding N samples in columns,
  where each column has M values. X can also be an object
  of class \code{dist}, in which case the distance is not
  re-calculated, and \code{method} is ignored.}

  \item{method}{is a string that is passed to the
  \code{\link[amap]{Dist}} function (default="spearman")}

  \item{D}{is a dist object. If provided then \code{X} and
  \code{method} are ignored (default=NULL,
  and is calculated from \code{X})}

  \item{knn}{a numeric value indicating the number of nearest
  neighbors to use in the initial neighborhood
  (default=NULL uses nclass.Sturges)}

  \item{G}{an N-by-N adjacency matrix on which modularity will be maximized.
  If provided \code{X}, \code{method}, \code{D}, and \code{knn} are ignored
  (default=NULL, and is calculated from \code{D})}

  \item{C}{an optional vector of length N with the initial
  numeric cluster assignments (default is for each sample
  to belong to its own singleton cluster)}

  \item{repeats}{a numeric value for how many times to run
  the Louvain algorithm in order to find a maximal local
  maximum (default=50).}

  \item{verbose}{a boolean vaue indicating whether to print
  progress messages (default=TRUE)}
}
\value{
  a list holding three values: \code{C} - a numeric vector
  of length N indicating the cluster number. \code{G} - an
  N-by-N matrix representing the Hadamrd graph. \code{Q} -
  the modularity provided by \code{G} and \code{C}.
}
\description{
  This method clusters sample data using a two-step
  approach. In the first step it converts the high
  dimensional data into a knn-nearest neighbor graph, thus
  connecting samples that have a similar profile. The edges
  in knn-nearest neighbors graph are given weights according
  to the Jaccard similarity coefficient of the neighbors of
  each node. In the second step the Louvaine algorithm,
  which was developed to find communities in social
  networks, is applied to this graph in order to find
  samples that share similar spaces in multi-dimensional
  space.
}
\examples{
data(iris)
C <- phenoClust(X=t(iris[,1:4]), verbose=FALSE, method="euclid")
print(table(C$C,iris[,5]))
}
\references{
  Levine et al, "Data-Driven Phenotypic Dissection of AML
  Reveals Progenitor-like Cells that Correlate with
  Prognosis", Cell, Volume 162, Issue 1, 2 July 2015, Pages
  184-197
}

