\name{modularity}
\alias{modularity}
\title{Calculate modularity of a weighted graph matrix}
\usage{
  modularity(G, C)
}
\arguments{
  \item{G}{a symmetric N-by-N numeric matrix representing
  the weights of edges between the N nodes}

  \item{C}{a vector or factor of length N with numeric or
  factor assignments}
}
\value{
  A numeric value of the the modularity
}
\description{
  Modularity is defined as \code{Q=1/(2m) * sum_(i,j)
  (Aij-ki*kj/(2m)delta(ci,cj)}. The implementation here
  only uses matrix multiplications and is quite efficient.
  However, when checking the changes in the modularity it
  is more efficient to use \code{deltaModularity}
}
\seealso{
  \code{\link{deltaModularity}}
}

