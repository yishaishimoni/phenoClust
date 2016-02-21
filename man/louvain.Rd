\name{louvain}
\alias{louvain}
\title{Louvain algorithm}
\usage{
  louvain(G, C=1:ncol(G), maxreps = 100, verbose = TRUE)
}
\arguments{
  \item{G}{a symmetric N-by-N numeric matrix representing
  the weights of edges between the N nodes}

  \item{maxreps}{a numeric value indicating the maximal
  number of iteration in the Louvain algorithm
  (default=100)}

  \item{verbose}{a boolean vaue indicating whether to print
  progress messages (default=TRUE)}

  \item{C}{an optional vector of length N with the initial
  numeric cluster assignments (default is for each sample
  to belong to its own singleton cluster)}
}
\value{
  a list containing the a vector of length N (clusters)
  with cluster assignments that maximizes modularity, and
  the new modularity (newQ)
}
\description{
  Implementation of the Louvain method for local
  maximization of modularity.
}

