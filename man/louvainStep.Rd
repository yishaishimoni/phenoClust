\name{louvainStep}
\alias{louvainStep}
\title{Run a single step of the Louvain algorithm}
\usage{
  louvainStep(G, C, O, Q = modularity(G, C), verbose = TRUE)
}
\arguments{
  \item{G}{a symmetric N-by-N numeric matrix representing
  the weights of edges between the N nodes}

  \item{C}{a numeric vector of length N with cluster
  assignments}

  \item{O}{a vecotor indicating the order by which the
  nodes are evaluated}

  \item{Q}{the initial modularity value, since the method
  evaluates only the change in modularity so as not to
  calculate the modularity each time. If no initial value
  is provided, the current modularity value is calculated.}

  \item{verbose}{a boolean vaue indicating whether to print
  progress messages and a progression bar (default=TRUE)}
}
\value{
  a list containing the a vector of length N (clusters)
  with cluster assignments that maximizes modularity, and
  the new modularity (newQ)
}
\description{
  In each step each node in the network is moved to a
  neighboring cluster (or moved to a new cluster including
  only itself) if it increases the modularity.
}

