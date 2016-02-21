\name{deltaModularity}
\alias{deltaModularity}
\title{Change in modularity by moving a node}
\usage{
  deltaModularity(G, C, i, newC = NULL, m = sum(G), Kj=rowSums(G))
}
\arguments{
  \item{G}{a symmetric N-by-N numeric matrix representing
  the weights of edges between the N nodes}

  \item{C}{a vector of length N with numeric cluster
  assignments}

  \item{i}{the index (or name in C) indicating which node
  is to be removed or added}

  \item{newC}{the new assignment for node i. If NULL
  (default) calculates the change from removing i from its
  current assignment}

  \item{m}{the sum of all the edge weights in the network.
  By default the value will be calculated. It is more
  efficient to provide the value if looping over many nodes.}

  \item{Kj}{the ranks of all nodes in G.
  By default the value will be calculated,
  but it is more efficient to provide the value if looping
  over many nodes.}
}
\value{
  If newC==NULL, the change in modularity from removing
  node i from its current assignment. If newC is another
  factor, the change in modularity by adding node i to the
  newC class, without calculating the change due to the
  removal from its current class (this is done for
  efficiency)
}
\description{
  Calculates the change in modularity by adding or removing
  a node to a cluster
}

