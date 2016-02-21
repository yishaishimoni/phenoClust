\docType{package}
\name{phenoGraph}
\alias{phenoGraph}
\alias{phenoGraph-package}
\title{Find clusters by network modularity}
\description{
  This package implements the PhenoGraph algorithm, which
  first creates an adjacancy matrix according to k-nearest
  neighbours, and then uses the Louvain algorithm to
  identify "communities" by maximizing the modularity.

  In general, the only function that needs to be run is
  \code{\link{phenoClust}}. Other functions are available
  as utilities.
}

