\name{hadamard}
\alias{hadamard}
\title{Create a matrix of hadamard indices from a directed adjacancy matrix}
\usage{
  hadamard(A)
}
\arguments{
  \item{A}{an N-by-N adjacency matrix holding TRUE or 1
  values for edges}
}
\value{
  an N-by-N matrix with the Hadamard coefficient of
  neighbor overlap
}
\description{
  This function converts a directed adjacency matrix to a
  matrix of Hadamard distances based on the overlap of
  neighbors. The implementation is efficiently based on
  matrix multiplications.
}

