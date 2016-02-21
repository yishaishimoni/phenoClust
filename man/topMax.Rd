\name{topMax}
\alias{topMax}
\alias{bottomMin}
\title{Find the Nth largest / smallest value}
\usage{
  topMax(x, N)
}
\arguments{
  \item{x}{a numeric list of length L}

  \item{N}{the Nth largest / smallest value to be found}
}
\value{
  the Nth largest / smallest value
}
\description{
  find the Nth largest / smallest value in an unsorted list
  without sorting. The functions use splitting and
  therefore run in \code{O(N)}.
}

