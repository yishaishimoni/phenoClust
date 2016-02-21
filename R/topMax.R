#' @title Find the Nth largest / smallest value
#'
#' @description find the Nth largest / smallest value in an unsorted list without sorting.
#' The functions use splitting and therefore run in \code{O(N)}.
#'
#' @param x a numeric list of length L
#' @param N the Nth largest / smallest value to be found
#' @return the Nth largest / smallest value

topMax <- function(x,N){
  L <- length(x)
  if (N > L)
    stop("N cannot be larger than the length of the list")

  repeat{
    if (L==1)
      break
    initialGuess <- x[1]
    toplist <- x[x>initialGuess]
    bottomlist <- x[x<initialGuess]

    topL <- length(toplist)
    bottomL <- length(bottomlist)

    if ((topL < N ) & (L - bottomL >= N)){
      x <- initialGuess
      break
    }

    if (topL >= N) {
      x <- toplist
    } else {
      x <- bottomlist
      N <- N - L + bottomL
    }
    L <- length(x)
  }
  return (x[1])
}

#' @export
bottomMin <- function(x,N){
  topMax(x, length(x) - N + 1)
}
