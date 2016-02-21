#' @title Create a matrix of hadamard indices from a directed adjacancy matrix
#'
#' @description
#' This function converts a directed adjacency matrix to a matrix of Hadamard distances
#' based on the overlap of neighbors. The implementation is efficiently based
#' on matrix multiplications.
#' @param A an N-by-N adjacency matrix holding TRUE or 1 values for edges
#' @return an N-by-N matrix with the Hadamard coefficient of neighbor overlap
#' @export

hadamard <- function(A){
  # multiplying A gives the size of the intersection of neighbors
  common <- t(A) %*% A
  # now we build a matrix of the rank of each node
  ranks <- outer(diag(common),rep(1,nrow(A)))
  neighborUnion <- ranks + t(ranks) - common
  # finally, calculate the Hadamard index
  return(common/neighborUnion)
}
