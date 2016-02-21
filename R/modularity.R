#' @title Calculate modularity of a weighted graph matrix
#'
#' @description
#' Modularity is defined as
#' \code{Q=1/m * sum_(i,j) (Aij-ki*kj/m)*delta(ci,cj)},
#' where m is the sum of all edges, ki is the rank of node i,
#' ci is the cluster assignment of node i, and delta is the delta dunction.
#' The implementation here only uses matrix multiplications and is
#' quite efficient. However, when checking the changes in the modularity
#' it is more efficient to use \code{deltaModularity}
#' @param G a symmetric N-by-N numeric matrix representing the weights of edges
#' between the N nodes
#' @param C a vector or factor of length N with numeric or factor assignments
#' @return A numeric value of the the modularity
#' @seealso \code{\link{deltaModularity}}
#' @export

modularity <- function(G, C){
  m=sum(G)
  Ki <- matrix(data = rowSums(G), ncol = ncol(G), nrow = nrow(G))
  Kj <- t(Ki)
  deltaFunction <- outer(X = C, Y = C, FUN = "==")
  Q <- (G - Ki*Kj/m)*deltaFunction/m
  return(sum(Q))
}
