#' @title Change in modularity by moving a node
#'
#' @description
#' Calculates the change in modularity by adding or removing a node to a cluster
#' @param G a symmetric N-by-N numeric matrix representing the weights of edges
#' between the N nodes
#' @param C a vector of length N with numeric cluster assignments
#' @param i the index (or name in C) indicating which node is to be removed or added
#' @param newC the new assignment for node i. If NULL (default) calculates the
#' change from removing i from its current assignment
#' @param m the sum of all the edge weights in the network. By default the value
#' will be calculated. It is more efficient to provide the value if looping over
#' many nodes.
#' @param Kj the ranks of all nodes in G. By default the value will be calculated,
#' but it is more efficient to provide the value if looping over many nodes.
#' @return If newC==NULL, the change in modularity from removing node i from
#' its current assignment. If newC is another factor, the change in modularity
#' by adding node i to the newC class, without calculating the change due to
#' the removal from its current class (this is done for efficiency)
#' @seealso modularity
#' @export

deltaModularity <- function(G, C, i, newC=NULL, m=sum(G), Kj=rowSums(G)){
  if (is.null(newC)){ # removing a node from C
    # removing a solitary node from a cluster does not change the modularity
    if (sum(C == C[i]) == 1)
      return(0)
    newC <- C[i]
    C[i] <- max(C)+1
  }
  I <- C==newC
  Ki <- sum(G[i, ]) # the sum of all the edges connected to node i
  KiIn <- sum(G[i,I]) # the sum of all the edges conneting node i to C
  Kjs <- sum(Kj[I])
  deltaQ <- KiIn/m - Ki*Kjs/m^2
  return(deltaQ*2)
}
