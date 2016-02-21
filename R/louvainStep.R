#' @title Run a single step of the Louvain algorithm
#'
#' @description
#' In each step each node in the network is moved to a neighboring cluster
#' (or moved to a new cluster including only itself)
#' if it increases the modularity.
#' @param G a symmetric N-by-N numeric matrix representing the weights of edges
#' between the N nodes
#' @param C a numeric vector of length N with cluster assignments
#' @param O a vecotor indicating the order by which the nodes are evaluated
#' @param Q the initial modularity value, since the method evaluates only the
#' change in modularity so as not to calculate the modularity each time.
#' If no initial value is provided, the current modularity value is calculated.
#' @return a list containing the a vector of length N (clusters) with cluster assignments
#' that maximizes modularity, and the new modularity (newQ)
#' @export

louvainStep <- function(G, C, O, Q=modularity(G, C), verbose=TRUE){
  m <- sum(G)
  Kj <- rowSums(G)

  # print progress bar for long processes
  t0 <- proc.time()
  pb <- NULL
  counter <- 0
  for (i in O){
    if (verbose){
      counter = counter + 1
      if (is.null(pb)){
        t1 <- proc.time()
        if ((t1-t0)["elapsed"]>10){
          pb <- txtProgressBar(min=0, max=length(O), initial = counter, style = 3)
        }
      } else {
        setTxtProgressBar(pb, counter)
      }
    }

    reassign <- setdiff(unique(C[G[i, ] > 0]), C[i])
    if (length(reassign)==0)
      # return(list(clusters=C, newQ=Q))
      next
    deltaRemove <- deltaModularity(G, C, i)
    deltaQ <- sapply(X = reassign, FUN = function(newC, G, C, i, m, Kj){
      deltaModularity(G, C, i, newC, m, Kj)
    }, G, C, i, m, Kj)
    if ((max(deltaQ) - deltaRemove) > 0){
      I <- which.max(deltaQ)
      Q <- Q + max(deltaQ) - deltaRemove
      C[i] <- reassign[I]
    } else if (deltaRemove < 0){
      C[i] <- max(C) + 1
      Q <- Q - deltaRemove
    }
  }
  return(list(clusters=C, newQ=Q))
}
