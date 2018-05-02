#' Compute the optimal epsilon using the median 
#' @param D distance matrix
#' @param percentage distance of neighbors from all n 
#' @export
medianEpsilon <- function(D,
                          percentage = .01){ # neighbors taken as % of whole 
  n <- dim(D)[1]
  k <- ceiling(percentage*n)
  D.sorted <- apply(D, 2, sort)
  scale <- D.sorted[k+1, ]
  epsilon <- 2*median(scale)^2
  return(epsilon)
}

#' Compute the optimal epsilon using kth nearest neighbor for each pt 
#' @param D distance matrix 
#' @param k kth nearest neighbor 
#' @export
localEpsilon <- function(D,
                         k = 2){
  D.sorted <- apply(D, 2, sort)
  scale <- D.sorted[k+1, ]
  K <- scaledist(D^2,scale)
  return(K)
}
