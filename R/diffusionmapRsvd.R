#' Diffusion maps
#'
#' @param D distance matrix
#' @param numeigen number of diffusion coordinates
#' @param t length of markov chain run
#' @param maxdim default number of coordinates if numeigen NULL
#' @param epsilon value to use in kernel
#' @param rsvd true or false parameter to use rsvd or not
#' @param local true or false parameter to use local scaling or not
#' @param k the number of nearest neighbors to eval for local scaling
#' @param p oversampling param for randomized SVD
#' @param delta for introducing sparsity 
#' @export
diffusionmapRsvd <- function(D,
                             numeigen = 8,
                             t = 0,
                             maxdim = 50,
                             rsvd = TRUE, 
                             local = TRUE,
                             k = 50, 
                             p = 5, 
                             delta = 10^(-6)){
  # make markov matrix
  D <- as.matrix(D)
  n <- dim(D)[1]
  if(local){
    K <- localEpsilon(D, k)
  } else { 
    epsilon <- medianEpsilon(D)
    K <- exp(-D^2 / epsilon)
  }
  # normalize
  v <- sqrt(apply(K, 1, sum))
  A <- K/(v %*% t(v))
  # make A sparse
  ind <- which(A > delta, arr.ind = TRUE)
  A.sparse <- Matrix::sparseMatrix(i = ind[,1], 
                                   j = ind[,2], 
                                   x = A[ind], 
                                   dims = c(n, n))
  # calculate effective number eigenvalues
  if(is.null(numeigen)){
    neff <- min(maxdim + 1,n)
  } else{
    neff <- numeigen + 1
  }
  # get SVD
  if(rsvd){
    output <- randomizedSVD(A, .05*n, p)
  } else{ # use the traditional method
    output <- svd(A, nu = neff)
  }
  eigenvals <- output$d
  psi <- output$u
  # phi <- output$v
  if(t <= 0){
    lambda <- eigenvals[-1]/(1 - eigenvals[-1])
    lambda <- rep(1, n) %*% t(lambda)
    X <- psi[,2:(numeigen+1)]*lambda[,1:(numeigen)]
  }
  else{
    lambda <- eigenvals[-1]^t
    lambda <- rep(1,n) %*% t(lambda)
    X <- psi[,2:(numeigen+1)]*lambda[,1:(numeigen)]
  }
  return(list(X = X, 
              largest_eigen = max(eigenvals), 
              eigenvals = eigenvals[-1]))
}
