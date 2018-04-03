
Randomized_SVD <- function(A,k,q){
  
  n <- ncol(A)
  #Stage A
  #1
  gauss_mat <- matrix( rnorm(n*(2*k),mean=0,sd=1), n, 2*k)
  #2
  Y <-  A %*% gauss_mat
  #3
  if( q > 0 ) {
    for( i in 1:q) {
      Y <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
      Z <- crossprod(A , Y )
      Z <- qr.Q( qr(Z, complete = FALSE) , complete = FALSE )
      Y <- A %*% Z
    }#End for
    remove(Z)
  }#End if
  Q <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
  remove(Y)
  
  
  #Stage B
  #1
  B <- t(Q) %*% A
  #2
  svd_B <- svd(B)
  #3
  U <- Q %*% svd_B$u
  diags <- svd_B$d
  
  ret_r_Svd <- list(U,diags)
  return(ret_r_Svd)
}


