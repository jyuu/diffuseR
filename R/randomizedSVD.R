#' Randomized SVD
#'
#' @param A the matrix we want a low rank approx of
#' @param k dominant # eigenfunctions of an svd
#' @param p oversampling parameter
#' @param q number of subspace iterations
#' @param sampling_dist distribution to use with random matrix generation
#' @export
randomizedSVD <- function(A, # make sure this is fed in with # rows > # cols
                          k = 6,
                          p = 5, # set to 5 as Halko paper recommends
                          q = 2,
                          sdist = "normal"){ # set to generate gaussian random matrix
  A <- as.matrix(A)
  n <- ncol(A)
  # calculate the target rank
  l <- round(k) + round(p)
  if(l > n) l <- n
  if(l < 1) stop("Rank not valid!")
  # create a random sampling matrix
  if(sdist == "normal"){
    G <- matrix(rnorm(l*n), n, l)
  } else if(dist == "unif"){
    G <- matrix(rnorm(l*n), n, l)
  } else {
    stop("Sampling distribution is not an option!")
  }
  # building the sampling matrix
  Y <- A %*% G
  rm(G) # free space
  # subspace iterations
  for(i in 1:q) {
    Y <- qr.Q(qr(Y, complete = FALSE), complete = FALSE)
    W <- crossprod(A, Y)
    W <- qr.Q(qr(W, complete = FALSE), complete = FALSE)
    Y <- A %*% W
  }
  rm(W)
  Q <- qr.Q(qr(Y, complete = FALSE) , complete = FALSE)
  rm(Y)
  # project to lower dim subspace
  B <- crossprod(Q, A)
  rsvd_output <- rARPACK::svds(B, k) # quick SVD computation
  #rsvd_output <- svd(B)
  # recover left singular
  rsvd_output$u <- Q %*% rsvd_output$u
  # truncate trailing p terms
  if( dim(rsvd_output$u)[2] > k) {
    rsvd_output$u <- rsvd_output$u[, 1:k]
    rsvd_output$v <- rsvd_output$v[, 1:k]
    rsvd_output$d <- rsvd_output$d[1:k]
  } 
  return(rsvd_output)
}
