# trying to copy randomized matrix decompositions 
# relies variety of matrix decomposition methods that seek to exploit low-rank features
# svd is the most ubi method for dimensionality reduction, data process and compresstion 
# assume the data has low-rank structure

# load --------------------------------------------------------------------
library(rsvd)
library(diffusionMap)

# ex ----------------------------------------------------------------------
# hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
# hilbert <- cbind(H, H)
# H <- hilbert(n=50)
# k=10
# s <- rsvd(H, k=k)
# Hre <- s$u %*% diag(s$d) %*% t(s$v) # matrix approximation
# print(100 * norm( H - Hre, 'F') / norm( H,'F')) # percentage error
# 

# implementation ----------------------------------------------------------
data("annulus")
H <- dist(annulus)


# set params 
A <- H
k <- 6 # for target rank 
nu <- NULL
nv <- NULL
p <- 10 # for oversampling param
q <- 2 
sdist <- "normal"

# change matrix 
A <- as.matrix(A)
m <- nrow(A)
n <- ncol(A)
# flip matix if wide
if (m < n){ 
  A <- t(A) # need to worry about complex inputs ? 
  m <- nrow(A)
  n <- ncol(A)
  flipped <- TRUE
} else {
  flipped <- FALSE
}

# set target rank 
if(is.null(k)) k <- n 
if(k > n) k <-n
if(is.character(k)) stop("Target rank is not valid!")
if(k < 1) stop("Target rank is not valid!")

# set oversampling parameter? 
l <- round(k) + round(p)
if(l > n) l <- n
if(l < 1) stop("Target rank is not valid!")

# check array -- check results if this is an image
if(is.complex(A)) {
  isreal <- FALSE
} else {
  isreal <- TRUE
}

# set number of singular vectors 
if(is.null(nu)) nu <- k
if(is.null(nv)) nv <- k
if(nu < 0) nu <- 0
if(nv < 0) nv <- 0
if(nu > k) nu <- k
if(nv > k) nv <- k
if(flipped==TRUE) {
  temp <- nu
  nu <- nv
  nv <- temp
}


# create a random sampling matrix -----------------------------------------
O <- switch(sdist,
            normal = matrix(stats::rnorm(l*n), n, l),
            unif = matrix(stats::runif(l*n), n, l),
            rademacher = matrix(sample(c(-1,1), (l*n), replace = TRUE, prob = c(0.5,0.5)), n, l),
            stop("Selected sampling distribution is not supported!"))

if(isreal==FALSE) {
  O <- O + switch(sdist,
                  normal = 1i * matrix(stats::rnorm(l*n), n, l),
                  unif = 1i * matrix(stats::runif(l*n), n, l),
                  rademacher = 1i * matrix(sample(c(-1,1), (l*n), replace = TRUE, prob = c(0.5,0.5)), n, l),
                  stop("Selected sampling distribution is not supported!"))
}


# build sample matrix -----------------------------------------------------
Y <- A %*% O # should approximate the range of A 
remove(O)

# orthogonalize Y using QR decomp -----------------------------------------
# q determines the number of subspace iterations 
if (q >0){
  for(i in 1:q) {
    Y <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
    Z <- crossprod(A , Y) # doesn't account for complex case
    Z <- qr.Q( qr(Z, complete = FALSE) , complete = FALSE )
    Y <- A %*% Z
  } 
  remove(Z)
} 

Q <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
remove(Y)


# project to lower dim subspafce ------------------------------------------
B <- crossprod(Q, A)


# singular value decomp ---------------------------------------------------
rsvdObj <- svd(B, nu=nu, nv=nv) # Compute SVD
rsvdObj$d <- rsvdObj$d[1:k] # Truncate singular values

if(nu != 0) rsvdObj$u <- Q %*% rsvdObj$u # Recover left singular vectors


# if flipped because it was wide ------------------------------------------
if(flipped == TRUE) {
  u_temp <- rsvdObj$u
  rsvdObj$u <- rsvdObj$v
  rsvdObj$v <- u_temp
}

# return ------------------------------------------------------------------
if(nu == 0){ rsvdObj$u <- NULL}
if(nv == 0){ rsvdObj$v <- NULL}
