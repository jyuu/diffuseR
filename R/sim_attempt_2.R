rm(list=ls())

# subfunction for randomized svd ------------------------------------------
Randomized_SVD <- function(A, k = 6, nu = NULL, nv = NULL, p = 10, q = 2, sdist = "normal"){
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
  # NOTE the computation of Y is vulnerable to round off errors!! 
  
  # project to lower dim subspafce ------------------------------------------
  B <- crossprod(Q, A)
  
  
  # singular value decomp ---------------------------------------------------
  #rsvdObj <- svd(B, nu=nu, nv=nv) # Compute SVD
  Bs <- as(B, "dgCMatrix")
  matmul = function(A, B, transpose = FALSE)
  {
    if(transpose) as.numeric(crossprod(A, B)) else as.numeric(A %*% B);
  }
  rsvdObj <- irlba::irlba(Bs, nu=nu, nv=nv, matmul=matmul) # Compute SVD
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
  
  return(rsvdObj)
}

# subfunction for bandwidth for kernel ------------------------------------

epsilonCompute <- function(D,p=.01){
  
  D = as.matrix(D)
  n = dim(D)[1]
  k = ceiling(p*n)
  k = ifelse(k<2,2,k) # use k of at least 2
  D.sort = apply(D,1,sort)
  dist.knn = D.sort[(k+1),] # find dists. to kth nearest neighbor
  epsilon = 2*median(dist.knn)^2
  
  return(epsilon)
}


# subfunction for diffusion map -------------------------------------------
diffusionMapper <- function(D, neigen = 8, t = 0, maxdim = 50, delta = 10^{-5}){
  # calc parms
  eps.val <- epsilonCompute(D) # .122
  
  # make matrix structure 
  D <- as.matrix(D)
  n <- dim(D)[1]
  K <- exp(-D^2)/eps.val # kernal
  v <- sqrt(apply(K, 1, sum)) # normalize
  A <- K/(v%*%t(v))
  
  # make A sparse
  ind <- which(A > delta, arr.ind = TRUE)
  Asp <- Matrix::sparseMatrix(i = ind[,1], j = ind[,2], x = A[ind], dims=c(n,n))
  
  # # see:http://li.mit.edu/Archive/Activities/Archive/CourseWork/Ju_Li/MITCourses/18.335/Doc/ARPACK/Lehoucq97.pdf
  # f <- function(x, A = NULL){ # matrix multiplication for ARPACK
  #   as.matrix(A %*% x)
  # }
  
  cat('Performing eigendecomposition\n') # eigendecomposition
  if(is.null(neigen)){ 
    neff = min(maxdim+1,n)  
  }else{
    neff = neigen+1
    #neff =  min(neigen+1, n)
  }
  
  # eigendecomposition using ARPACK
  # decomp = igraph::arpack(f,extra=Asp,sym=TRUE,
  #                         options=list(which='LA',nev=neff,n=n,ncv=max(min(c(n,4*neff)))))
  # psi = decomp$vectors/(decomp$vectors[,1]%*%matrix(1,1,neff))#right ev
  # phi = decomp$vectors * (decomp$vectors[,1]%*%matrix(1,1,neff))#left ev
  # old_eigenvals = decomp$values #eigenvalues
  # 
  #  Use Randomized SVD instead
  output <- Randomized_SVD(A, k=neff, q=2)
  eigenvals <- output$d
  #print(eigenvals)
  psi <- output$u
  phi <- output$v
    
  cat('Computing Diffusion Coordinates\n')
  if(t<=0){# use multi-scale geometry
    lambda=eigenvals[-1]/(1-eigenvals[-1])
    lambda=rep(1,n)%*%t(lambda)
    if(is.null(neigen)){#use no. of dimensions corresponding to 95% dropoff
      lam = lambda[1,]/lambda[1,1]
      neigen = min(which(lam<.05)) # default number of eigenvalues
      neigen = min(neigen,maxdim)
      eigenvals = eigenvals[1:(neigen+1)]  
      cat('Used default value:',neigen,'dimensions\n')
    }
    X = psi[,2:(neigen+1)]*lambda[,1:(neigen)] #diffusion coords. X
  }
  else{# use fixed scale t
    lambda=eigenvals[-1]^t
    lambda=rep(1,n)%*%t(lambda)
    
    if(is.null(neigen)){#use no. of dimensions corresponding to 95% dropoff
      lam = lambda[1,]/lambda[1,1]
      neigen = min(which(lam<.05)) # default number of eigenvalues
      neigen = min(neigen,maxdim)
      eigenvals = eigenvals[1:(neigen+1)]  
      cat('Used default value:',neigen,'dimensions\n')
    }
    X = psi[,2:(neigen+1)]*lambda[,1:(neigen)] #diffusion coords. X
  }
  
  # Let's also calculate the deviation from the original data set we fed into function 
  #m <- dim(phi)[2]
  #X_new = X %*% t(phi[,1:30])
  #print(dim(X_new))
  # Difference
  #err <- Matrix::norm(X_new-D)/Matrix::norm(D)
  #print(err)
 # print(eigenvals)
  return(list(X = X, largest_eigen = max(eigenvals)))
}


# try it  -----------------------------------------------------------------
library(diffusionMap)
data("annulus")
some_data <- dist(annulus)
some_results <- diffusionMapper(some_data)
some_results <- diffuse(some_data, neigen = 10)

# where is the bottle neck ?
profvis({
  some_results <- diffuse(some_data, neigen = 10)
})

profvis({
  Randomized_SVD(some_data, k=18, q=2)
})

profvis({
  f = function(x, A = NULL){ # matrix multiplication for ARPACK
    as.matrix(A %*% x)
  }
  igraph::arpack(f,extra=some_data,sym=TRUE,
                  options=list(which='LA',nev=50,n=50,ncv=max(min(c(1000,4*50)))))
})



# simulation --------------------------------------------------------------

# matrices to try 
#matrix_sims <- c(10, 20)
matrix_sims <- seq(1000, 10000, by=1000) # so 30 iterations in total
output <- data.frame(n = matrix_sims, 
                     object_size = NA,
                     time_old = NA, 
                     time_new = NA, 
                     #largest_eigen_actual = NA, 
                     largest_eigen_old = NA, 
                     largest_eigen_new = NA)
for(i in 1:length(matrix_sims)){
  cat("working on sim #: ", i)
  n <- matrix_sims[i] 
  t <- runif(n)^.7*10 
  al <- .15
  bet <- .5
  # ex with noisy spiral 
  x1 <- bet*exp(al*t)*cos(t)+rnorm(n,0,.1)
  y1 <- bet*exp(al*t)*sin(t)+rnorm(n,0,.1)
  plot(x1, y1, pch=20, main="Noisy spiral")
  D <- as.matrix(dist(cbind(x1, y1)))
  output$object_size[i] <- object.size(D)
  # the original 
  original <- system.time(dmap_old <- diffusionMap::diffuse(D, neigen=5))
  output$time_old[i] <- as.numeric(original[3])
  eigen_va <- dim(dmap_old$X)[2]
  new <- system.time(dmap_new <- diffusionMapper(D, neigen = eigen_va))
  output$time_new[i] <- as.numeric(new[3])
  # viz 
  filename <- paste0("~/workspace/diffuser/output/image_",n,".png") 
  png(filename)
  par(mfrow=c(1,2))
  plot(dmap_old$X, main = "Traditional Diffusion")
  plot(dmap_new$X, main = "With Randomized SVD")
  dev.off()
  # store largest eigen 
  #out <- svd(D)$d[1]
  #output$largest_eigen_actual[i] <- out
  output$largest_eigen_old[i] <- max(dmap_old$eigenvals) 
  output$largest_eigen_new[i] <- max(dmap_new$largest_eigen)
}

readr::write_csv(output, path = "~/workspace/diffuser/output/results_040918.csv")
# compare against old 
compare <- readr::read_csv("~/workspace/diffuser/output/results.csv")
# clean up ----------------------------------------------------------------
library(dplyr)
library(xtable)
output <- read.csv("~/workspace/diffuser/output/results.csv")
new_results <- output %>%
  filter(n < 29000) %>%
  mutate(object_size = object_size / 1.25e8, 
         difference = 1000*(largest_eigen_new - largest_eigen_old)) %>%
  select(n, object_size, time_old, time_new, difference)
new_results$n <- round(new_results$n)
new_results$object_size <- round(new_results$object_size, 4)
xtable(new_results)

library(ggplot2)
library(ggthemes)
results <- new_results %>%
  select(n, time_old, time_new) %>%
  rename(`Traditional Time` = time_old, 
         `Adjusted Time` = time_new) %>%
  tidyr::gather(type, time, -n) 
ggplot(results, aes(x = n, y = time, group=type)) + geom_line(aes(linetype=type)) + geom_point(aes(shape=type)) + 
  ggthemes::theme_base() + 
  theme(legend.position="bottom")
