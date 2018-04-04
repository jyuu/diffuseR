
# subfunction for randomized svd ------------------------------------------

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
diffusionMapper <- function(D, neigen = 15, t = 0, maxdim = 50, delta = 10^{-5}){
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
  eigenvals <- output[[2]]
  psi <- output[[1]]
  
    
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
    X = psi[,2:(neigen+1)]*lambda[,1:neigen] #diffusion coords. X
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
    
    X = psi[,2:(neigen+1)]*lambda[,1:neigen] #diffusion coords. X
  }
  
  # Let's also calculate the deviation from the original data set we fed into function 
  
  return(X)
}


# try it  -----------------------------------------------------------------
library(diffusionMap)
data("annulus")
some_data <- dist(annulus)
some_results <- diffusionMapper(some_data)


# simulation --------------------------------------------------------------

# matrices to try 
matrix_sims <- seq(1000, 30000, by=1000) # so 30 iterations in total
output <- data.frame(n = matrix_sims, 
                     object_size = NA,
                     time_old = NA, 
                     time_new = NA, 
                     )
for(i in 1:length(matrix_sims)){
  n <- matrix_sims[i] 
  t <- runif(n)^.7*10 
  al <- .15
  bet <- .5
  # ex with noisy spiral 
  x1 <- bet*exp(al*t)*cos(t)+rnorm(n,0,.1)
  y1 <- bet*exp(al*t)*sin(t)+rnorm(n,0,.1)
  plot(x1, y1, pch=20, main="Noisy spiral")
  D <- dist(cbind(x1, y1))
  output$object_size[i] <- object_size(D)
  # the original 
  original <- system.time(dmap_old <- diffusionMap::diffuse(D, neigen=10))
  output$time_old[i] <- as.numeric(original[3])
  new <- system.time(diffusionMapper(D, neigen = 10))
  new$time_new[i] <- as.numeric(new[3])
  # viz 
  
}

