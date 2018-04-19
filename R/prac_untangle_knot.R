rm(list=ls())

# load --------------------------------------------------------------------
#library(diffusionMap)
library(scatterplot3d)
library(destiny)

# trefoil knot ------------------------------------------------------------
n <- 400
t <- runif(n, 0, 2*pi)
x <- (2 + cos(3*t))*cos(2*t)
y <- (2 + cos(3*t))*sin(2*t)
z <- sin(3*t)
# viz
scatterplot3d(x,y,z)

# untangle it -------------------------------------------------------------
D <- as.matrix(dist(cbind(x,y,z)))
esp_range <- seq(0.01, 1, by = .01)
norms <- array(NA, c(length(esp_range)))
for(i in 1:length(esp_range)){
  dm <- diffuse(D, eps.val = esp_range[i])
  plot(dm$X[,1], dm$X[,2])
  norms[i] <- norm(dm$X)
}
# see the choice of the eps value is important 
# show that it is data dependent 

# density -----------------------------------------------------------------
# n <- 1000
# r <- sqrt(runif(n))
# theta <- 2*pi*runif(n)
# # convert to cartesian
# x <- r*cos(theta)
# y <- r*sin(theta)
# # viz
# plot(x,y) 
# # x5 <- x
# # y5 <- y
# x1 <- x
# y1 <- y*(1-cos(pi*x))
# D <- as.matrix(dist(cbind(x1,y1)))
# D <- as.matrix(dist(cbind(x,y)))
# output <- DiffusionMap(D)
# plot(output)


# work with trefoil -------------------------------------------------------
# imputed_data <- D
# k <- 20
# n_local = seq(to = min(k, 7L), length.out = min(k, 3L))
# nn_dist <- t(apply(D, 1, function(row) sort(row)[2:k]))
# 
# sig_mat <- nn_dist[, n_local, drop = FALSE]
# sigma <- rowSums(sig_mat) / length(n_local) / 2 
# eps <- min(sigma)

# finding epsilon ---------------------------------------------------------

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

