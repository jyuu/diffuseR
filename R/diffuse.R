

# Load --------------------------------------------------------------------
library(Matrix)
library(igraph)

# compute eps for kernel --------------------------------------------------
epsc <- function(D, p=.01){
  D <- as.matrix(D)
  n <- dim(D)[1]
  k <- ceiling(p*n)
  k <- ifelse(k<2,2,k) # set min dis of 2 
  D.sort <- apply(D,1,sort) # columns are sorted rows
  dist.knn <- D.sort[(k+1),] # find dist to kth nearest neighbor
  epsilon <- 2*median(dist.knn)^2
  return(epsilon)
}

# sloppy ------------------------------------------------------------------
eps.val <- epsc(D)
# set params 
neigen <- NULL
t <- 0
maxdim <- 50
delta <- 10^-5
# start <- proc.time()
D <- as.matrix(D)
n <- dim(D)[1]
K <- exp(-D^2/eps.val)
v <- sqrt(apply(K,1,sum))
A <- K/(v %*% t(v)) # symmetric graph laplacian 

# make A matrix sparse 
ind <- which(A>delta, arr.ind=TRUE)
Asp <- sparseMatrix(i = ind[,1], j = ind[,2], x = A[ind], dims = c(n,n))

f = function(x, A = NULL){ # matrix multiplication for ARPACK
  as.matrix(A %*% x)
}

cat('Performing eigendecomposition\n') # eigendecomposition
if(is.null(neigen)){ 
  neff = min(maxdim+1,n)  
}else{
  neff =  min(neigen+1, n)
}

# eigendecomposition using ARPACK
decomp = arpack(f,extra=Asp,sym=TRUE,
                options=list(which='LA',nev=neff,n=n,ncv=max(min(c(n,4*neff)))))
psi = decomp$vectors/(decomp$vectors[,1]%*%matrix(1,1,neff))#right ev
phi = decomp$vectors * (decomp$vectors[,1]%*%matrix(1,1,neff))#left ev
eigenvals = decomp$values #eigenvalues
