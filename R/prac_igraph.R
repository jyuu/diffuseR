# learning how to use ARPACK

# load --------------------------------------------------------------------
library(igraph)
library(diffusionMap)

# snippet from diffuser ---------------------------------------------------

# set up 
data("annulus")
D <- as.matrix(dist(annulus))
neigen <- 15
t <- 0
maxdim <- 50
delta <- 10^{-5}

# calc parms
eps.val <- epsilonCompute(D) # .122

# make matrix structure 
n <- dim(D)[1]
K <- exp(-D^2)/eps.val # kernal
v <- sqrt(apply(K, 1, sum)) # normalize
A <- K/(v%*%t(v))

# make A sparse
ind <- which(A > delta, arr.ind = TRUE)
Asp <- Matrix::sparseMatrix(i = ind[,1], j = ind[,2], x = A[ind], dims=c(n,n))

f = function(x, A = NULL){ # matrix multiplication for ARPACK
  as.matrix(A %*% x)
}

cat('Performing eigendecomposition\n') # eigendecomposition
if(is.null(neigen)){ 
  neff = min(maxdim+1,n)  
}else{
  neff =  min(neigen+1, n)
}

# nev: numeric scalar, the number of eigenvalues to be computed 
# n: numeric scolar, dimension of the eigenproblem 
# ncv: number of lanczosa vectors to be generated 
# eigendecomposition using ARPACK
decomp = arpack(f,extra=Asp,sym=TRUE,options=list(which='LA',nev=neff,n=n,ncv=max(min(c(n,4*neff)))))
psi = decomp$vectors/(decomp$vectors[,1]%*%matrix(1,1,neff))#right ev
phi = decomp$vectors * (decomp$vectors[,1]%*%matrix(1,1,neff))#left ev
eigenvals = decomp$values #eigenvalues

compare <- svd(Asp, nu =16)

