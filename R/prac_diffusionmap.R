rm(list=ls())
# Load library ------------------------------------------------------------
library(diffusionMap)
library(scatterplot3d)

# Ex 1 --------------------------------------------------------------------
t <- seq(-pi, pi, .001)
x <- cbind(cos(t), sin(t))
y <- cos(3*t) + rnorm(length(t), 0, .1)
tcol <- topo.colors(32)
colvec <- floor((y-min(y))/(max(y)-min(y))*32)
colvec[colvec==0] <- 1
scatterplot3d(x[,1], x[,2], y, highlight.3d = T, pch=20, 
              xlab = NULL, ylab=NULL, main =NULL, grid=FALSE, 
              tick.marks = FALSE,
              axis=FALSE)
#D <- as.matrix(dist(cbind(x,y)))
D <- as.matrix(dist(x))
#scatterplot3d(x[,1], x[,2], y, color=tcol[colvec], pch=20, 
              main="Cosine function supported on circle", angle=55,
              cex.main=2, col.axis="gray", cex.symbols=2, cex.lab=2,
              xlab=expression("x"[1]), ylab=expression("x"[2]), zlab="y")
# Now use diffusion map 
D <- as.matrix(dist(x))
AR <- adapreg(D, y, mmax=5, nfolds=2, nrep=2)
print(paste("optimal model size:",AR$mopt,"; optimal epsilon:",
            round(AR$epsopt,4),"; min. CV risk:",round(AR$mincvrisk,5)))
plot(y,AR$y.hat,ylab=expression(hat("y")),cex.lab=1.5,cex.main=1.5,
     main="Predictions")
abline(0,1,col=2,lwd=2)
# plot.dmap(AR$y.hat)
# print.dmap(AR)
# diffussionK

# General diffuse ---------------------------------------------------------
# for data parameterization, that exploits the natural geometry of a data 
# set, diffusion map uses local interactions between data points, propagated
# to larger scales, to construct a global representation of the data.
# computation requires singular value decomposition of the normalized
# graph laplacian. the operation is optimized for speed by exploiting the
# sparseness of the graph laplacian and by using ARPACK for fast matrix
# decomposition, increasing the sparseness parameter, delta, will speed
# up the algo
n <- 000 
t <- runif(n)^.7*10 
al <- .15
bet <- .5
# ex with noisy spiral 
x1 <- bet*exp(al*t)*cos(t)+rnorm(n,0,.1)
y1 <- bet*exp(al*t)*sin(t)+rnorm(n,0,.1)
plot(x1, y1, pch=20, main="Noisy spiral")
D <- dist(cbind(x1, y1))
system.time(dmap <- diffuse(D, neigen=10))
# results 
# Performing eigendecomposition
# Computing Diffusion Coordinates
# Elapsed time: 1.005 seconds
# user  system elapsed 
# 0.977   0.029   1.005 

par(mfrow=c(2,1))
plot(t, dmap$X[,1], pch=20, axes=FALSE, xlab="spiral parameter",
     ylab="1st diffusion coefficient")
box()
plot(1:10, dmap$eigenmult, typ="h", xlab="diffusion map dimension", 
     ylab="eigen-multipliers")


# max matrix size ---------------------------------------------------------
P <- 20000
temp <- matrix(rep(0, P*P), nrow=P)
temp_size <- object.size(temp) 
print(temp_size/(1.25e8)) # 25.6 bytes
# try it again
n <- 20000 
t <- runif(n)^.7*10 
al <- .15
bet <- .5
# ex with noisy spiral 
x1 <- bet*exp(al*t)*cos(t)+rnorm(n,0,.1)
y1 <- bet*exp(al*t)*sin(t)+rnorm(n,0,.1)
plot(x1, y1, pch=20, main="Noisy spiral")
D <- dist(cbind(x1, y1))
system.time(dmap <- diffuse(D, neigen=10))
# Performing eigendecomposition
# Computing Diffusion Coordinates
# Elapsed time: 105.3 seconds
# user  system elapsed 
# 83.009   9.800 105.308

# ideas -------------------------------------------------------------------
# track time for various sizes of D

# try it with annulus -----------------------------------------------------
rm(list=ls())
data("annulus")
plot(annulus, main="Annulus Data", pch=20, cex=.7)
D <- dist(annulus)
dmap <- diffuse(D, eps.val=.1)
print(dmap)
plot(dmap)

# redo their diffuse function ---------------------------------------------
# based off annulus D
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

# set parms
eps.val <- epsilonCompute(D) # .122
neigen <- NULL
t <- 0
maxdim <- 50
delta <- 10^{-5}

# make matrix structure 
D <- as.matrix(D)
n <- dim(D)[1]
K <- exp(-D^2)/eps.val # kernal
v <- sqrt(apply(K, 1, sum)) # normalize
A <- K/(v%*%t(v))

# make A sparse
ind <- which(A > delta, arr.ind = TRUE)
Asp <- Matrix::sparseMatrix(i = ind[,1], j = ind[,2], x = A[ind], dims=c(n,n))

# see:http://li.mit.edu/Archive/Activities/Archive/CourseWork/Ju_Li/MITCourses/18.335/Doc/ARPACK/Lehoucq97.pdf
f <- function(x, A = NULL){ # matrix multiplication for ARPACK
  as.matrix(A %*% x)
}

cat('Performing eigendecomposition\n') # eigendecomposition
if(is.null(neigen)){ 
  neff = min(maxdim+1,n)  
}else{
  neff =  min(neigen+1, n)
}

# eigendecomposition using ARPACK
decomp = igraph::arpack(f,extra=Asp,sym=TRUE,
                options=list(which='LA',nev=neff,n=n,ncv=max(min(c(n,4*neff)))))
psi = decomp$vectors/(decomp$vectors[,1]%*%matrix(1,1,51))#right ev
phi = decomp$vectors * (decomp$vectors[,1]%*%matrix(1,1,51))#left ev
eigenvals = decomp$values #eigenvalues

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
  lambda=rep(1,n)%*%t(lambda) # sum all the lambda values raised to t 
  
  if(is.null(neigen)){#use no. of dimensions corresponding to 95% dropoff
    lam = lambda[1,]/lambda[1,1]
    neigen = min(which(lam<.05)) # default number of eigenvalues
    neigen = min(neigen,maxdim)
    eigenvals = eigenvals[1:(neigen+1)]  
    cat('Used default value:',neigen,'dimensions\n')
  }
  X = psi[,2:(neigen+1)]*lambda[,1:neigen] #diffusion coords. X
}

# results 
X
psi
phi
eigenval <- eigenvals[-1]
eigenmult <- lambda[1,1:neigen]
neigen
eps.val
