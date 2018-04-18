
# load --------------------------------------------------------------------
library(diffusionMap)
library(scatterplot3d)

# params ------------------------------------------------------------------
n <- 1000

# generate points uniformly on disk ---------------------------------------
# # naive
# r <- runif(n, 0, 1)
# theta <- runif(n, 0, 2*pi)
# # convert to x, y
# x <- r*cos(theta)
# y <- r*sin(theta)
# # viz 
# plot(x,y) # too concentrated at the center 
r <- sqrt(runif(n))
theta <- 2*pi*runif(n)
# convert to cartesian
x <- r*cos(theta)
y <- r*sin(theta)
# viz
plot(x,y) 
# use above for x = 5 
x5 <- x
y5 <- y

# define maps h, v --------------------------------------------------------
x1 <- x
y1 <- y*(1-cos(pi*x))
x9 <- x*(1-cos(pi*y))
y9 <- y 
# viz 
plot(x1, y1) # yay, barbell
plot(x9, y9) # yay another barbell

# intermediate deformations -----------------------------------------------
others <- c(2, 3, 4)
for(i in others){
  name_x <- paste0("x", i)
  assign(name_x, (x5-x1)*(i-1)/4 + x1)
  name_y <- paste0("y", i)
  assign(name_y, (y5-y1)*(i-1)/4 + y1)
}
more <- c(6, 7, 8)
for(i in more){
  name_x <- paste0("x", i)
  assign(name_x, (x9-x5)*(i-5)/4 + x5)
  name_y <- paste0("y", i)
  assign(name_y, (y9-y5)*(i-5)/4 + y5)
}
# viz 
plot(x1, y1)
plot(x2, y2)

# compile it  -------------------------------------------------------------
X <- array(NA, c(n, 2, 9))
for(i in 1:9){
  X[,1,i] <- eval(parse(text = paste0("x",i)))
  X[,2,i] <- eval(parse(text = paste0("y",i)))
}
# trajectory for a single point
one <- t(X[100,,])
plot(one)
lines(one)

# compute diffusion map ---------------------------------------------------
c1 <- diffuse(dist(cbind(y1,x1)))
plot(c1)
c1 <- diffuse(dist(cbind(y1,x1)))
plot(c1)

c5 <- diffuse(dist(cbind(y5,x5)))
plot(c5)
c9 <- diffuse(dist(cbind(y9,x9)))
plot(c9)

# not getting the same figures as the paper... 

# Try doing it manually ---------------------------------------------------
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

temp <- X1[,,1]
D <- as.matrix(dist(temp))
eps <- epsilonCompute(D)
K <- exp(-(D^2)/eps)
W <- sqrt(apply(K, 1, sum)) # normalize
A <- K/(W %*% t(W))
# output <- svd(A, 20)
# dm <- output$u %*% diag(output$d[1:20])
# scatterplot3d(dm[,2:4])

# make A sparse
delta <- 10e-5
ind <- which(A > delta, arr.ind = TRUE)
Asp <- Matrix::sparseMatrix(i = ind[,1], j = ind[,2], x = A[ind], dims=c(n,n))
f <- function(x, A = NULL){ # matrix multiplication for ARPACK
  as.matrix(A %*% x)
}
cat('Performing eigendecomposition\n') # eigendecomposition
neigen <- 20  
neff <- neigen+1
# eigendecomposition using ARPACK
decomp = igraph::arpack(f,extra=Asp,sym=TRUE,
                        options=list(which='LA',nev=neff,n=n,ncv=max(min(c(n,4*neff)))))
psi = decomp$vectors/(decomp$vectors[,1]%*%matrix(1,1,neff))#right ev
phi = decomp$vectors * (decomp$vectors[,1]%*%matrix(1,1,neff))#left ev
eigenvals = decomp$values #eigenvalues

t <- 1
lambda=eigenvals[-1]^t
lambda=rep(1,n)%*%t(lambda)
X = psi[,2:(neigen+1)]*lambda[,1:neigen] #diffusion coords. X

scatterplot3d(X[,1:3])

# concatenate up to time i  -----------------------------------------------
for(i in 1:9){
  name <- paste0("X", i)
  assign(name, array(NA, c(n, 2, i)))
}
for(i in 1:9){
  for(j in 1:i){
    temp_y <- eval(parse(text = paste0("y",j)))
    eval(parse(text = paste0("X", i, "[,1,", j, "] <- temp_y"))) 
    temp_x <- eval(parse(text = paste0("x",j)))
    eval(parse(text = paste0("X", i, "[,2,", j, "] <- temp_x"))) 
  }
}
# # sanity check 
# c5 <- diffuse(dist(X[,,5]))
# plot(c5)

# try new diffusion kernal on concatenated --------------------------------
D <- as.matrix(dist(X9))
n <- dim(D)[1]
K <- exp(-D^2)/eps.val # kernal
v <- sqrt(apply(K, 1, sum)) # normalize
A <- K/(v%*%t(v))

