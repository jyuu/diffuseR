rm(list=ls())


# load --------------------------------------------------------------------
library(diffusionMap)
library(MASS)
library(matrixcalc)
library(expm)

# create three clusters ---------------------------------------------------
n <- 20
r <- sqrt(runif(n, 0, .3))
theta <- 2*pi*runif(n)
# convert to cartesian
x1 <- r*cos(theta)
y1 <- r*sin(theta)
x2 <- x1 + 1
y2 <- y1 - 2.5
x3 <- x2 - .5
y3 <- y2 + 1.5
# combine 
x <- c(x1, x3, x2)
y <- c(y1, y3, y2)
# viz
plot(x,y) 
# use diffusion maps
d <- dist(cbind(x,y))
dm <- diffuse(d, neigen=2)
output <- as.matrix(dist(dm$X))
heatmap(output, Rowv=NA, Colv=NA, xlab=NULL, ylab=NULL)
# t = 10 
dm <- diffuse(d, t=10, neigen=2)
output <- as.matrix(dist(dm$X))
heatmap(output, Rowv=NA, Colv=NA, xlab=NULL, ylab=NULL)
# t = 50
dm <- diffuse(d, t=50, neigen=2)
output <- as.matrix(dist(dm$X))
heatmap(output, Rowv=NA, Colv=NA, xlab=NULL, ylab=NULL)


# manual ------------------------------------------------------------------
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
X <- cbind(x1,y1)
dists <- as.matrix(dist(X))
eps <- epsilonCompute(dists)
W <- exp(-dists^2/eps)
D <- diag(apply(W,1,sum))
K <- ginv(D) %*% W
sub <- diag(diag(D)^(-.5))
#sub <- expm(-.5*logm(D))
S <-  sub %*% W %*% sub
output <- svd(S)
V <- sub %*% output$u
V <- V/V[,1]
scatterplot3d(V[,2:4])
