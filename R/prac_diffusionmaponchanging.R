
# attempt at diffussion mapping over time

# load --------------------------------------------------------------------
#library(R.matlab)

# pics --------------------------------------------------------------------
#file <- "/Users/jcahoon/Downloads/Sec_6.1/Patch/"
#ex <- "npatchJ August.mat"
#readMat(paste0(file, ex))
#file <- "/Users/jcahoon/workspace/diffuser/raw/"
file <- "/home/joyce/workspace/diffuser/raw/"
ex <- "oct_chg.dat"
ex2 <- "nov.dat"
X <- read.table(paste0(file, ex), sep=",")
Y <- read.table(paste0(file, ex2), sep=",")

# eps values --------------------------------------------------------------
ep.aug = 0.233972639306330;
ep.sep = 0.242537079542777;
ep.oct = 0.151503508527451;
ep.octc = 0.138122144033416; 
ep.nov = 0.161087707605546;
# number of eigenvectors to keep 
num_eig = 20;
# diffusion time 
t = 1;

# try it with built in R package ------------------------------------------
library(diffusionMap)
D.X <- as.matrix(dist(X))
d.X <- diffuse(D.X, eps.val= ep.octc, neigen = num_eig, t = t)
rm(D.X)
rm(X)
# now apply to second one
D.Y <- as.matrix(dist(Y))
d.Y <- diffuse(D.Y, eps.val= ep.nov, neigen = num_eig, t = t)
# now apply the paper stuff
VX <- d.X$X
EX <- diag(d.X$eigenvals)
VY <- d.Y$X
EY <- diag(d.Y$eigenvals)

# omap function -----------------------------------------------------------
omap <- function(DMs, Vs, Ve){
  O <- t(Vs) %*% Ve 
  DM_e <- DMs %*% O
  return(DM_e)
}
# their stuff
t <- Inf
if (t < Inf){
  DMX <- VX %*% (EX^t);
  DMY_I <- VY %*% (EY^t)
  DMY <- omap(DMY_I,VY,VX);
}

# Try visualizing mid way -------------------------------------------------
x <- DMX[,1]
y <- DMX[,2]
# x <- VX[,1]
# y <- VY[,2]
image(matrix(x, 100, 100), axes=FALSE)
image(matrix(y, 100, 100), axes=FALSE)

# change map 
if (t < Inf){
  chg_map = sqrt(apply((DMX-DMY)^2, 1, sum));
} else {
  # joyce <- (VX[,1]-VY[,1])^2
  # joyce <- VX[,1] * VY[,1]*norm(as.matrix(VX[,1]-VY[,1]),type="2")^2
  # joyce <-(VX[,1]-VY[,1])^2 + VX[,1] * VY[,1]*norm(as.matrix(VX[,1]-VY[,1]),type="2")^2
  chg_map <- sqrt((VX[,1]-VY[,1])^2 + VX[,1] * VY[,1]*norm(as.matrix(VX[,1]-VY[,1]),type="2")^2);
}

chg_map = matrix(chg_map,100,100);

# display
# Display exponent for visualization
if (t < Inf){
  alpha = (1/2000)*t + 0.5*(1-1/1000);
}else{
  alpha = 1;
}

# Plot
# imagesc(chg_map.^alpha);
# axis('image');
# after!! 
image(chg_map^alpha, axes=FALSE) # NO CHANGE DETECTED!!
image(chg_map^.1, axes=FALSE) # NO CHANGE DETECTED!!

# lets compare to before 
octc <- X[,1] # one spectrum
oct.pic <- matrix(octc, 100, 100)
image(oct.pic, axes=FALSE)
nov <- Y[,1] # one spectrum
nov.pic <- matrix(nov, 100, 100)
image(nov.pic, axes=FALSE)


# make the symmetric diffusion kernals 
D.X <- as.matrix(dist(X))
K.X <- exp((-D.X^2)/ep.octc^2) # kernal
v.X <- sqrt(apply(K.X, 1, sum)) # normalize
A.X <- K.X/(v.X%*%t(v.X))
# make the same symmetric diffusion kernals for Y 
D.Y <- as.matrix(dist(Y))
K.Y <- exp((-D.Y^2)/ep.nov^2) # kernal
v.Y <- sqrt(apply(K.Y, 1, sum)) # normalize
A.Y <- K.Y/(v.Y%*%t(v.Y))
# set epsilon to make into sparse matrix
delta <- 10^{-5}
n <- dim(D)[1]
ind <- which(A > delta, arr.ind = TRUE)
Asp <- Matrix::sparseMatrix(i = ind[,1], j = ind[,2], x = A[ind], dims=c(n,n))
# sub fun 
f <- function(x, A = NULL){ # matrix multiplication for ARPACK
  as.matrix(A %*% x)
}
decomp = igraph::arpack(f,extra=Asp,sym=TRUE,
                        options=list(which='LA',nev=51,n=dim(A)[1],ncv=200))
psi = decomp$vectors/(decomp$vectors[,1]%*%matrix(1,1,51))#right ev
phi = decomp$vectors * (decomp$vectors[,1]%*%matrix(1,1,51))#left ev
eigenvals = decomp$values #eigenvalues
# diffusion map for t < Inf 
DMX <- VX %*% diag(eigenvals_X) 


DMY <- O