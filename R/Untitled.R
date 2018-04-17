
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