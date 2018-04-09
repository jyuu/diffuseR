
# attempt at diffussion mapping over time

# load --------------------------------------------------------------------
#library(R.matlab)


# pics --------------------------------------------------------------------
#file <- "/Users/jcahoon/Downloads/Sec_6.1/Patch/"
#ex <- "npatchJ August.mat"
#readMat(paste0(file, ex))
file <- "/Users/jcahoon/workspace/diffuser/raw/"
ex <- "oct_chg.dat"
ex2 <- "nov.dat"
X <- read.table(paste0(file, ex), sep=",")
Y <- read.table(paste0(file, ex2), sep=",")

# need to convert to mat 6 ------------------------------------------------
# done 


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
D <- as.matrix(dist(X))
K <- exp(-D^2)/ep.octc # kernal
v <- sqrt(apply(K, 1, sum)) # normalize
A <- K/(v%*%t(v))
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