rm(list=ls())

# Notes -------------------------------------------------------------------
# 290 and 315 could be possible frames to stop at
# the frames start at 0023
path <- "/home/joyce/Downloads/frames_small/output_"
frames <- seq(23, 426, by = 1)
library(stringr)
frames <- str_pad(frames, 4, pad = 0)
# read in the first 10 images
library(png)
x <- readPNG(paste0(path, frames[133], ".png"))
y <- readPNG(paste0(path, frames[134], ".png"))

# functions ---------------------------------------------------------------
# # make the image into a matrix with each channel
# X <- array(NA, dim= c(dim(x)[1]*dim(x)[2], 3))
# for(i in 1:3){
#   X[,i] <- c(x[,,i])
# }
# # try reconstructing 
# joyce <- X[,1]
# ross <- matrix(joyce, dim(x)[1])
# image(ross, axes=FALSE) 
# # this works! so functionalize 
create_tensor <- function(image){
  I <- array(NA, dim= c(dim(image)[1]*dim(image)[2], 3))
  for(i in 1:3){
    I[,i] <- c(image[,,i])
  }
  return(I)
}
# # try it
# y <- create_tensor(x)
X <- create_tensor(x)
Y <- create_tensor(y)
# sanity check 
image(matrix(X[,1], 100, 100), axes=FALSE)
image(matrix(Y[,1], 100, 100), axes=FALSE)

# set params --------------------------------------------------------------
num_eig <- 20
t <- 1

# Feed into our diff change -----------------------------------------------
D.X <- as.matrix(dist(X))
d.X <- diffuse(D.X, neigen = num_eig, t = t)
D.Y <- as.matrix(dist(Y))
d.Y <- diffuse(D.Y, neigen = num_eig, t = t)
# rm heavy
# rm(X)
# rm(D.X)
# rm(Y)
# rm(D.Y)
# now apply the paper stuff
VX <- d.X$X
EX <- diag(d.X$eigenvals)
VY <- d.Y$X
EY <- diag(d.Y$eigenvals)
# sanity check 
image(matrix(VX[,1], 100, 100), axes=FALSE)
image(matrix(VY[,2], 100, 100), axes=FALSE)

# omap function -----------------------------------------------------------
omap <- function(DMs, Vs, Ve){
  O <- t(Vs) %*% Ve 
  DM_e <- DMs %*% O
  return(DM_e)
}
# their stuff
if (t < Inf){
  DMX <- VX %*% (EX^t);
  DMY_I <- VY %*% (EY^t)
  DMY <- omap(DMY_I,VY,VX);
}
# sanity check
image(matrix(DMY[,1], 100, 100), axes=FALSE)
image(matrix(DMX[,20], 100, 100), axes=FALSE)


# chg map -----------------------------------------------------------------
chg_map = sqrt(apply((DMX-DMY)^2, 1, sum));
image(matrix(chg_map, 100, 100), axes=FALSE)

# Try visualizing mid way -------------------------------------------------
x <- DMX[,1]
y <- DMX[,2]
k <- 20
x <- VX[,k]
y <- VY[,k]
image(matrix(x, 100, 100), axes=FALSE)
image(matrix(y, 100, 100), axes=FALSE)


