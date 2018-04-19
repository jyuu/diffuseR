rm(list=ls())
# load --------------------------------------------------------------------
#library(R.matlab)
#library(reshape2)
#library(ggplot2)
library(diffusionMap)
library(png)
#library(spatstat)

clock <- readPNG("~/Desktop/clock_small.png")
clock <- as.matrix(clock[,,1])
pixel <- im(clock)
image(pixel)
angles <- sample(1:360, 20)
saveRDS(angles, "results/dm_clock_angles.rds")
data <- array(1, c(115, 115, length(angles)))
for(i in 1:length(angles)){
  y <- rotate(pixel, angle=angles[i])
  temp <- y$v 
  temp[is.na(temp)] <- 1
  data[1:(dim(temp)[1]), 1:(dim(temp)[2]), i] <- temp
}
# convert to tensor -------------------------------------------------------
X <- array(NA, c(dim(data)[1]*dim(data)[2], dim(data)[3]))
for(i in 1:dim(X)[2]){
  X[,i] <- c(data[,,i])
}
#sanity check
a_clock <- matrix(X[,15], dim(data)[1])
image(a_clock)
dist_X <- as.matrix(dist(X))
#saveRDS(dist_X, "results/dm_clock_dst.rds")
rm(list=ls())
dist_X <- readRDS("results/dm_clock_dst.rds")
dm <- diffuse(dist_X, neigen = 6, t= 0)
# save 
saveRDS(dm, "results/dm_clock.rds")

# # data --------------------------------------------------------------------
# joyce <- readMat("~/Downloads/face.mat")
# data <- joyce$Y
# for(j in 1:dim(data)[3]){
#   a_face <- data[,,j]
#   test <- melt(a_face)
#   image(a_face)
#   # ggplot(test, aes(x = Var2, y = Var1, fill = factor(value))) + 
#   #   geom_raster() + theme(legend.position="none") + 
#   #   scale_fill_grey(start = 0, end = .9)
#   # Sys.sleep(1)
# }
# # but i dont know the original angles of rotation

# # reverse images practice -------------------------------------------------
# m1 <- matrix(c(1,0,0,0,1,0,1,0,0), 3)
# image(m1) # want it facing up! 
# m2 <- apply(m1, 2, rev) 
# image(t(m2)) # now it is correct.. 
# # now try it on our data
# a_face <- data[,,10]
# temp <- apply(a_face, 2, rev)
# image(t(temp))

# # make data into tensor ---------------------------------------------------
# X <- array(NA, c(dim(data)[1]*dim(data)[2], dim(data)[3]))
# for(i in 1:dim(X)[2]){
#   X[,i] <- c(data[,,i])
# }
# # sanity check 
# # a_face <- matrix(X[,1], dim(data)[1]) 
# # image(a_face)
# 
# # apply diffusion map to it -----------------------------------------------
# dist_X <- dist(X)
# dm <- diffuse(dist_X)
# 
# # viz ---------------------------------------------------------------------
# a_face <- data[,,1]
# image(a_face)
# # actual
# second_face <- data[,,2]
# image(second_face)
# # from algo 
# rot <- matrix(dm$X[,2], dim(data)[1])
# test <- a_face %*% t(rot)
# image(test)
