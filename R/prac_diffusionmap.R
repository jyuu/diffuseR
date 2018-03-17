
# Load library ------------------------------------------------------------
library(diffusionMap)
library(scatterplot3d)

# Ex 1 --------------------------------------------------------------------
t <- seq(-pi, pi, .01)
x <- cbind(cos(t), sin(t))
y <- cos(3*t) + rnorm(length(t), 0, .1)
tcol <- topo.colors(32)
colvec <- floor((y-min(y))/(max(y)-min(y))*32)
colvec[colvec==0] <- 1
scatterplot3d(x[,1], x[,2], y, color=tcol[colvec], pch=20, 
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


# General diffuse ---------------------------------------------------------
# for data parameterization, that exploits the natural geometry of a data 
# set, diffusion map uses local interactions between data points, propagated
# to larger scales, to construct a global representation of the data.
# computation requires singular value decomposition of the normalized
# graph laplacian. the operation is optimized for speed by exploiting the
# sparseness of the graph laplacian and by using ARPACK for fast matrix
# decomposition, increasing the sparseness parameter, delta, will speed
# up the algo
n <- 2000 
t <- runif(n)^.7*10 
al <- .15
bet <- .5
# ex with noisy spiral 
x1 <- bet*exp(al*t)*cos(t)+rnorm(n,0,.1)
y1 <- bet*exp(al*t)*sin(t)+rnorm(n,0,.1)
plot(x1, y1, pch=20, main="Noisy spiral")
D <- dist(cbind(x1, y1))
dmap <- diffuse(D, neigen=10)
par(mfrow=c(2,1))
plot(t, dmap$X[,1], pch=20, axes=FALSE, xlab="spiral parameter",
     ylab="1st diffusion coefficient")
box()
plot(1:10, dmap$eigenmult, typ="h", xlab="diffusion map dimension", 
     ylab="eigen-multipliers")


# try it with annulus -----------------------------------------------------
data("annulus")
plot(annulus, main="Annulus Data", pch=20, cex=.7)
D <- dist(annulus)
dmap <- diffuse(D, eps.val=.1)
print(dmap)
plot(dmap)
