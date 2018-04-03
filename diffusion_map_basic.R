#' diffusion map
#'
#' @param X input matrix to reduce dimension for
#' @param d number of dimension for matrix with reduced columns
#' @param h kernel bandwidth
#' @export

diffusion_map_basic <- function(X,d,h){
  
  # 'dist' computes the distances between the rows of a data matrix.
  dist_L <- dist(X)
 
  kern_mat <- dnorm(dist_L,sd=h)
  
  full_KM <- as.matrix(kern_mat)
  
  norm_KM <- t(apply(full_KM,1,norm<-function(x){return (x/sum(x))}))
  
  eig <- eigen(norm_KM, symmetric = TRUE)
  
  df_map <- matrix(0, ncol=d, nrow=nrow(X) ) 
  for(i in 1:d){
    df_map[,i] = eig$values[i]*eig$vectors[,i]
  }
  
  return(df_map)
}