rm(list=ls())

# load --------------------------------------------------------------------
library(dplyr)
library(raster)
library(diffusionMap)

# # data --------------------------------------------------------------------
# filepath <- "/media/joyce/fractals/raw-data/raw-BCDR/BCDR-D02_dataset/"
# labels <- readr::read_csv(paste0(filepath, "bcdr_d02_features.csv"))
# 
# # outlines ----------------------------------------------------------------
# outlines <- readr::read_csv(paste0(filepath, "bcdr_d02_outlines.csv"))
# 
# # adjust to include filepaths to images -----------------------------------
# mapper <- function(id){
#   if(id == 1){
#     type <- "RCC"
#   } else if(id == 2){
#     type <- "LCC"
#   } else if(id == 3){
#     type <- "RO"
#   } else{
#     type <- "LO"
#   }
#   return(type)
# }
# labels <- labels %>%
#   rowwise() %>%
#   mutate(image_clean = mapper(image_view))
# with_names <- labels %>%
#   mutate(folder_names = paste0(filepath, "patient_", patient_id, "/study_", study_id,
#                                "/img_", patient_id, "_", study_id, "_", series, "_", 
#                                image_clean, ".tif"))
# saveRDS(with_names, "/home/joyce/workspace/diffuser/rdata/features.rds")
# saveRDS(outlines, "/home/joyce/workspace/diffuser/rdata/outlines.rds")

# read in some ------------------------------------------------------------
features <- readRDS("/home/joyce/workspace/diffuser/rdata/features.rds")
outlines <- readRDS("/home/joyce/workspace/diffuser/rdata/outlines.rds")

# # lazy way
# i <- 2
# ex1 <- raster(features$folder_names[i])
# plot(ex1)
# outx1 <- as.numeric(unlist(strsplit(outlines$lw_x_points[i], "\\s+")))
# outy1 <- as.numeric(unlist(strsplit(outlines$lw_y_points[i], "\\s+")))
# lines(outx1, outy1, new =T)

# use diffusionmap --------------------------------------------------------
i <- 2
ex1 <- raster(features$folder_names[i])
ex <- as.matrix(ex1)
ex <- ex[1:min(dim(ex)), 1:min(dim(ex))]
dmap <- diffuse(ex, eps.val=.1)
# try diffusion kmeans 
# Diffusion K-means is a special form of spectral clustering. It is a unique algorithm because the eigenvectors of the symmetric Laplacian are weighted in such a way to guarantee that Euclidean distance in diffusion space will be approximately equal to the diffusion distance between objects. Clustering by Euclidean distance in diffusion space exploits this fact.
output <- diffusionKmeans(dmap, K = 3, Niter = 20)

dim(dmap$X)
dmap$eigenvals

test <- dmap$X %*% diag(dmap$eigenvals) %*% t(dmap$X)


# use ggplot2 -------------------------------------------------------------
imgDm <- dim(ex1)
ex <- ex1
img <- data.frame(
  x = rep(1:imgDm[2], each = imgDm[1]), 
  y = rep(imgDm[1]:1, imgDm[2]), 
  z = as.vector(ex[,,1])
)
plotter <- function(){
  theme(
    panel.background = element_rect(
      size = 3,
      colour = "black",
      fill = "white"),
    axis.ticks = element_line(
      size = 2),
    panel.grid.major = element_line(
      colour = "gray80",
      linetype = "dotted"),
    panel.grid.minor = element_line(
      colour = "gray90",
      linetype = "dashed"),
    axis.title.x = element_text(
      size = rel(1.2),
      face = "bold"),
    axis.title.y = element_text(
      size = rel(1.2),
      face = "bold"),
    plot.title = element_text(
      size = 20,
      face = "bold",
      vjust = 1.5)
  )
}
ggplot(data = img, aes(x=x, y=y)) +
  geom_point(aes(img[c("z")]))

# locate microcalcification mass ------------------------------------------


# notes -------------------------------------------------------------------
# 
