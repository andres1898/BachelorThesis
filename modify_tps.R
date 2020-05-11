## Andres Ordoñez & Daniel Rafael Miranda-Esquivel (Laboratorio de Sistemática y Biogeografía, Universidad Industrial de Santander, Bucaramanga, Colombia)
## Bachelor Thesis: Composition and morphometric analysis of Canthon genus in Santander and south of Bolívar localities 
## Part: Make scaled and align tps files
## R version 3.6.3
## May 2020

# 1. Load packeges ####

library(geomorph)
library(shapes)
library(Morpho)
library(stringr)
library(NbClust)
library(cluster)
library(colorspace)

# 0. Function to scale the the forma based in coordinates #####
Manual_scale <- function(tps, landmark_1, landmark_2, scale_size) {
  
  for (i in 1:dim(tps)[3]){
    scb <- tps[c(landmark_1, landmark_2),,i] #extract references coordinates
    scale <- sqrt((scb[1,1]-scb[2,1])^2+(scb[1,2]-scb[2,2])^2)/scale_size  #Cartesian plane displacement
    tps[,,i] <- tps[,,i]/scale
  }
  
  tps <- tps[-c(landmark_1, landmark_2),,] # remove the landmarks
  
  return(tps)
  
}


# 2. Create directories for each step ######
tmp <- c("Head", "Metaesternum", "Pronotum", "Protibia")
for (i in 1:length(tmp)) {
  dir.create("./data/Scaled/", showWarnings = FALSE) # Scaled coordinates
  dir.create("./data/Aligned/", showWarnings = FALSE) # Aligned coordinates
  dir.create(paste0("./data/Scaled/", tmp[i]), showWarnings = FALSE) # Scaled coordinates
  dir.create(paste0("./data/Aligned/", tmp[i]), showWarnings = FALSE) # Aligned coordinates
}


# 3. Head ####
# read data
files_names <- list.files(path = "./data/Raw/head/") #list of files names
files_names <- files_names[-8]

# configurations have diferents landmarks order
tmp1.2 <- files_names[c(2,3)]  
tmp11.12 <- files_names[-c(2,3)]

# read as tps
tps1.2 <- lapply(paste0("./data/Raw/head/", tmp1.2), specID = "imageID", readland.tps)
tps11.12 <- lapply(paste0("./data/Raw/head/", tmp11.12), specID = "imageID", readland.tps)

# scale data
scale1.2 <- lapply(tps1.2, Manual_scale, 1, 2, 0.4)
scale11.12 <- lapply(tps11.12, Manual_scale, 11, 12, 0.4)

scaled_all <- c(scale1.2, scale11.12)

# align data
aligned <- lapply(scaled_all, gpagen)

# save data
for (i in 1:length(files_names)) {
  
  writeland.tps(aligned[i][[1]][[1]], file = paste0("./data/Aligned/Head/", files_names[i]))
  writeland.tps(scaled_all[i][[1]], file = paste0("./data/Scaled/Head/", files_names[i]))
  
}


# 4. Metaesternum ####
# read data
files_names <- list.files(path = "./data/Raw/metasternum/") #list of files names
files_names <- files_names[-8]

# read as tps
tps <- lapply(paste0("./data/Raw/metasternum/", files_names), specID = "imageID", readland.tps)

# scale data
scale <- lapply(tps, Manual_scale, 10, 11, 0.4)

# align data
aligned <- lapply(scale, gpagen)

# save data
for (i in 1:length(files_names)) {
  
  writeland.tps(aligned[i][[1]][[1]], file = paste0("./data/Aligned/Metaesternum/", files_names[i]))
  writeland.tps(scale[i][[1]], file = paste0("./data/Scaled/Metaesternum/", files_names[i]))
  
}

# 5. Pronotum ####
# read data
files_names <- list.files(path = "./data/Raw/pronotum/") #list of files names
files_names <- files_names[-8]

# read as tps
tps <- lapply(paste0("./data/Raw/pronotum/", files_names), readland.tps, specID = "imageID")

# scale data
scale <- lapply(tps, Manual_scale, 1, 2, 0.4)

# align data
aligned <- lapply(scale, gpagen)

# save data
for (i in 1:length(files_names)) {
  
  writeland.tps(aligned[i][[1]][[1]], file = paste0("./data/Aligned/Pronotum/", files_names[i]))
  writeland.tps(scale[i][[1]], file = paste0("./data/Scaled/Pronotum/", files_names[i]))
  
}

# 6. Protibia ####
# read data
files_names <- list.files(path = "./data/Raw/protibia/") #list of files names
files_names <- files_names[-8]

# read as tps
tps <- lapply(paste0("./data/Raw/protibia/", files_names), specID = "imageID", readland.tps)

# scale data
scale <- lapply(tps, Manual_scale, 11, 12, 0.4)

# align data
aligned <- lapply(scale, gpagen)

# save data
for (i in 1:length(files_names)) {
  
  writeland.tps(aligned[i][[1]][[1]], file = paste0("./data/Aligned/Protibia/", files_names[i]))
  writeland.tps(scale[i][[1]], file = paste0("./data/Scaled/Protibia/", files_names[i]))
  
}
