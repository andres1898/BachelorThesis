## Andres Ordoñez & Daniel Rafael Miranda-Esquivel (Laboratorio de Sistemática y Biogeografía, Universidad Industrial de Santander, Bucaramanga, Colombia)
## Bachelor Thesis: Composition and morphometric analysis of Canthon genus in Santander and south of Bolívar localities 
## Part: Make measurements for traditional morphometric analyses
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

# 2. Create directories for Geometrical Morphometric ######
dir.create("./Head_results/Traditional_Morpho_results", showWarnings = FALSE)
dir.create("./Metaesterno_results/Traditional_Morpho_results", showWarnings = FALSE)
dir.create("./Pronoto_results/Traditional_Morpho_results", showWarnings = FALSE)
dir.create("./Protibia_results/Traditional_Morpho_results", showWarnings = FALSE)


# 3. Head ####
# read data
head <- readland.tps("data/Scaled/Head/head_scaled.tps", specID = "ID")
information <- read.csv("Head_results/info_head.csv")

measurements <- matrix(nrow = dim(head)[3], ncol = 5)

colnames(measurements) <- c("ID", "eye-clypeo", "clypeal sutures", "extenal edge eye", "inner edge eye")

# measurements based in landmarks
for (i in 1:dim(head)[3]){
  
  single <- head[,,i]
  
  measurements[i,2] <- sqrt((single[1,1]-single[4,1])^2+(single[1,2]-single[4,2])^2)
  
  measurements[i,3] <- sqrt((single[11,1]-single[16,1])^2+(single[11,2]-single[16,2])^2)
  
  measurements[i,4] <- sqrt((single[5,1]-single[4,1])^2+(single[5,2]-single[4,2])^2)
  
  measurements[i,5] <- sqrt((single[6,1]-single[3,1])^2+(single[6,2]-single[3,2])^2)
  
}

measurements <- as.data.frame(measurements)

measurements[,1] <- attributes(head)[2][[1]][[3]]

write.csv(x = measurements, file = "Head_results/Traditional_Morpho_results/measurements_head.csv", row.names = F)



# 4. Metaesternum ####
# read data
metaesternum <- readland.tps("data/Scaled/Metaesternum/metaesternum_scaled.tps", specID = "ID")
information <- read.csv("Metaesterno_results/info_metaesternum.csv")

measurements <- matrix(nrow = dim(metaesternum)[3], ncol = 5)

colnames(measurements) <- c("ID", "anterior width", "middle width", "posterior width", "long")

for (i in 1:dim(metaesternum)[3]){
  
  single <- metaesternum[,,i]
  
  measurements[i,2] <- sqrt((single[2,1]-single[9,1])^2+(single[2,2]-single[9,2])^2)
  
  measurements[i,3] <- sqrt((single[8,1]-single[3,1])^2+(single[8,2]-single[3,2])^2)
  
  measurements[i,4] <- sqrt((single[7,1]-single[4,1])^2+(single[7,2]-single[4,2])^2)
  
  measurements[i,5] <- sqrt((single[1,1]-single[5,1])^2+(single[1,2]-single[5,2])^2)
  
}

measurements <- as.data.frame(measurements)

measurements[,1] <- attributes(metaesternum)[2][[1]][[3]]

write.csv(x = measurements, file = "Metaesterno_results/Traditional_Morpho_results/measurements_metaesternum.csv", row.names = F)

# 5. Pronotum ####
# read data
Pronotum <- readland.tps("data/Scaled/Pronotum/pronotum_scaled.tps", specID = "ID")
information <- read.csv("Pronoto_results/info_pronoto.csv")

measurements <- matrix(nrow = dim(Pronotum)[3], ncol = 6)

colnames(measurements) <- c("ID", "anterior angles", "posterior angles", "middle length",  "rigth length", "middle width")

for (i in 1:dim(Pronotum)[3]){
  
  single <- Pronotum[,,i]
  
  measurements[i,2] <- sqrt((single[1,1]-single[13,1])^2+(single[1,2]-single[13,2])^2)
  
  measurements[i,3] <- sqrt((single[9,1]-single[5,1])^2+(single[9,2]-single[5,2])^2)
  
  measurements[i,4] <- sqrt((single[16,1]-single[7,1])^2+(single[16,2]-single[7,2])^2)
  
  measurements[i,5] <- sqrt((single[1,1]-single[5,1])^2+(single[1,2]-single[5,2])^2)
  
  measurements[i,6] <- sqrt((single[11,1]-single[3,1])^2+(single[11,2]-single[3,2])^2)
  
}

measurements <- as.data.frame(measurements)

measurements[,1] <- attributes(Pronotum)[2][[1]][[3]]

write.csv(x = measurements, file = "Pronoto_results/Traditional_Morpho_results/measurements_pronotum.csv", row.names = F)

# 6. Protibia ####
# read data
protibia <- readland.tps("data/Scaled/Protibia/protibia_scaled.tps", specID = "ID")
information <- read.csv("Protibia_results/info_protibia.csv")

measurements <- matrix(nrow = dim(protibia)[3], ncol = 4)

colnames(measurements) <- c("ID", "total length", "distance between teeth",  "anterior width")

for (i in 1:dim(protibia)[3]){
  
  single <- protibia[,,i]
  
  measurements[i,2] <- sqrt((single[1,1]-single[9,1])^2+(single[1,2]-single[9,2])^2)
  
  measurements[i,3] <- sqrt((single[1,1]-single[8,1])^2+(single[1,2]-single[8,2])^2)
  
  measurements[i,4] <- sqrt((single[1,1]-single[10,1])^2+(single[1,2]-single[10,2])^2)
  
}

measurements <- as.data.frame(measurements)

measurements[,1] <- attributes(protibia)[2][[1]][[3]]

write.csv(x = measurements, file = "Protibia_results/Traditional_Morpho_results/measurements_protibia.csv", row.names = F)


