## Andres Ordoñez & Daniel Rafael Miranda-Esquivel (Laboratorio de Sistemática y Biogeografía, Universidad Industrial de Santander, Bucaramanga, Colombia)
## Bachelor Thesis: Composition and morphometric analysis of Canthon genus in Santander and south of Bolívar localities 
## Part: Disparity between localities
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
library(vegan)
library(philentropy)
library(colorspace)
library(iNEXT)
library(ggplot2)

# 1. Read data ####

######## read as tps and align
head <- readland.tps("./data/Aligned/Head/head_alg.tps"); head_alg <- gpagen(head)
metaesterno <- readland.tps("./data/Aligned/Metaesternum/metaesternum_alg.tps"); metaesterno_alg <- gpagen(metaesterno)
pronoto <- readland.tps("./data/Aligned/Pronotum/pronotum_alg.tps"); pronoto_alg <- gpagen(pronoto)
protibia <- readland.tps("./data/Aligned/Protibia/protibia_alg.tps"); protibia_alg <- gpagen(protibia)

# create a list that contains tps
list_structures <- list(head = head_alg, metaesterno = metaesterno_alg,
                        pronoto = pronoto_alg, protibia = protibia_alg)

######## read information matrix
info_head <- read.csv("./Head_results/info_head.csv", header = TRUE)
info_metaesterno <- read.csv("./Metaesterno_results/info_metaesternum.csv", header = TRUE)
info_pronoto <- read.csv("./Pronoto_results/info_pronoto.csv", header = TRUE)
info_protibia <- read.csv("./Protibia_results/info_protibia.csv", header = TRUE)

# create a list that contains info matrices
list_information <- list(head = info_head, metaesterno = info_metaesterno,
                         pronoto = info_pronoto, protibia = info_protibia)


# 2. Create directories  ######
dir.create("./Head_results/Geometrical_Morpho_results/Disparity", showWarnings = FALSE)
dir.create("./Metaesterno_results/Geometrical_Morpho_results/Disparity", showWarnings = FALSE)
dir.create("./Pronoto_results/Geometrical_Morpho_results/Disparity", showWarnings = FALSE)
dir.create("./Protibia_results/Geometrical_Morpho_results/Disparity", showWarnings = FALSE)

# list with the names of results directories
resultados <- c("./Head_results/Geometrical_Morpho_results/Disparity/", "./Metaesterno_results/Geometrical_Morpho_results/Disparity/", "./Pronoto_results/Geometrical_Morpho_results/Disparity/", "./Protibia_results/Geometrical_Morpho_results/Disparity/")

# 3. Calculate partial disparity ####

list_disparity <- list(head = NA, metaesterno = NA,
                         pronoto = NA, protibia = NA)

for (i in 1:4) {
  
  tmp_DF <- geomorph.data.frame(coord = list_structures[i][[1]]$coords, species = list_information[i][[1]]$Morfotipo, 
                                localidad = list_information[i][[1]]$Localidad, municipio = list_information[i][[1]]$Municipio,
                                altura500 = list_information[i][[1]]$X500mts, altura1000 = list_information[i][[1]]$X1000mts)
  
  MD_localidad <- morphol.disparity(coord ~ 1, groups= ~ localidad, partial = T,
                                    data = tmp_DF, iter = 500, print.progress = FALSE)
  
  MD_municipio <- morphol.disparity(coord ~ 1, groups= ~ municipio, partial = T,
                                    data = tmp_DF, iter = 500, print.progress = FALSE)
  
  MD_500mts <- morphol.disparity(coord ~ 1, groups= ~ altura500, partial = T,
                                    data = tmp_DF, iter = 500, print.progress = FALSE)
  
  MD_1000mts <- morphol.disparity(coord ~ 1, groups= ~ altura1000, partial = T,
                                    data = tmp_DF, iter = 500, print.progress = FALSE)
  
  list_disparity[i][[1]] <-  list(locality = MD_localidad, munic = MD_municipio,
                                  alt_500 = MD_500mts, alt_1000 = MD_1000mts)
   
}  

# 3.1 plots and csv from disparity ####

list_ProVar <- list(head = list(localidad = NA, municipio = NA,
                                altura500 = NA, altura1000 = NA), 
                    metaesterno = list(localidad = NA, municipio = NA,
                                       altura500 = NA, altura1000 = NA),
                    pronoto = list(localidad = NA, municipio = NA,
                                      altura500 = NA, altura1000 = NA), 
                    protibia = list(localidad = NA, municipio = NA,
                                    altura500 = NA, altura1000 = NA)) #save procrustes Var for final figures

for (i in 1:4) { #for each structure
  
  for (j in 1:4) { #for each category in disparity comparation (localities, munic, 500mts and 1000mts)
    
    list_ProVar[i][[1]][j][[1]] <- list_disparity[i][[1]][j][[1]]$Procrustes.var #save procrustes Var for final figures
    
    for (k in 1:3) { #for each matrix calculated (var.Proc, Disparity distance and Significance differnts)
      
      tmp_name <- paste0(resultados[i], names(list_structures)[i], "_", names(list_disparity[i][[1]])[j]) #temporal name for outputs
      
      pdf(file = paste0(tmp_name, ".pdf")) #save pdf files with all images
        #barplot
        tmp_data <- as.data.frame(list_disparity[i][[1]][j][[1]]$Procrustes.var)
        names(tmp_data) <- "disparity"
        
        p <- ggplot(data=tmp_data, aes(y=disparity, x=row.names(tmp_data))) +
          geom_bar(stat="identity", position=position_dodge()) +
          labs(x = "", y = "Procrustes Var", fill = "") +
          coord_flip() +
          theme(legend.position = "right",
                axis.text=element_text(size=14),
                axis.title=element_text(size=18),
                panel.background = element_rect("white"),
                panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                                colour = "gray90"), 
                panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                colour = "gray90"))
        
        plot(p)
        
        #cluster
        dist.tmp <- as.dist(list_disparity[i][[1]][j][[1]]$PV.dist)
        
        clust.tmp <- hclust(dist.tmp, method = "average")
        
        plot(clust.tmp ,horiz=F,axes=T, ylab= "Distancia", cex.axis=1.3,cex.lab=1.3)
        
      dev.off()
      
      #csv
      tmp_name <- paste0(resultados[i], names(list_structures)[i], "_", names(list_disparity[i][[1]])[j], "_", names(list_disparity[i][[1]][j][[1]])[k]) #temporal name for outputs
      write.csv(x = list_disparity[i][[1]][j][[1]][k], file = paste0(tmp_name, ".csv"))
      
    }
    
  }
  
}

# 4. final figures ####

dir.create("./Diversity", showWarnings = FALSE) #create a directory
dir.create("./Diversity/Disparity", showWarnings = FALSE) #create a directory

list_ProVar

for (i in 1:4) {
  
  tmp_list <- c(list_ProVar$head[i][[1]], list_ProVar$metaesterno[i][[1]],
              list_ProVar$pronoto[i][[1]], list_ProVar$protibia[i][[1]])
  
  tmp_length <- length(list_ProVar$head[i][[1]])
  
  tmp_DF <- data.frame(cate = names(tmp_list), proVar = tmp_list, str = c(rep("Head", tmp_length), rep("Metaesternum", tmp_length), 
                                                                          rep("Pronotum", tmp_length), rep("Protibia", tmp_length)))
  
  
  bar_plot <- ggplot(data = tmp_DF, aes(x = cate, y = proVar, fill = str)) +
    geom_bar(stat="identity", position=position_dodge()) +
    labs(x = "", y = "Disparidad", fill = "") +
    coord_flip() + 
    theme(legend.position = "right",
          axis.text=element_text(size=14),
          axis.title=element_text(size=18),
          panel.background = element_rect("white"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "gray90"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "gray90"))
    
  
  png(paste0("./Diversity/Disparity/", names(list_ProVar$head[i]), "_BarDisparity.png"), width = 420*3, height = 290*3, res = 200)
    plot(bar_plot)
  dev.off()
}
