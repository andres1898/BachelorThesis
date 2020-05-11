## Andres Ordoñez & Daniel Rafael Miranda-Esquivel (Laboratorio de Sistemática y Biogeografía, Universidad Industrial de Santander, Bucaramanga, Colombia)
## Bachelor Thesis: Composition and morphometric analysis of Canthon genus in Santander and south of Bolívar localities 
## Part: Geometrical morphometric analyses
## R version 3.6.3
## May 2020

# 1. Load packeges ####

library(geomorph)
library(shapes)
library(ggplot2)
library(Morpho)
library(ellipse)
library(stringr)
library(factoextra)
library(FactoMineR)
library(fpc)
library(dendextend)
library(NbClust)
library(cluster)
library(colorspace)
library(phangorn)


# 2. Create directories for Geometrical Morphometric ######
dir.create("./Head_results/Geometrical_Morpho_results", showWarnings = FALSE)
dir.create("./Metaesterno_results/Geometrical_Morpho_results", showWarnings = FALSE)
dir.create("./Pronoto_results/Geometrical_Morpho_results", showWarnings = FALSE)
dir.create("./Protibia_results/Geometrical_Morpho_results", showWarnings = FALSE)

# 3. Read data ####

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

# list with the names of results directories
resultados <- c("./Head_results/Geometrical_Morpho_results/", "./Metaesterno_results/Geometrical_Morpho_results/", "./Pronoto_results/Geometrical_Morpho_results/", "./Protibia_results/Geometrical_Morpho_results/")

# 4. Generalized Procrustes Analyses ####
for (i in 1:4) { #for each strucutre into the list. Saves the results, doesn't plot it
  png(filename = paste0(resultados[i], names(list_structures)[i], "_GPA.png"), width = 320*5, height = 241*5, res = 300)
    plotAllSpecimens(list_structures[i][[1]]$coords) #ver el gpa y todos align
  dev.off()
}

# 5.1 Principal Component Analyses with deformation plates ####

#create an empty list to fill with the result behind
pca_structures <- list(head = NA, metaesterno = NA, pronoto = NA, protibia = NA)

for (i in 1:4) { #for each strucutre into the list. Saves the results, doesn't plot it
  png(filename = paste0(resultados[i], names(list_structures)[i], "_PCA_deformations.png"), width = 320*10, height = 241*10, res = 300)
  pca_structures[i][[1]] <- plotTangentSpace(list_structures[i][[1]]$coords, warpgrids = TRUE, legend = TRUE, groups = list_information[i][[1]][,"Morfotipo"])
  dev.off()
}

# 5.2 Principal Component Analyses with ellipse ####

for (i in 1:4) { #for each strucutre into the list. Saves the results, doesn't plot it
  
  pca_scores <- as.data.frame(pca_structures[i][[1]]$pc.scores) #individual pca scores as dataframe

  pdf(file = paste0(resultados[i], names(list_structures)[i], "_PCA_ellipses.pdf"))
  
  for (j in c(2,3,5)) {  #corresponde to category in info matrices
    p <- ggplot(pca_scores, aes(PC1, PC2, color = list_information[i][[1]][,j])) +
      geom_point() +
      stat_ellipse(type = "t") +
      xlab(paste0("PC1 = ", 
                  round(pca_structures[i][[1]]$pc.summary$importance[2]*100, digits = 1), "%")) +
      ylab(paste0("PC2 = ", 
                  round(pca_structures[i][[1]]$pc.summary$importance[5]*100, digits = 1), "%")) +
      labs(color = colnames(list_information[i][[1]])[j]) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14))
    plot(p)
  }
  dev.off()
  
  # png individual image for morphotypes category
  png(filename = paste0(resultados[i], names(list_structures)[i], "_PCA_deformation_morphotypes.png"), width = 420*4, height = 290*3, res = 200)
  p <- ggplot(pca_scores, aes(PC1, PC3, color = list_information[i][[1]][,2])) +
    geom_point() +
    stat_ellipse(type = "t") +
    xlab(paste0("PC1 = ", 
                round(pca_structures[i][[1]]$pc.summary$importance[2]*100, digits = 1), "%")) +
    ylab(paste0("PC2 = ", 
                round(pca_structures[i][[1]]$pc.summary$importance[5]*100, digits = 1), "%")) +
    labs(color = colnames(list_information[i][[1]])[2]) +
    theme(#legend.position = "none",
          axis.text=element_text(size=14),
          axis.title=element_text(size=18),
          panel.background = element_rect("white"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "gray90"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "gray90")) 
  
  plot(p)
  dev.off()
  
  # save components
  write.csv(file = paste0(resultados[i], names(list_structures)[i], "_PCAscores.csv"), x = pca_scores)

}





# 6. Deformation grills ####

# separate and analyses morphotypes
list_grills <- list(head = NA, metaesterno = NA, pronoto = NA, protibia = NA) #list to fill

for (i in 1:4) { 
  structures <- list_structures[i][[1]]
  reference <- mshape(structures$coords) #for each structure
  #separate morphotypes
  C_juvencus <- mshape(structures$coords[,,list_information[i][[1]][,2]=="C. juvencus"]) 
  C_septemmaculatus <- mshape(structures$coords[,,list_information[i][[1]][,2]=="C. septemmaculatus"])
  C_subhyalinus <- mshape(structures$coords[,,list_information[i][[1]][,2]=="C. subhyalinus"])
  C_cyanellus <- mshape(structures$coords[,,list_information[i][[1]][,2]=="C. cyanellus"])
  CGm1 <- mshape(structures$coords[,,list_information[i][[1]][,2]=="C. (G.) morfo1"])
  CGm2 <- mshape(structures$coords[,,list_information[i][[1]][,2]=="C. (G.) morfo2"])
  CGm3 <- mshape(structures$coords[,,list_information[i][[1]][,2]=="C. (G.) morfo1"])
  
  list_grills[i][[1]] <- list(Reference = reference, C_juvencus = C_juvencus, 
                         C_septemmaculatus = C_septemmaculatus, C_subhyalinus = C_subhyalinus,
                         C_cyanellus = C_cyanellus, CGm1 = CGm1, CGm2 = CGm2, CGm3 = CGm3)
}

# plot deformations for each structures and morphotypes respect to the reference
for (i in 1:4) { #for each structure
  
  pdf(file = paste0(resultados[i], names(list_structures)[i], "_deformation_grills.pdf"))
  
  for (j in 2:8) { #for each morphotype
    
    TPSgrid(TT = list_grills[i][[1]][1][[1]], YY = list_grills[i][[1]][j][[1]], 
            xbegin=-999, ybegin = -999, xwidth = -999,
            opt=1, ext=0.1, ngrid = 20, cex=1, pch=20, 
            col=1, mag=2.5, axes3=FALSE) #deformation for reference to each morphotype
    title(names(list_grills[i][[1]][j]))
    
  }
  
  dev.off()
  
}

# 7. Allometry ####

#make a list with all information
list_dataFrames_size_shape <- list(head = NA, metaesterno = NA, pronoto = NA, protibia = NA) #list to fill

for (i in 1:4) { #for each strucutre
    
  tmp_size <- list_structures[i][[1]]$Csize #extract the centroide size
  tmp_shape <- list_structures[i][[1]]$coords #extract the shape coordinates
  list_dataFrames_size_shape[i][[1]] <- geomorph.data.frame(size = tmp_size, shape = tmp_shape,
                                                            groups = list_information[i][[1]]$Morfotipo)
}

#regresion between shape and other variables
for (i in 1:4) {
  cab_anova_size_spp <- procD.lm(shape ~ size*groups, data = list_dataFrames_size_shape[i][[1]]) #ANOVA-procrustes size vs shape, groups and shape*groups
  tmp <- summary(cab_anova_size_spp) #the result
  
  write.csv(file = paste0(resultados[i], names(list_structures)[i], "_allometry.csv"), as.data.frame(tmp$table)) #save the statistic values
  #
  pdf(file = paste0(resultados[i], names(list_structures)[i], "_allometry.pdf"))
      
  plot(cab_anova_size_spp, type = "regression", predictor = as.numeric(list_dataFrames_size_shape[i][[1]]$groups), reg.type = "RegScore")
  title("Shape vs Morphotypes")
  
  plot(cab_anova_size_spp, type = "regression", predictor = list_dataFrames_size_shape[i][[1]]$size, reg.type = "RegScore")
  title("Shape vs Size")
  
  plot(cab_anova_size_spp, type = "regression", predictor = as.numeric(list_dataFrames_size_shape[i][[1]]$groups)*list_dataFrames_size_shape[i][[1]]$size, reg.type = "RegScore")
  title("Shape vs Size*Morphotypes")
  
  dev.off()
  #
  
  plotlm <- plot(cab_anova_size_spp, type = "regression", predictor = list_dataFrames_size_shape[i][[1]]$size, reg.type = "RegScore")
  
  DF_ggplot <- data.frame(Reg_shape = plotlm$PredLine,
                          spp = list_dataFrames_size_shape[i][[1]]$groups,
                          CS = list_dataFrames_size_shape[i][[1]]$size)
  
  p <- ggplot(DF_ggplot, aes(x = CS, y = Reg_shape, colour = spp)) + 
    geom_point()
  
  pdf(file = paste0(resultados[i], names(list_structures)[i], "_Reg_shapeVSsize.pdf"))
      plot(p)
  dev.off()
  
}

# 8. Canonical Variable Analyses -CVA- #####

#CVA for test morphotype and locality recognition
for (i in 1:4) {
  
  for (j in c(2,3,5)) {
    
    tmp_cantagallo <- str_which(string = list_information[i][[1]]$Municipio, pattern = "Cantagallo")
    
    #make the CVA, elimnate the group with one entry (individuo 19 from Santo Domingo, Cantagallo)
    cva_structure <- CVA(list_structures[i][[1]]$coords[,,-tmp_cantagallo], list_information[i][[1]][j][[1]][-tmp_cantagallo], 
                      weighting = T, plot = TRUE, 
                      rounds = 0, cv = TRUE, p.adjust.method = "none")
    
    #plot the CVA space and save
    p <- ggplot(as.data.frame(cva_structure$CVscores), 
                aes(x = cva_structure$CVscores[,1], 
                    y = cva_structure$CVscores[,2],
                    color = as.factor(row.names(cva_structure$CVscores)))) + 
      geom_point() +
      stat_ellipse(type = "t") +
      xlab(paste0("CV1 = ", 
                  round(cva_structure$Var[1,2], digits = 1), "%")) +
      ylab(paste0("CV2 = ", 
                  round(cva_structure$Var[2,2], digits = 1), "%")) +
      labs(color = names(list_information[i][[1]][j])) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14))
    
    png(filename = paste0(resultados[i], names(list_structures)[i], "_CVA_", names(list_information[i][[1]][j]), ".png"), width = 420*4, height = 290*3, res = 200)
      plot(p)
    dev.off()
    
    #save the table of frecuencies hits
    tmp_cva <- print(cva_structure)
    
    write.csv(file = paste0(resultados[i], names(list_structures)[i], "_CVA_", names(list_information[i][[1]][j]), ".csv"), x = tmp_cva$table )
  } 
}

# 9. Cluster UPGMA ####
# colLab function of Joris Meys (http://stackoverflow.com/) modified by Ambrosio Torres to coloring edgePar and labels of a dendrogram ####

colLab <- function(n) {
  if(is.leaf(n)) {
    a <- attributes(n)
    # clusMember - a vector designating leaf grouping
    # labelColors - a vector of colors for the above grouping
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "edgePar") <- list(col = labCol)
    attr(n, "label") <- NULL ;
  }
  n
}

# 9.1 by morphotypes ####
for (i in 1:4){
  datavalores <- pca_structures[i][[1]]$pc.scores
  row.names(datavalores) <- list_information[i][[1]]$ID
  
  distvalores <- dist(datavalores) 
  clust_valores <- hclust(distvalores, method = "average") 
  labelColors <- c("#F8766D", "#BB9D00", "#00B81F", "#00C0B8", "#00A5FF", "#E76BF3", "#FF6C90")
  
  clusMember <- rep(NA,length(rownames(datavalores)))
  
  for (j in 1:length(list_information[i][[1]]$Morfotipo)) {
    
    if (list_information[i][[1]]$Morfotipo[j] == "C. cyanellus") {
      clusMember[j] <- 1
    }
    if (list_information[i][[1]]$Morfotipo[j] == "C. juvencus") {
      clusMember[j] <- 2
    }
    if (list_information[i][[1]]$Morfotipo[j] == "C. septemmaculatus") {
      clusMember[j] <- 3
    }
    if (list_information[i][[1]]$Morfotipo[j] == "C. subhyalinus") {
      clusMember[j] <- 4
    }
    if (list_information[i][[1]]$Morfotipo[j] == "C. (G.) morfo1") {
      clusMember[j] <- 5
    }
    if (list_information[i][[1]]$Morfotipo[j] == "C. (G.) morfo2") {
      clusMember[j] <-6
    }
    if (list_information[i][[1]]$Morfotipo[j] == "C. (G.) morfo3") {
      clusMember[j] <- 7
    }
    
  }
  
  names(clusMember) <- rownames(datavalores)
  
  clusDendro <- as.dendrogram(as.hclust(clust_valores))  
  
  clusDendro <- dendrapply(clusDendro, colLab)
  
  #make graph and save it
  
  png(filename = paste0(resultados[i], names(list_structures)[i], "_Dendrogram_morphotypes.png"), width = 420*4, height = 290*3, res = 200)
    par(mar=c(0.5, 4.5, 1, 1))
    par(oma= c(0,0,0,0))
    par(lwd=2)
    plot(clusDendro,horiz=F,axes=T, ylab= "Euclidean distance", cex.axis=1.3,cex.lab=1.7)
    par(lwd=1)
    legend("topright", pch= 21, pt.bg=labelColors, 
    legend=expression(italic("C. cyanellus"), italic("C. juvencus"),italic("C. septemmaculatus"), italic("C. subhyalinus"), "C. (G.) morfo1", "C. (G.) morfo2", "C. (G.) morfo2"), pt.cex = .9, cex=.6)
  dev.off()
}

# 9.2 by municipality ####
for (i in 1:4){
  datavalores <- pca_structures[i][[1]]$pc.scores
  row.names(datavalores) <- list_information[i][[1]]$ID
  
  distvalores <- dist(datavalores) 
  clust_valores <- hclust(distvalores, method = "average") 
  labelColors <- sequential_hcl(4, palette = "Terrain")
  
  clusMember <- rep(1,length(datavalores))
  
  for (j in 1:length(list_information[i][[1]]$Municipio)) {
    
    if (list_information[i][[1]]$Municipio[j] == "Cantagallo") {
      clusMember[j] <- 1
    }
    if (list_information[i][[1]]$Municipio[j] == "Carmen del Chucurí") {
      clusMember[j] <- 2
    }
    if (list_information[i][[1]]$Municipio[j] == "Cimitarra") {
      clusMember[j] <- 3
    }
    if (list_information[i][[1]]$Municipio[j] == "Santa Bárbara") {
      clusMember[j] <- 4
    }
  }
  
  names(clusMember) <- rownames(datavalores)
  
  clusDendro <- as.dendrogram(as.hclust(clust_valores))  
  
  clusDendro <- dendrapply(clusDendro, colLab)
  
  #make graph and save it
  
  png(filename = paste0(resultados[i], names(list_structures)[i], "_Dendrogram_munic.png"), width = 420*4, height = 290*3, res = 200)
  par(mar=c(0.5, 4.5, 1, 1))
  par(oma= c(0,0,0,0))
  par(lwd=2)
  plot(clusDendro,horiz=F,axes=T, ylab= "Euclidean distance", cex.axis=1.3,cex.lab=1.7)
  par(lwd=1)
  legend("topright", pch= 21, pt.bg=labelColors, 
         legend=expression("Cantagallo", "Carmen de Chucurí", "Cimitarra", "Santa Bárbara"), pt.cex = .9, cex=.6)
  dev.off()
}

# 9.3 by localities ####
for (i in 1:4){
  datavalores <- pca_structures[i][[1]]$pc.scores
  row.names(datavalores) <- list_information[i][[1]]$ID
  
  distvalores <- dist(datavalores) 
  clust_valores <- hclust(distvalores, method = "average") 
  labelColors <- sequential_hcl(7, palette = "Plasma")
  
  clusMember <- rep(1,length(datavalores))
  
  for (j in 1:length(list_information[i][[1]]$Localidad)) {
    
    if (list_information[i][[1]]$Localidad[j] == "Bocas del Carare") {
      clusMember[j] <- 1
    }
    if (list_information[i][[1]]$Localidad[j] == "Santo Domingo") {
      clusMember[j] <- 2
    }
    if (list_information[i][[1]]$Localidad[j] == "El Tigre") {
      clusMember[j] <- 3
    }
    if (list_information[i][[1]]$Localidad[j] == "La Belleza") {
      clusMember[j] <- 4
    }
    if (list_information[i][[1]]$Localidad[j] == "La Dorada") {
      clusMember[j] <- 5
    }
    if (list_information[i][[1]]$Localidad[j] == "Salina") {
      clusMember[j] <- 6
    }
    if (list_information[i][[1]]$Localidad[j] == "Vereda Salinas") {
      clusMember[j] <- 7
    }
  }
  
  names(clusMember) <- rownames(datavalores)
  
  clusDendro <- as.dendrogram(as.hclust(clust_valores))  
  
  clusDendro <- dendrapply(clusDendro, colLab)
  
  #make graph and save it
  
  png(filename = paste0(resultados[i], names(list_structures)[i], "_Dendrogram_local.png"), width = 420*4, height = 290*3, res = 200)
  par(mar=c(0.5, 4.5, 1, 1))
  par(oma= c(0,0,0,0))
  par(lwd=2)
  plot(clusDendro,horiz=F,axes=T, ylab= "Euclidean distance", cex.axis=1.3,cex.lab=1.7)
  par(lwd=1)
  legend("topright", pch= 21, pt.bg=labelColors, 
         legend=expression("Bocas del Carare", "Santo Domingo", "El Tigre", "La Belleza", "La Dorada", "La Bodega", "La Salina"), pt.cex = .9, cex=.6)
  dev.off()
}

# 10. K-means analyses ####
############# optimal cluster by Gap statistic
for (i in 1:4) {
  tmp_data <- two.d.array(list_structures[i][[1]]$coords)
  png(filename = paste0(resultados[i], names(list_structures)[i], "_GAP_kmeans.png"), width = 420*4, height = 290*3, res = 200)
  plot(fviz_nbclust(tmp_data, kmeans, method = "gap_stat", nboot = 500))
  dev.off()
}
# head has 5 k-groups 
# Metaesternum has 3 k-groups
# Pronotum has 5 k-groups
# Protibia has 7 k-groups

############# K-groups plots
kgroups <- c(head = 5, metaesternun = 3, pronoto = 5, protibia = 7) #groups recognition before

for (i in 1:4) {
  
  tmp_data <- two.d.array(list_structures[i][[1]]$coords)
  km.res <- kmeans(tmp_data, centers = kgroups[i], nstart = 25, iter.max = 500) #calculate the data.frame
  
  TMP_ROW <- nrow(list_information[i][[1]])
  morfos <- vector()
  for (j in 0:TMP_ROW) {
    morfos[j] <- paste0(str_sub(list_information[i][[1]]$Morfotipo[j], 1, 7), j)
  }
  rownames(tmp_data) <- morfos
  
  #### in PCA space 
  p <- fviz_cluster(km.res, data = tmp_data,
               ellipse.type = "norm",
               palette = "jco",
               ggtheme = theme_minimal(),
               ellipse = TRUE, shape = NULL, 
               ellipse.level = 0.95, ellipse.alpha = 0.2)
  
  png(filename = paste0(resultados[i], names(list_structures)[i], "_PCA_Kgroups.png"), width = 420*4, height = 290*3, res = 200)
    plot(p)
  dev.off()
  
  #### in cluster (UPGMA) 
  res.hc <- tmp_data %>%
    scale() %>%                    # Scale the data
    dist(method = "euclidean") %>% # Compute dissimilarity matrix
    hclust(method = "average")     # Compute hierachical clustering UPGMA
  
  p <- fviz_dend(res.hc, k = kgroups[i], # Cut in four groups
            cex = 0.3, # label size
            color_labels_by_k = TRUE, # color labels by groups
            rect = F )# Add rectangle around groups
  
  png(filename = paste0(resultados[i], names(list_structures)[i], "_cluster_Kgroups.png"), width = 420*4, height = 290*3, res = 200)
    plot(p)
  dev.off()
}
