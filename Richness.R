## Andres Ordoñez & Daniel Rafael Miranda-Esquivel (Laboratorio de Sistemática y Biogeografía, Universidad Industrial de Santander, Bucaramanga, Colombia)
## Bachelor Thesis: Composition and morphometric analysis of Canthon genus in Santander and south of Bolívar localities 
## Part: Diversity and composition
## R version 3.6.3
## May 2020

# 1. Load packeges ####

library(vegan)
library(philentropy)
library(ggplot2)
library(iNEXT)


# 2. Create directories ######
dir.create("./Diversity/Richness", showWarnings = FALSE) #create a directory
dir.create("./Diversity/Abundance", showWarnings = FALSE) #create a directory

# 3. read Data, just need one information table of any structure
info <- read.csv("Head_results/info_head.csv")

# 4. calculate and plot richness and abundance for localities, municipality, 500mts, 1000mts and lands type
for (i in c(3, 5, 7:9)) {

  ######## abundance
  abundance <- ftable(info$Morfotipo ~ info[,i])
  abundance <- as.table(abundance)
  ######## bar plot
  for_gg <- data.frame(abundance = apply(abundance, 1, sum),
                       localidades = row.names(abundance))
  
  p <- ggplot(data = for_gg, aes(x = localidades, y = abundance)) +
    geom_col(fill = "#56B4E9", position = "stack") +
    xlab("Localidades") +
    ylab("Riqueza") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          panel.background = element_rect("white"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "gray90"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "gray90"))
  png(paste0("./Diversity/Abundance/abundance_", names(info)[i], ".png"), width = 420*3, height = 290*3, res = 200)
  plot(p)
  dev.off()
  
  ######## richness
  richness <- abundance
  richness[richness > 1] <- 1
  
  ######## bar plot
  for_gg <- data.frame(riqueza = apply(richness, 1, sum),
                       localidades = row.names(richness))
  
  p <- ggplot(data = for_gg, aes(x = localidades, y = riqueza)) +
    geom_col(fill = "#56B4E9", position = "stack") +
    xlab("Localidades") +
    ylab("Riqueza") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          panel.background = element_rect("white"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "gray90"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "gray90"))
  png(paste0("./Diversity/Richness/rich_", names(info)[i], ".png"), width = 420*3, height = 290*3, res = 200)
  plot(p)
  dev.off()
  
  ######## Jaccard and Sorensen distances
  jaccard <- vegdist(richness, method = "jaccard")
  sorensen <- vegdist(richness, method = "bray")
  jac_clust <- hclust(jaccard)
  png(paste0("./Diversity/Richness/", "Jaccard_", names(info)[i], ".png"), width = 420*3, height = 290*3, res = 200)
  plot(jac_clust)
  dev.off()
}


# Extrapolation for municipality #####
abundancia <- ftable(info$Morfotipo ~ info$Altura)
abundancia <- as.table(abundancia)
riqueza <- abundancia

riqueza[riqueza > 1] <- 1

row.names(riqueza)

Cantagallo_r <- rbind(riqueza["105",])
Cimitarra_r <-  rbind(riqueza[c(1,2,4:16),])
Carmen_r <- rbind(riqueza[c(17:37),])
StaBarbara_r <- rbind(riqueza[c(38:47),])

row.names(Cantagallo_r) <- NULL
row.names(Cimitarra_r) <- NULL
row.names(Carmen_r) <- NULL
row.names(StaBarbara_r) <- NULL

# volverlas vectores
Cimitarra_v <- apply(Cimitarra_r, 1, sum)
Cantagallo_v <- apply(Cantagallo_r, 1, sum)
Carmen_v <- apply(Carmen_r, 1, sum)
StaBarbara_v <- apply(StaBarbara_r, 1, sum)

# rellenar hasta completar el num de trampas usadas en campo, tabla en el informe
Cimitarra_c <- c(63, Cimitarra_v)
Cantagallo_c <- c(44, Cantagallo_v)
Carmen_c <- c(44, Carmen_v)
StaBarbara_c <- c(20, StaBarbara_v)

#hacer la lista que los reuna a todos
list.riqueza <- list(Cimitarra = Cimitarra_c,
                     #Cantagallo = Cantagallo_c,
                     Carmen = Carmen_c,
                     StaBarbara = StaBarbara_c)
# iNEXT
t <- seq(1, 100, by=5)
out.inc <- iNEXT(list.riqueza, q=0, datatype="incidence_freq", size = t)

png("./Diversity/rarefaccion.png", width = 420*3, height = 290*3, res = 200)
ggiNEXT(out.inc, type=1) + 
  ylab("Diversidad de especies") +
  xlab("Unidades de muestreo") +
  theme(#legend.position="none", 
        axis.text=element_text(size=14),
        axis.title=element_text(size=18),
        panel.background = element_rect("white"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                        colour = "gray90"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "gray90")) +
  scale_shape_manual(values = c(20,20,20))
dev.off()

