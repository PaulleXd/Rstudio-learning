# Establecer un directorio de trabajo, depende de la computadora de cada persona y donde quiera establecer la ruta para trabjar
# sus archivos
setwd("C:/Users/Paulw/Desktop/Arbol filogenetico/cordyceps data")
# instalar los paquetes necesarios
install.packages("remotes")
remotes::install_github("GuangchuangYu/treeio")
install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("DECIPHER")
install.packages("ape")

library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)

# Establecer el archivo a partir del cual tomaremos las secuencias
# usar la funcion readRNAStringSet() o readAAStringSet() en caso de ser necesario 
seq_name <- file.choose()

seqs <- readDNAStringSet(seq_name, format = "fasta")

# visualizar las secuencias
seqs

# Orienta las secuencias en una misma direccion
seqs <- OrientNucleotides(seqs)

# realizar un alineamiento de las secuencias
aligned <- AlignSeqs(seqs)

# Vizualisar las secuencias en el explorador de internet
BrowseSeqs(aligned,colorPatterns = TRUE, highlight=0)

# guardar las secuencias alineadas en un archivo .fasta
writeXStringSet(aligned,
                file="Cordyceps.fasta")

# a partir de aqui usaremos el archivo guardado en formato .fasta para trabajar
dna <- read.alignment("Cordyceps.fasta", format = "fasta")

# Matriz de distanciamiento, 
D <- dist.alignment(dna, matrix = "similarity")


temp <- as.data.frame(as.matrix(D))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5)+

# we can start to see a pattern because the data is ordered by year, 
# but we can't really make any conclusions yet
print(D)
class(D)
  
tre <- nj(D)
class(tre) #all trees created using {ape} package will be of class phylo

tre <- ladderize(tre)

# ~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~Base R plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot(tre, cex = 0.6)
title("Similitud en Cordyceps de region (ITS)")


# or 
h_cluster <- hclust(D, method = "average", members = NULL) # method = average is used for UPGMA, members can be equal to NULL or a vector with a length of size D
plot(h_cluster, cex = 0.6)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tree Plotting in ggtree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# you can fan it out 
ggtree(tre, yscale = "NA")+
  geom_tiplab(hjust = -0.3, size=4, align = TRUE)+
  xlim(0,0.5) 

# or whatever this thing does???
ggtree(tre,layout = "daylight")+
  geom_tiplab(hjust = -0.3, size=4, align = TRUE)+
  xlim(0,0.5) 

# plot a basic tree
ggtree(tre) + 
  geom_tiplab(hjust = -0.3, size=4, align = TRUE)+
  xlim(0,0.5)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~ Customize your trees ~~~~~~~~~~~~~~~~~~~~~~~~
# plot using ggtree and highlight clusters
# change the node values for your own data
ggtree(tre) + 
  geom_tiplab(hjust = -0.3, size=4, align = TRUE) + 
  geom_hilight(node=19, fill="purple", alpha = 0.2) + 
  geom_hilight(node=17, fill="dark green", alpha = 0.2) +
  geom_hilight(node=20, fill="gold", alpha = 0.2) +
  xlim(0,0.5) 

# highlight clusters and add a vertical line to group clusters
# change the node values for your own data
ggtree(tre) + 
  geom_tiplab(hjust = -0.3, size=4, align = TRUE) + 
  geom_hilight(node=19, fill="purple", alpha = 0.2) + 
  geom_hilight(node=17, fill="dark green", alpha = 0.2) +
  geom_hilight(node=20, fill="gold", alpha = 0.2) +
  geom_cladelabel(node=19, label=" Cluster 1", 
                  color="purple", offset=.1, barsize = 2,
                  fontsize = 5, align=TRUE, alpha = 0.5) + 
  geom_cladelabel(node=17, label=" Cluster 2", 
                  color="dark green", offset=.1, barsize = 2,
                  fontsize = 5, align=TRUE, alpha = 0.5) + 
  geom_cladelabel(node=20, label=" Cluster 3", 
                  color="gold", offset=.1, barsize = 2,
                  fontsize = 5, align=TRUE, alpha = 0.5) + 
  xlim(0,0.5) 



# ~~~~~~~~~~~~~~~~~~~~~~~~~~ Plot the allignment with the tree ~~~~~~~~~~~~~~~~~

# lets plot the alignment with the tree, to do this we first have to
# match the names to the tip labels
# set our tree into a new name
tre.new <- tre
# change tip labels to full alignment names
tre.new$tip.label <- aligned@ranges@NAMES

# plot the alignment 
  msaplot(p=ggtree(tre.new), fasta="Cordyceps.fasta", 
        window=c(150, 175))+
  scale_fill_viridis_d(alpha = 0.8)

