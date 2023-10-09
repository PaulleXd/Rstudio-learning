#Haciendo un analisis de expresion de microarreglos usando un heatmap
#cargar librerias

BiocManager::install("affy")
BiocManager::install("oligo")
BiocManager::install("Biobase")
BiocManager::install("GEOquery")
BiocManager::install("arrayQualityMetrics")
BiocManager::install("a")
BiocManager::install("simpleaffy")
BiocManager::install("limma")
install.packages("splitstackshape")
install.packages("caret")
install.packages("pheatmap")

help("repositories",package = "BiocManager")

library("pheatmap")
require("ggplot2")
require("colorspace")
library("grid")
library(limma)

library("GEOquery")
library(oligo)
library(Biobase)
library(affy)
library("splitstackshape")
library("tidyr")
library("dplyr")
library("arrayQualityMetrics")
library(affyPLM)
library(affy)

library(DESeq2)




#Setear el directorio donde trabjaremos los datos
setwd("C:\\Users\\Paulw\\Desktop\\TODO\\rstudioarchivos\\Cancer_Data")

# Funcion que lee todos los archivos .CEL en el directorio seleccionado anteriormente
#celFiles <- list.celfiles()
# Funcion para leer los archivos .CEL que guardamos anteriormente
affyRaw <- read.celfiles(list.celfiles())

#normalizacion de los datos cargados en la variable affyRaw (caracterizacion de los datos, normalizacion y calculacion de la expresion)
eset <- oligo::rma(affyRaw)


#guardar los datos en un archivo
write.exprs(eset, file = "data1.txt")

#guardar los datos de el archivo GSE19697
mydata <- cSplit(data.frame(read.delim("GSE19697_family.soft",check.names = FALSE)),"Gene.Symbol","//")

# transformar los datos a un df
# df_mydata <- data.frame(mydata)

# separar los "gene symbol" con //
a <- cSplit(df_mydata,"Gene.Symbol","//")

#guardarlo en un archivo
write.table(a, file = "importData1.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# guardar los dataframes exportados
data_a <- read.delim("importData1.txt",check.names = FALSE)
data_b <- read.delim("data1.txt",check.names = FALSE)

# combinar los datos
data_combined <- left_join(data_a, data_b, by = "ID" )

# guardar la data en un csv
write.csv(data_combined, "combinedData.csv")


boxplot(eset,col=rainbow(8),las=2,ylab="datos")

hist(eset)

expression.level <- exprs(eset)

head(expression.level)
dim(expression.level)



sampleID <- c("PRBT 7411","PRBT 7412","PRBT 7415","PRBT 7421","PRBT 7422","PRBT 8159",
              "PRBT 8170","PRBT 8189","PRBT 8191","PRBT 8199","PRBT 8200","PRBT 8209",
              "PRBT 8210","PRBT 8413","PRBT 11519","PRBT 11524","PRBT 11526","PRBT 11529","PRBT 11530",
              "PRBT 11532","PRBT 11571","PRBT 11582","PRBT 11585", "PRBT final")

colnames(expression.level)


colnames(expression.level) <- sampleID

colnames(expression.level)

mean.treatment1 <- (expression.level[,1]+expression.level[,2])/2
mean.treatment2 <- (expression.level[,3]+expression.level[,4])/2
mean.treatment3 <- (expression.level[,5]+expression.level[,6])/2
mean.treatment4 <- (expression.level[,7]+expression.level[,8])/2
mean.treatment5 <- (expression.level[,9]+expression.level[,10])/2
mean.treatment6 <- (expression.level[,11]+expression.level[,12])/2
mean.treatment7 <- (expression.level[,13]+expression.level[,14])/2
mean.treatment8 <- (expression.level[,15]+expression.level[,16])/2
mean.treatment9 <- (expression.level[,17]+expression.level[,18])/2
mean.treatment10 <- (expression.level[,19]+expression.level[,20])/2
mean.treatment11 <- (expression.level[,21]+expression.level[,22])/2
mean.treatment12 <- (expression.level[,23]+expression.level[,24])/2


mean.treatment <- (expression.level[,1] + expression.level[,2] + expression.level[,3] + expression.level[,4] + expression.level[,5] + 
                     expression.level[,6] + expression.level[,7] + expression.level[,8] + expression.level[,9] + expression.level[,10] + 
                     expression.level[,11] + expression.level[,12] + expression.level[,13] + expression.level[,14] + expression.level[,15] + 
                     expression.level[,16] + expression.level[,17] + expression.level[,18] + expression.level[,19] + expression.level[,20] + 
                     expression.level[,21] + expression.level[,22] + expression.level[,23] + expression.level[,24])/24




head(mean.treatment)






plot(PRBT_1_2,PRBT_3_4,xlab = "PRBT 1", ylab = "PRBT 2",pch=19,cex=0.5)


# SELECCION DE GENES DIFERENCIALMENTE EXPRESADO
# Seleccion basada en el Factor de proporcionalidad donde se calcula el nivel de expresion medio

 

experimental.desing <- model.matrix(~ -1+factor(c(1:24)))

colnames(experimental.desing) <- sampleID

colnames(experimental.desing)


linear.fit <- lmFit(eset,experimental.desing)

constrat.matrix <- makeContrasts(mean.treatment1-mean.treatment2,
                                 mean.treatment3-mean.treatment4,
                                 mean.treatment5-mean.treatment6,
                                 mean.treatment7-mean.treatment8,
                                 mean.treatment9-mean.treatment10,
                                 mean.treatment11-mean.treatment12,
                                 mean.treatment13-mean.treatment14,
                                 mean.treatment15-mean.treatment16,
                                 mean.treatment17-mean.treatment18,
                                 mean.treatment19-mean.treatment20,
                                 mean.treatment21-mean.treatment22,
                                 mean.treatment23-mean.treatment24,
                                 levels = c("Treatment1","Treatment2","Treatment3","Treatment4","Treatment5","Treatment6",
                                            "Treatment7","Treatment8","Treatment9","Treatment10","Treatment11","Treatment12"))

constrat.matrix <- makeContrasts(mean.treatment,levels = "Tratamiento")





