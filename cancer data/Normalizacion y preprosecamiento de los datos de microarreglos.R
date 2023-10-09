setwd("C:\\Users\\Paulw\\Desktop\\TODO\\rstudioarchivos\\Cancer_Data")

install.packages("matrixStats")
install.packages("genefilter")
install.packages("car")

library(tidyverse)
library(oligo)
library(matrixStats)
library(genefilter)
library(oligoClasses)
library(magrittr)
library("car")

celFiles <- list.celfiles()

rawData <- read.celfiles(celFiles)

class(rawData)

print(rawData)

image(rawData[,24],col=rainbow(100))

hist(rawData, col=rainbow(30), main="microarreglos")

boxplot(rawData, col=rainbow(30), las=2, main="miroarreglos")

normalData <- rma(rawData)

class(normalData)

image(normalData[,1],col=rainbow(100))

hist(normalData, col=rainbow(30), main="microarreglos")

boxplot(normalData, col=rainbow(30), las=2, main="miroarreglos")


for (i in 1:24) {
  nombre_salida <- paste(celFiles[i], ".png", sep="")
  png(nombre_salida)
  image(rawData[,i], col=rainbow(100))
  dev.off()
  
  cat("Imagen", i, "procesada y guardada como", nombre_salida, "\n")
}


pmSeq <- pmSequence(rawData)

pmSeq[1:24]

counts <- Biostrings::alphabetFrequency(pmSeq, baseOnly=TRUE)

GCcontent <- ordered(counts[, "G"]+counts[, "C"])


fit1 <- fitProbeLevelModel(rawData)

coef(fit1)[1:4, 1:2]
se(fit1)[1:4, 1:2]

RLE(fit1, type='stats')[, 1:2]

RLE(fit1, type='values')[1:4, 1:2]

NUSE(fit1, type='stats')[, 1:2]

NUSE(fit1, type='values')[1:4, 1:2]

nuseData <- NUSE(fit1, type='values')
class(nuseData)
class(rawData)
for (i in 1:24) {
  nombre_salida <- paste(celFiles[i], ".png", sep="")
  png(nombre_salida)
  image(as.matrix(nuseData[,]), col=rainbow(100))
  dev.off()
  
  cat("Imagen", i, "procesada y guardada como", nombre_salida, "\n")
}


class(normalData)


dim(normalData)

rownames(normalData)

summary(normalData)

colnames(normalData)

names(normalData)

dfData <- as.data.frame(normalData)

image(as.matrix(dfData[1,]),col=rainbow(100))

boxplot(as.matrix(dfData), col=rainbow(30), las=2, main="miroarreglos")

image(as.matrix(dfData[1,]),col=rainbow(100))


rawData1 <- rawData[,1]

image(rawData1,col=rainbow(100))

plot(rawData1[1,], y = NULL)
class(rawData1)
exprs(rawData)[1,]

exprs(normalData)[1,]

rawdfData1 <- as.matrix(exprs(normalData)[,1:24])

rawdfData <-  as.matrix(exprs(normalData))

plot(rawdfData1[,1], rawdfData1[,2],
     main = "grafico de distribucion",
     xlab = "GSM491610",
     ylab = "GSM491611",
     pch = 19, 
     frame = FALSE,
     col = rainbow(100)
)

value <- rawdfData1[,1]

ggplot(as.data.frame(rawdfData1), aes(x = rawdfData1[,1], y = rawdfData1[,2], color = value)) +
  geom_point() +
  scale_color_gradient(low = "green", high = "red") +
  labs(title = "Grafico de distribucion",
       x = "GSM491610",
       y = "GSM491611")


ggplot(as.data.frame(rawdfData), aes(x = rawdfData1[,1], y = rawdfData1[,2], color = value)) +
  geom_point() +
  scale_color_gradient(low = "green", high = "red") +
  labs(title = "Grafico de distribucion",
       x = "GSM491610",
       y = "GSM491611")
