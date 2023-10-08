install.packages("tkWidgets")
library(affy)
library(affyPLM)
library(limma)

#Setear el directorio donde trabjaremos los datos
setwd("C:\\Users\\Paulw\\Desktop\\TODO\\rstudioarchivos\\Cancer_Data")

# Funcion que lee todos los archivos .CEL en el directorio seleccionado anteriormente
#celFiles <- list.celfiles()
# Funcion para leer los archivos .CEL que guardamos anteriormente
affyRaw <- read.celfiles(list.celfiles())

microarray.raw.data <- ReadAffy(affyRaw)

boxplot(microarray.raw.data,col= rainbow(8))

eset <- rma(microarray.raw.data)

boxplot(eset,col= rainbow(8))

affyRaw

image(microarray.raw.data[,1],col=rainbow(100))

microarray.processed.data <- rma(microarray.raw.data)

boxplot(microarray.processed.data,col=rainbow(8),las = 2, ylab = "expresion")


expression.level <- exprs(eset)
head(expression.level)
dim(expression.level)
sampleID <- c("PRBT 7411","PRBT 7412","PRBT 7415","PRBT 7421","PRBT 7422","PRBT 8159",
              "PRBT 8170","PRBT 8189","PRBT 8191","PRBT 8199","PRBT 8200","PRBT 8209",
              "PRBT 8210","PRBT 8413","PRBT 11519","PRBT 11524","PRBT 11526","PRBT 11529","PRBT 11530",
              "PRBT 11532","PRBT 11571","PRBT 11582","PRBT 11585", "PRBT final")

colnames(expression.level) <- sampleID

colnames(expression.level)
head(expression.level)


treat1 <- (expression.level[,1] + 
             expression.level[,2])/2
treat2 <- (expression.level[,3] + 
             expression.level[,4])/2
treat3 <- (expression.level[,5] + 
             expression.level[,6])/2
treat4 <- (expression.level[,7] + 
             expression.level[,8])/2
treat5 <- (expression.level[,9] + 
             expression.level[,10])/2
treat6 <- (expression.level[,11] + 
             expression.level[,12])/2
treat7 <- (expression.level[,13] + 
             expression.level[,14])/2
treat8 <- (expression.level[,15] + 
             expression.level[,16])/2
treat9 <- (expression.level[,17] + 
             expression.level[,18])/2
treat10 <- (expression.level[,19] + 
             expression.level[,20])/2
treat11 <- (expression.level[,21] + 
             expression.level[,22])/2
treat12 <- (expression.level[,23] + 
             expression.level[,24])/2

mean.treats <- rowMeans(expression.level)
conditions.id <- c("treat1","treat2","treat3","treat4","treat5","treat6","treat7","treat8","treat9","treat10","treat11","treat12")

mean.expression <- matrix(mean.treats, ncol= 25)

mean.expression <- matrix(c(treat1,treat2,treat3,treat4,treat5,treat6,treat7,treat8
                            ,treat9,treat10,treat11,treat12),ncol = 12)

matrix.mean.expression <- matrix(mean.expression,ncol = 24)


rownames(mean.expression) <- names(treat1)
colnames(mean.expression) <- conditions.id
colnames(matrix.mean.expression) <- sampleID

head(mean.expression)
head(matrix.mean.expression)

plot(treat1,treat2,xlab = "Treat 1", ylab = "Treat 2", pch= 19,cex= 0.5)



experimental.design <- model.matrix(~ -1+factor(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)))

colnames(experimental.design) <- sampleID


linear.fit <- lmFit(expression.level,experimental.design)
head(linear.fit)

contrast.matrix <- makeContrasts(treat1-treat2,
                                 treat3-treat4,
                                 treat5-treat6,
                                 treat7-treat8,
                                 treat9-treat10,
                                 treat11-treat12,
                                 levels = conditions.id)

dim(contrast.matrix)
dim(linear.fit)

contrast.linear.fit <- contrasts.fit(linear.fit,contrast.matrix)
contrast.result <- eBayes(contrast.linear.fit)

