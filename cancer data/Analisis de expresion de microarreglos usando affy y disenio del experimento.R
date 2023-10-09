library(affy)
library(limma)
library(affy)



setwd("C:\\Users\\Paulw\\Desktop\\TODO\\rstudioarchivos\\Cancer_Data")

cel_files <- list.celfiles()

eset <- ReadAffy(filenames = cel_files)

#preprocesamiento de datos, que incluye background correction, normalización, y resumen de las sondas
eset <- rma(eset)

#Diseño experimental
tratamientos_originales <- c("tratamiento1", "tratamiento2","tratamiento3","tratamiento4","tratamiento5","tratamiento6","tratamiento7",
                  "tratamiento8","tratamiento9","tratamiento10","tratamiento11","tratamiento12","tratamiento13","tratamiento14",
                  "tratamiento15","tratamiento16","tratamiento17","tratamiento18","tratamiento19","tratamiento20","tratamiento21",
                  "tratamiento22","tratamiento23","tratamiento24")
tratamientos <- make.names(tratamientos_originales)

design <- model.matrix(~0 + factor(tratamientos))

colnames(design) <- tratamientos_originales

class(design)



#Ajuste del modelo con los datos
fit <- lmFit(eset, design)

# Contraste de interés
contraste <- makeContrasts(tratamiento1 - tratamiento2, levels = design)

# Ajuste del modelo con el contraste
fit2 <- contrasts.fit(fit, contrasts = contraste)

fit2 <- eBayes(fit2)

ma_plot <- plotMA(fit2)
write.csv(fit, file = "fit")

ma_plot <- limma::plotMA(fit2, main = "Grafico")
pdf("ma_plot.pdf")
print(ma_plot)      
dev.off()           
