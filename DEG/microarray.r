library(affy)
#library(simpleaffy)
library(affyPLM)
library(limma)
library(annaffy)
library(gplots)
library(arrayQualityMetrics)
library(ath1121501.db)



rawData <- ReadAffy(verbose=TRUE)

image(rawData[,1],col=rainbow(100))

arrayQualityMetrics(expressionset = rawData, outdir = "Quality_Metrics", force = TRUE, do.logtransform =  TRUE)

boxplot(rawData,col=rainbow(30),las=2,ylab="Fluorescencia")

hist(rawData,col=rainbow(30))

eset <- rma(rawData)

boxplot(eset,col=rainbow(30),las=2,ylab="Fluorescencia")

hist(eset,col=rainbow(30))


expression.level <- exprs(eset)
head(expression.level)
dim(expression.level)

samplesNames <- c("WT_with_Fe_1","WT_with_Fe_2","WT_no_Fe_1","WT_no_Fe_2",
                  "pye_with_Fe_1","pye_with_Fe_2","pye_no_Fe_1","pye_no_Fe_2")
colnames(expression.level) <- samplesNames
head(expression.level)

wt.with.fe <- (expression.level[,"WT_with_Fe_1"] + 
                 expression.level[,"WT_with_Fe_2"])/2

wt.no.fe <- (expression.level[,"WT_no_Fe_1"] + 
               expression.level[,"WT_no_Fe_2"])/2

pye.with.fe <- (expression.level[,"pye_with_Fe_1"] + 
                  expression.level[,"pye_with_Fe_2"])/2

pye.no.fe <- (expression.level[,"pye_no_Fe_1"] + 
                expression.level[,"pye_no_Fe_2"])/2

mean.expression <- matrix(c(wt.with.fe,wt.no.fe,
                            pye.with.fe,pye.no.fe),ncol=4)

conditions <- c("WT_with_Fe","WT_no_Fe","pye_with_Fe","pye_no_Fe")
rownames(mean.expression) <- names(wt.with.fe)
colnames(mean.expression) <- conditions
head(mean.expression)

plot(wt.with.fe,wt.no.fe,xlab="WT with Fe",ylab="WT no Fe",pch=19,cex=0.5)
plot(pye.with.fe,pye.no.fe,xlab="pye with Fe",ylab="pye no Fe",pch=19,cex=0.5)
plot(wt.with.fe,pye.with.fe,xlab="WT with Fe",ylab="pye with Fe",pch=19,cex=0.5)
plot(wt.no.fe,pye.no.fe,xlab="WT no Fe",ylab="pye no Fe",pch=19,cex=0.5)


eDesign <- model.matrix(~ -1+factor(c(1,1,2,2,3,3,4,4)))
colnames(eDesign) <- c("WT_with_Fe","WT_no_Fe",
                       "pye_with_Fe","pye_no_Fe")

linear.fit <- lmFit(eset,eDesign)

contrast.matrix <- makeContrasts(WT_no_Fe-WT_with_Fe,
                                 pye_no_Fe-pye_with_Fe,
                                 pye_with_Fe-WT_with_Fe,
                                 pye_no_Fe-WT_no_Fe,
                                 levels=c("WT_with_Fe","WT_no_Fe",
                                          "pye_with_Fe","pye_no_Fe"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)


WT.with.no.Fe <- topTable(contrast.results, number=22810,coef=1,sort.by="logFC")
head(WT.with.no.Fe)

pye.with.no.Fe <- topTable(contrast.results, number=22810,coef=2,sort.by="logFC")
head(pye.with.no.Fe)

fold.change.WT.with.no.Fe <- WT.with.no.Fe$logFC
genes.ids.WT.with.no.Fe <- rownames(WT.with.no.Fe)

activated.genes.WT.with.no.Fe.1 <- genes.ids.WT.with.no.Fe[fold.change.WT.with.no.Fe > 1]
repressed.genes.WT.with.no.Fe.1 <- genes.ids.WT.with.no.Fe[fold.change.WT.with.no.Fe < - 1]

fold.change.pye.with.no.Fe <- pye.with.no.Fe$logFC
genes.ids.pye.with.no.Fe <- rownames(pye.with.no.Fe)

activated.genes.pye.with.no.Fe.1 <- genes.ids.pye.with.no.Fe[fold.change.pye.with.no.Fe > 1]
repressed.genes.pye.with.no.Fe.1 <- genes.ids.pye.with.no.Fe[fold.change.pye.with.no.Fe < - 1]


plot(wt.with.fe,wt.no.fe,pch=19,cex=0.5,col="grey",xlab="WT with Fe",ylab="WT no Fe")

points(wt.with.fe[activated.genes.WT.with.no.Fe.1],
       wt.no.fe[activated.genes.WT.with.no.Fe.1],pch=19,cex=0.5,col="red")

points(wt.with.fe[repressed.genes.WT.with.no.Fe.1],
       wt.no.fe[repressed.genes.WT.with.no.Fe.1],pch=19,cex=0.5,col="blue")

text(wt.with.fe["254550_at"]+0.3,wt.no.fe["254550_at"]+0.3,"IRT1", col="black", cex=0.7)
text(wt.with.fe["252427_at"]+0.3,wt.no.fe["252427_at"]+0.3,"PYE", col="black", cex=0.7)
text(wt.with.fe["257062_at"]+0.3,wt.no.fe["257062_at"]+0.3,"BTS", col="black", cex=0.7)
text(wt.with.fe["251109_at"]+0.3,wt.no.fe["251109_at"]+0.3,"FER1", col="black", cex=0.7)



activated.genes.WT.with.no.Fe.1.table <- aafTableAnn(activated.genes.WT.with.no.Fe.1, 
                                                     "ath1121501.db", aaf.handler())
saveHTML(activated.genes.WT.with.no.Fe.1.table, 
         file="activated_genes_WT_with_no_Fe_1.html")
saveText(activated.genes.WT.with.no.Fe.1.table, 
         file="activated_genes_WT_with_no_Fe_1.txt")

repressed.genes.WT.with.no.Fe.1.table <- aafTableAnn(repressed.genes.WT.with.no.Fe.1, 
                                                     "ath1121501.db", aaf.handler())
saveHTML(repressed.genes.WT.with.no.Fe.1.table, 
         file="repressed_genes_WT_with_no_Fe_1.html")
saveText(repressed.genes.WT.with.no.Fe.1.table, 
         file="repressed_genes_WT_with_no_Fe_1.txt")



activated.genes.pye.with.no.Fe.1.table <- aafTableAnn(activated.genes.pye.with.no.Fe.1, 
                                                      "ath1121501.db", aaf.handler())
saveHTML(activated.genes.pye.with.no.Fe.1.table, 
         file="activated_genes_pye_with_no_Fe_1.html")
saveText(activated.genes.pye.with.no.Fe.1.table, 
         file="activated_genes_pye_with_no_Fe_1.txt")

repressed.genes.pye.with.no.Fe.1.table <- aafTableAnn(repressed.genes.pye.with.no.Fe.1, 
                                                      "ath1121501.db", aaf.handler())
saveHTML(repressed.genes.pye.with.no.Fe.1.table, 
         file="repressed_genes_pye_with_no_Fe_1.html")
saveText(repressed.genes.WT.with.no.Fe.1.table, 
         file="repressed_genes_pye_with_no_Fe_1.txt")

fold.change.WT.with.no.Fe <- WT.with.no.Fe$logFC
p.value.WT.with.no.Fe <- WT.with.no.Fe$adj.P.Val
genes.ids.WT.with.no.Fe <- rownames(WT.with.no.Fe)

activated.genes.WT.with.no.Fe.2 <- genes.ids.WT.with.no.Fe[fold.change.WT.with.no.Fe > 1 & 
                                                             p.value.WT.with.no.Fe < 0.05]

repressed.genes.WT.with.no.Fe.2 <- genes.ids.WT.with.no.Fe[fold.change.WT.with.no.Fe < -1 & 
                                                             p.value.WT.with.no.Fe < 0.05]

names(fold.change.WT.with.no.Fe) <- genes.ids.WT.with.no.Fe
log.p.value.WT.with.no.Fe <- -log10(p.value.WT.with.no.Fe)
names(log.p.value.WT.with.no.Fe)  <- genes.ids.WT.with.no.Fe

plot(fold.change.WT.with.no.Fe,log.p.value.WT.with.no.Fe,
     pch=19,cex=0.5,col="grey",ylab="-log10(p value)",xlab="log2 fold change",xlim=c(-6,6))

points(fold.change.WT.with.no.Fe[activated.genes.WT.with.no.Fe.2],
       log.p.value.WT.with.no.Fe[activated.genes.WT.with.no.Fe.2],
       pch=19,cex=0.5,col="red")

points(fold.change.WT.with.no.Fe[repressed.genes.WT.with.no.Fe.2],
       log.p.value.WT.with.no.Fe[repressed.genes.WT.with.no.Fe.2],
       pch=19,cex=0.5,col="blue")

text(fold.change.WT.with.no.Fe["254550_at"],
     log.p.value.WT.with.no.Fe["254550_at"]+0.3,"IRT1", col="black")

text(fold.change.WT.with.no.Fe["251109_at"]+0.3,
     log.p.value.WT.with.no.Fe["251109_at"],"FER1", col="black")

DEGs.fold.change <- function(contrast.results,
                             gene.number,
                             comparison.number,
                             fold.change.threshold)
{
  ## Extraer con topTable la información de DEGs
  transcriptome.comparison <- topTable(contrast.results, number=gene.number,coef=comparison.number,sort.by="logFC")
  ## Extraer el fold change de cada gen y los nombres de las sondas (genes).
  current.fold.changes <- transcriptome.comparison[["logFC"]]
  gene.ids <- rownames(transcriptome.comparison)
  #gene.ids <- transcriptome.comparison[[1]]
  
  ## Extraer el nombre de los genes cuyo fold change excede el umbral prefijado (genes
  ## activados) o que es menor que menos el umbral prefijado (genes reprimidos).
  current.activated.genes <- gene.ids[current.fold.changes > fold.change.threshold]
  current.repressed.genes <- gene.ids[current.fold.changes < - fold.change.threshold]
  
  ## Devolver una lista cuyas componentes se identifican con los nombres "activated.genes"
  ## y "repressed.genes" y que contienen los correspondientes vectores. 
  return(list(activated.genes = current.activated.genes, repressed.genes = current.repressed.genes))
}

diff.genes.1 <- DEGs.fold.change(contrast.results,gene.number=22810,comparison.number=1,fold.change.threshold=1) 
diff.genes.2 <- DEGs.fold.change(contrast.results,gene.number=22810,comparison.number=2,fold.change.threshold=1) 
diff.genes.3 <- DEGs.fold.change(contrast.results,gene.number=22810,comparison.number=3,fold.change.threshold=1) 
diff.genes.4 <- DEGs.fold.change(contrast.results,gene.number=22810,comparison.number=4,fold.change.threshold=1) 

complete.DEGs <- c(diff.genes.1$activated.genes,
                   diff.genes.1$repressed.genes,
                   diff.genes.2$activated.genes,
                   diff.genes.2$repressed.genes,
                   diff.genes.3$activated.genes,
                   diff.genes.3$repressed.genes,
                   diff.genes.4$activated.genes,
                   diff.genes.4$repressed.genes)

complete.DEGs <- unique(complete.DEGs)


DEG.expression <- mean.expression[complete.DEGs,
                                  c("WT_with_Fe","WT_no_Fe","pye_with_Fe","pye_no_Fe")]

normalized.DEG.expression <- t(scale(t(DEG.expression)))


help(heatmap)
heatmap(normalized.DEG.expression,Colv=FALSE,dendrogram="row",
        labRow=c(""),density.info="none",trace="none",
        col=redgreen(100),margins = c(8,8),cexCol=1.2)
