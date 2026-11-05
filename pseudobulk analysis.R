library(DESeq2)
library(Seurat)
library(Matrix)
library(ggplot2)

#importar el objeto seurat con el analisis single cell

single_cell_data <- readRDS("results/single_cell_object.rds")

# Seleccionar la población celular de interés
Drd2 <- subset( single_cell_data, subset = CellType == "Drd2-MSN" )

# Ver cuántos núcleos tenemos por muestra
table(Drd2$Sample)

# Seleccionar el ensayo RNA
DefaultAssay(Drd2) <- "RNA"

# Sumar los conteos UMI de los núcleos de cada muestra
pseudobulk <- AggregateExpression( Drd2,  assays = "RNA",
                                   group.by = "Sample",
                                   return.seurat = FALSE)

# Extraer la matriz de conteos
counts_pb <- pseudobulk$RNA

count.matrix <- as.matrix(counts_pb)

# Examinar dimensiones y nombres de las muestras
dim(count.matrix)
colnames(count.matrix)

# creamos las etiquetas del analisis
labels <- ifelse(grepl("Cocaine", colnames(count.matrix)), "Cocaine","Saline" )

labels <- factor(labels, levels = c("Saline", "Cocaine"))

# Crear una tabla con la información de las muestras
sample.table <- data.frame(row.names = colnames(counts_pb),
                           group = labels )

sample.table

## Filter low-expression genes with edgeR
design <- model.matrix(~ labels)

keep <- edgeR::filterByExpr(edgeR::DGEList(counts = count.matrix), 
                            design = design)

count.matrix <- count.matrix[keep, ]

dim(count.matrix)

##################################################################
FDR <- 0.1 #setear

# Crear el objeto DESeq2
dds <- DESeqDataSetFromMatrix(countData = count.matrix,
                                colData = sample.table,
                                 design = ~ group )

#ANALISIS ESTADISTICO CON DESeq2

# add analysis to dds object
dds <- DESeq(dds)

#results
res.DESeq2 <- results(dds)

#vemos las primeras 5 filas de results
head(res.DESeq2)

#resumen estandar de los resultados que arroja DESeq2
summary(res.DESeq2)

#conteos de resultados
table(res.DESeq2$padj < FDR)

table(res.DESeq2$padj < FDR & res.DESeq2$log2FoldChange< -1) # downregulated
table(res.DESeq2$padj < FDR & res.DESeq2$log2FoldChange> 1) # upregulated

#seleccionamos genes UP y DOWN
UP = res.DESeq2[which(res.DESeq2$padj < 0.1 & res.DESeq2$log2FoldChange > 1),] %>% na.omit()
DOWN = res.DESeq2[which(res.DESeq2$padj < 0.1 & res.DESeq2$log2FoldChange< -1),] %>% na.omit()

#Unimos UP y DOWN en una sola tabla
DEGs.DESeq2 = rbind(UP,DOWN)
DEGs.DESeq2 <- as.data.frame(DEGs.DESeq2)

head(DEGs.DESeq2)

#search genes of interest
#res.DESeq2[which(res.DESeq2$symbol=="ACTA2"),]

#export results
write.csv2(DEGs.DESeq2, "results/DEGs Drd2 MSNs coc vs saline DESeq2.csv")

#-------------------------------------------------------------------------------
#visualizaciones
#-------------------------------------------------------------------------------
#aplicamos transformaciones para mejorar la visualizacion
#covert to regularized logarithm
rld <- rlog(dds)

#convert to variance stabilizing transformation
vst <- vst(dds)

#A PCA plot and a heatmap of the top genes: we need to use the rld object

#Transformación logarítmica regularizada
plotPCA(rld, intgroup="group")

#Transformación estabilizadora de varianza
plotPCA(vst, intgroup="group")

library(EnhancedVolcano)

#hacemos un volcano plot de los resultados
EnhancedVolcano(
  res.DESeq2,
  lab = rownames(res.DESeq2),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = FDR,
  FCcutoff = 1,
  title = "DRD2-MSN cocaine vs saline"
)

#heatmap de DEGs
library(pheatmap)

indexDEGs <- which(rownames(res.DESeq2) %in% rownames(DEGs.DESeq2))

#heatmap for DEGs
mat.DEseq2 <- assay(rld)[indexDEGs,] #opcion 1

#mat.DEseq2 <- assay(vst)[indexDEGs,] #opcion 2

mat.DEseq2 <- t(scale(t(mat.DEseq2)))  # z score

head(mat.DEseq2)

pheatmap(mat.DEseq2, 
         fontsize = 7.5)

#guardamos el grafico en pdf
pheatmap(mat.DEseq2, 
         fontsize = 7.5, #tamanio de la fuente 
         filename = "plots/heatmap_DESeq2.pdf")

graphics.off()
