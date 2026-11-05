library("dplyr")
library("Seurat")
library("svMisc")
library("cowplot")
library("ggplot2")
library("pheatmap")
library("reshape2")
library("gridExtra")
library("RColorBrewer")
library("dplyr")
library("Matrix")

#CellRanger output will be loaded into the local environment using the Read10X function

Fem_Sal_data  <- Read10X("data/Female_Saline")
Fem_Coc_data  <- Read10X("data/Female_Cocaine")
Male_Sal_data <- Read10X("data/Male_Saline")
Male_Coc_data <- Read10X("data/Male_Cocaine")

#Create Seurat Object function

Fem_Sal  <- CreateSeuratObject(counts = Fem_Sal_data,  min.cells = 1  ,min.features = 1) 
Fem_Coc  <- CreateSeuratObject(counts = Fem_Coc_data,  min.cells = 1,  min.features = 1) 
Male_Sal <- CreateSeuratObject(counts = Male_Sal_data, min.cells = 1,  min.features = 1)  
Male_Coc <- CreateSeuratObject(counts = Male_Coc_data, min.cells = 1,  min.features = 1) 

dim(Fem_Sal)          # Genes × núcleos
Layers(Fem_Sal)       # solo tenemos layer counts

# Extraer la matriz de conteos con GetAssayData()
Fem_Sal_counts <- GetAssayData( Fem_Sal, assay = "RNA",layer = "counts" ) 

# Ver los primeros 5 genes y 5 núcleos
as.matrix(Fem_Sal_counts[1:5, 1:5])

#removing cells with a high percentage of reads mapping to the mitochondrial genome
Fem_Sal  <- PercentageFeatureSet(Fem_Sal, pattern = "^Mt-", col.name = "percent_mito")
Fem_Coc  <- PercentageFeatureSet(Fem_Coc, pattern = "^Mt-", col.name = "percent_mito")
Male_Sal <- PercentageFeatureSet(Male_Sal, pattern = "^Mt-", col.name = "percent_mito")
Male_Coc <- PercentageFeatureSet(Male_Coc, pattern = "^Mt-", col.name = "percent_mito")

#quality control, aplicamos filtros
Fem_Sal  <- subset(x = Fem_Sal, subset  =  nFeature_RNA > 200 & percent_mito < 5)
Fem_Coc  <- subset(x = Fem_Coc, subset  =  nFeature_RNA > 200 & percent_mito < 5)  
Male_Sal <- subset(x = Male_Sal, subset =  nFeature_RNA > 200 & percent_mito < 5)  
Male_Coc <- subset(x = Male_Coc, subset =  nFeature_RNA > 200 & percent_mito < 5)

#violin plots
VlnPlot(Fem_Sal, features = c("nFeature_RNA","percent_mito"), ncol = 2,pt.size = 0)
VlnPlot(Fem_Coc, features = c("nFeature_RNA","percent_mito"), ncol = 2,pt.size = 0)
VlnPlot(Male_Sal, features = c("nFeature_RNA","percent_mito"), ncol = 2,pt.size = 0)
VlnPlot(Male_Coc, features = c("nFeature_RNA","percent_mito"), ncol = 2,pt.size = 0)

VlnPlot(Fem_Sal, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"),
                               ncol = 3,   pt.size = 0 )

#Los núcleos con mayor número de UMI suelen presentar más genes detectados.
FeatureScatter( Fem_Sal,  feature1 = "nCount_RNA", feature2 = "nFeature_RNA" )

### Data Normalization: Log-normalization method with a scaling factor of 1000
# se crea la layer data
Fem_Sal  <- NormalizeData(Fem_Sal, normalization.method = "LogNormalize", scale.factor = 10000)
Fem_Coc  <- NormalizeData(Fem_Coc, normalization.method = "LogNormalize", scale.factor = 10000)
Male_Sal <- NormalizeData(Male_Sal, normalization.method = "LogNormalize", scale.factor = 10000)
Male_Coc <- NormalizeData(Male_Coc, normalization.method = "LogNormalize", scale.factor = 10000)

# Extraer la matriz normalizada
Fem_Sal_expr <- GetAssayData( Fem_Sal,
                              assay = "RNA",
                              layer = "data"  ) #obtenemos counts normalizadas 

as.matrix(Fem_Sal_expr[1:5, 1:5])

# 2000 of the most variable features within each dataset are identified
Fem_Sal  <- FindVariableFeatures(Fem_Sal, selection.method = "vst", nfeatures = 2000)
Fem_Coc  <- FindVariableFeatures(Fem_Coc, selection.method = "vst", nfeatures = 2000)
Male_Sal <- FindVariableFeatures(Male_Sal, selection.method = "vst", nfeatures = 2000)
Male_Coc <- FindVariableFeatures(Male_Coc, selection.method = "vst", nfeatures = 2000)

#vector con los 2000 genes mas variables de Fem_Sal
Fem_Sal_features <- VariableFeatures(Fem_Sal) 

# Graficar variabilidad de los genes
VariableFeaturePlot(Fem_Sal)

LabelPoints( VariableFeaturePlot(Fem_Sal), 
             points = Fem_Sal_features[1:10], 
             repel = TRUE )

## Data Integration

#Add treatment information
Fem_Sal$Stim  <- "Saline"
Fem_Coc$Stim  <- "Cocaine"
Male_Sal$Stim <- "Saline"
Male_Coc$Stim <- "Cocaine"

#Add sex information
Fem_Sal$Sex   <- "Female"
Fem_Coc$Sex   <- "Female"
Male_Sal$Sex  <- "Male"
Male_Coc$Sex  <- "Male"

# Definimos numero de dimensiones para integracion

# PCA exploratoria en Female Saline
Fem_Sal <- ScaleData(Fem_Sal, verbose = FALSE)
Fem_Sal <- RunPCA(Fem_Sal, npcs = 30, verbose = FALSE)
ElbowPlot(Fem_Sal, ndims = 30)

# Female Cocaine
Fem_Coc <- ScaleData(Fem_Coc, verbose = FALSE)
Fem_Coc <- RunPCA(Fem_Coc, npcs = 30, verbose = FALSE)
ElbowPlot(Fem_Coc, ndims = 30)

# Male Saline
Male_Sal <- ScaleData(Male_Sal, verbose = FALSE)
Male_Sal <- RunPCA(Male_Sal, npcs = 30, verbose = FALSE)
ElbowPlot(Male_Sal, ndims = 30)

# Male Cocaine
Male_Coc <- ScaleData(Male_Coc, verbose = FALSE)
Male_Coc <- RunPCA(Male_Coc, npcs = 30, verbose = FALSE)
ElbowPlot(Male_Coc, ndims = 30)

# INTEGRACION

# Buscamos las anclas de integracion
All_Groups_log <- FindIntegrationAnchors(object.list = list(Fem_Sal, Fem_Coc,
                                                            Male_Sal, Male_Coc), 
                                         dims = 1:17  )

# integramos los 4 datasets
All_Groups_log  <- IntegrateData(anchorset = All_Groups_log, dims = 1:17)

##################
## standard Seurat workflow

All_Groups_log <- ScaleData(All_Groups_log,verbose = FALSE)

All_Groups_log <- RunPCA(All_Groups_log,npcs = 17 ,verbose = FALSE)

All_Groups_log <- RunUMAP(All_Groups_log, reduction = "pca", dims = 1:17)

All_Groups_log <- FindNeighbors(All_Groups_log, reduction = "pca", dims = 1:17)

All_Groups_log <- FindClusters(All_Groups_log, resolution = 0.2)

# clusters identificados por Seurat
table(All_Groups_log$seurat_clusters)

#UMAP visualization
DimPlot(All_Groups_log, reduction = "umap", label = TRUE) + NoLegend()

# inspeccionamos los ensayos que tenemos en el objeto integrado
# RNA=UMIs normalized , integrated=expresion corregida para integracion
Assays(All_Groups_log)

# Ver las capas actuales del assay RNA
Layers(All_Groups_log[["RNA"]])

# Ver las capas del ensayo integrated
Layers(All_Groups_log[["integrated"]])

# Ver qué ensayo está seleccionado por default
DefaultAssay(All_Groups_log)

####################################################
# Colorear por muestra de origen para chequear integracion

# Crear una nueva columna con la identidad de cada muestra
All_Groups_log$Sample <- paste(All_Groups_log$Sex, All_Groups_log$Stim,sep = "_")

# Comprobar las cuatro muestras
table(All_Groups_log$Sample)

#verificamos que los clusters tengan celulas de los distintos datasets
DimPlot( All_Groups_log, reduction = "umap", group.by = "Sample") +
  ggtitle("Distribución de las cuatro muestras")

##########################################################
#IDENTIFICACION DE LA IDENTIDAD CELULAR DE LOS CLUSTERS 

# Unir las capas para buscar marcadores por cluster
All_Groups_log <- JoinLayers(object = All_Groups_log, 
                              assay = "RNA")

# Comprobar que se unieron
Layers(All_Groups_log[["RNA"]])

# Utilizar la expresión normalizada original por default
DefaultAssay(All_Groups_log) <- "RNA"

# Recuperar las identidades numéricas originales
Idents(All_Groups_log) <- "seurat_clusters"

# Marcadores de tipos celulares
marcadores <- c(
  "Drd1", "Drd2", "Ppp1r1b",
  "Tac1", "Penk",
  "Aqp4", "Slc1a3",
  "Mbp", "Plp1",
  "P2ry12", "Cx3cr1",
  "Pdgfra", "Cspg4",
  "Pvalb", "Sst",
  "Rgs5" )

# Conservar los genes presentes en el dataset
marcadores <- intersect(marcadores, rownames(All_Groups_log[["RNA"]]))

# visualizamos expresion de genes marcadores de identidad celular
DotPlot(All_Groups_log, features = marcadores) + RotatedAxis()

# Buscar genes enriquecidos en cada cluster
markers <- FindAllMarkers(All_Groups_log, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25)

# Mostrar los primeros resultados
head(markers)

write.csv2(markers, "results/genes markers por cluster.csv")

# crear tabla con top10 markers por cluster
top10 <- markers %>% group_by(cluster) %>% slice_max(
          order_by = avg_log2FC,  n = 10, with_ties = FALSE ) %>%
          ungroup() %>% select(cluster, gene, avg_log2FC, pct.1, pct.2)

print(top10, n = 160)

# miramos un cluster el particular
ver_cluster <- 7

markers %>% filter(cluster == ver_cluster) %>% arrange(desc(avg_log2FC)) %>%
            select(gene, avg_log2FC, pct.1, pct.2) %>% head(30)


# Now that clusters have been identified, the arbitrary numbers can be replaced 
# with cell type identities

All_Groups_log <- RenameIdents( All_Groups_log,
                                "0"  = "Olig-1",
                                "1"  = "Drd1-MSN",
                                "2"  = "Drd2-MSN",
                                "3"  = "Olig-2",
                                "4"  = "Neuronal-undefined",
                                "5"  = "Astrocytes",
                                "6"  = "MSN-undefined",
                                "7"  = "Polydendrocyte",
                                "8"  = "Microglia",
                                "9"  = "Cluster-9-pending",
                                "10" = "MSN-10-pending",
                                "11" = "Sst-Interneuron",
                                "12" = "Drd2-MSN-2",
                                "13" = "Mural",
                                "14" = "Olig-3",
                                "15" = "Cluster-15-pending")

All_Groups_log$CellType <- Idents(All_Groups_log)

Idents(All_Groups_log) <- All_Groups_log$CellType

#UMAP visualization
DimPlot(object = All_Groups_log,reduction = "umap",label = TRUE) + NoLegend()

# guardar el objeto para analisis posteriores 
saveRDS( All_Groups_log, "results/single_cell_object.rds")

# Contar núcleos por muestra y tipo celular
tabla <- table( All_Groups_log$Sample, All_Groups_log$CellType )

# Calcular proporciones por muestra
prop <- prop.table(tabla, margin = 1) * 100

library(ggplot2)

# Convertir la tabla de proporciones a formato largo
prop.df <- as.data.frame(prop)

colnames(prop.df) <- c("Sample", "CellType", "Percentage")

# Graficar
ggplot(prop.df, aes(x = Sample, y = Percentage, fill = CellType)) +
  geom_col(width = 0.75) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    title = "Composición celular por muestra",
    x = "Muestra",
    y = "Porcentaje de núcleos",
    fill = "Tipo celular"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.text = element_text(size = 9)
  )

#vemos un gen de interes en todos los clusters
FeaturePlot(All_Groups_log, features = c("Fosb", "Drd1","Drd2"), pt.size =1.5, max.cutoff = 1, 
            cols = c("gray90", "red")) + coord_fixed(ratio = 1)

# Ver la expresión normalizada original de Fos
FeaturePlot(All_Groups_log,  features = "Fosb", split.by = "Stim" )

# Distribución de expresión por población celular
VlnPlot( All_Groups_log, features = "Fosb", pt.size = 0)

# Seleccionar la población celular de interés
Drd2 <- subset( All_Groups_log, subset = CellType == "Drd2-MSN" )

VlnPlot( Drd2, features = "Adora2a",  split.by = "Stim", pt.size = 0)

FeaturePlot(Drd2, features = c("Fosb","Adora2a"), pt.size =1.5, max.cutoff = 1, 
                  cols = c("gray90", "red")) + coord_fixed(ratio = 1)

FeaturePlot(Drd2, features = c("Fosb"), split.by = "Stim",
            pt.size = 1.5, 
            max.cutoff = 1, 
            cols = c("gray90", "red")) + coord_fixed(ratio = 1)
