
# ============================================================
# INSTALACION DE PAQUETES - CLASE snRNA-seq
# Seurat + Pseudobulk + DESeq2
# Ejecutar una sola vez ANTES de la clase
# ============================================================

# 1. PAQUETES DE CRAN

paquetes_cran <- c(
  "dplyr",
  "Seurat",
  "svMisc",
  "cowplot",
  "ggplot2",
  "pheatmap",
  "reshape2",
  "gridExtra",
  "RColorBrewer",
  "Matrix",
  "devtools",
  "BiocManager"
)

# Identificar paquetes faltantes
faltantes_cran <- setdiff(
  paquetes_cran,
  rownames(installed.packages())
)

# Instalar solo los faltantes
if (length(faltantes_cran) > 0) {
  install.packages(faltantes_cran, dependencies = TRUE)
}


# 2. PAQUETES DE BIOCONDUCTOR

paquetes_bioc <- c(
  "DESeq2",
  "edgeR",
  "EnhancedVolcano"
)

faltantes_bioc <- setdiff(
  paquetes_bioc,
  rownames(installed.packages())
)

if (length(faltantes_bioc) > 0) {
  BiocManager::install(
    faltantes_bioc,
    update = FALSE,
    ask = FALSE
  )
}


# 3. PAQUETE PRESTO DESDE GITHUB

if (!requireNamespace("presto", quietly = TRUE)) {
  
  devtools::install_github(
    "immunogenomics/presto",
    upgrade = "never"
  )
  
}


# 4. VERIFICACION FINAL

todos_paquetes <- c(
  paquetes_cran,
  paquetes_bioc,
  "presto"
)

estado <- sapply(
  todos_paquetes,
  requireNamespace,
  quietly = TRUE
)

resultado <- data.frame(
  Paquete = todos_paquetes,
  Instalado = estado
)

print(resultado, row.names = FALSE)

if (all(estado)) {
  cat("\nINSTALACION COMPLETA: todos los paquetes disponibles.\n")
} else {
  cat("\nATENCION: faltan paquetes por instalar:\n")
  print(names(estado)[!estado])
}

# Versiones de R y paquetes
sessionInfo()
