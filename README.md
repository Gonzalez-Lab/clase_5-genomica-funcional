# Single-nucleus RNA-seq of rat nucleus accumbens | Seurat v5 + pseudobulk DESeq2

**Reproducible analysis workflow · R / Bioconductor · Teaching repository**

## Overview

This repository implements a two-stage analysis of single-nucleus RNA-sequencing (snRNA-seq) data from rat nucleus accumbens, using four experimental libraries spanning sex and acute cocaine exposure. The workflow first resolves transcriptionally defined nuclear populations through Seurat integration and marker-based annotation, then estimates cocaine-associated expression changes within **Drd2 medium spiny neurons (Drd2-MSN)** using sample-level pseudobulk counts and DESeq2.

**Analysis scope:** Cell Ranger count matrices → nuclear QC → log-normalization → highly variable features → anchor-based integration → PCA/UMAP/SNN clustering → cell-type annotation → cell-type composition and marker visualization → Drd2-MSN pseudobulk → differential expression and diagnostic plots.

The repository uses the **classical Seurat anchor-based integration workflow** (`FindIntegrationAnchors`/`IntegrateData`), not Seurat v5 `IntegrateLayers`. Pseudotime/Monocle3 is **not included**.

## Experimental design

| Sample | Sex | Exposure | Input directory |
|---|---|---|---|
| Female_Saline | Female | Saline | `data/Female_Saline/` |
| Female_Cocaine | Female | Cocaine | `data/Female_Cocaine/` |
| Male_Saline | Male | Saline | `data/Male_Saline/` |
| Male_Cocaine | Male | Cocaine | `data/Male_Cocaine/` |

**Unit of observation:** individual nuclei for clustering/visualization; **library/sample** for pseudobulk differential expression. Each sex × treatment combination is represented by one library. The `~ group` DESeq2 model is intentionally simplified for instruction and does not adjust for sex. Differential-expression statistics are therefore **exploratory**, not confirmatory evidence of a reproducible treatment effect.

## Repository layout

```text
.
├── README.md
├── 00_instalar_paquetes.R           # Pre-class setup, if supplied
├── scripts/
│   ├── script single cell.R
│   └── pseudobulk analysis.R
├── data/
│   ├── Female_Saline/
│   ├── Female_Cocaine/
│   ├── Male_Saline/
│   └── Male_Cocaine/
├── results/
└── plots/
```

Each `data/<sample>/` directory must contain a Cell Ranger **filtered feature-barcode matrix** readable by `Read10X()`: `matrix.mtx[.gz]`, `features.tsv[.gz]` (or equivalent gene-feature file), and `barcodes.tsv[.gz]`. The data are not embedded in this README. The actual dataset accession and distribution URL must be verified and added before public release.

## Environment and dependencies

- R and RStudio; Seurat v5-compatible installation.
- CRAN: `Seurat`, `Matrix`, `dplyr`, `ggplot2`, `cowplot`, `pheatmap`, `reshape2`, `gridExtra`, `RColorBrewer`, `svMisc`.
- Bioconductor: `DESeq2`, `edgeR`, `EnhancedVolcano`.
- Optional performance dependency: `presto` from `immunogenomics/presto`.

Install packages using `00_instalar_paquetes.R` if provided. Record `sessionInfo()` for reproducibility. Package versions are not pinned in the supplied scripts; exact reproduction across installations is therefore not guaranteed.

**Working directory:** run both scripts from the repository root. Before execution:

```r
getwd()
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)
```

## Workflow 01 — Seurat preprocessing, integration and annotation

**Entry point:** `scripts/script single cell.R`  
**Primary output in memory:** `All_Groups_log` (Seurat object).

| Stage | Implementation | Parameters / output |
|---|---|---|
| Import | `Read10X`; `CreateSeuratObject` | `min.cells = 1`, `min.features = 1` |
| Nuclear QC | `PercentageFeatureSet`; `subset` | Mitochondrial pattern `^Mt-`; `nFeature_RNA > 200`; `percent_mito < 5` |
| Normalization | `NormalizeData` | `LogNormalize`, `scale.factor = 10000` |
| Feature selection | `FindVariableFeatures` | `vst`, `nfeatures = 2000` per sample |
| Preliminary dimensionality reduction | `ScaleData`; `RunPCA`; `ElbowPlot` | Up to 30 PCs per sample, inspected to choose dimensions |
| Integration | `FindIntegrationAnchors`; `IntegrateData` | `dims = 1:17`; default anchor reduction (CCA) |
| Joint embedding | `ScaleData`; `RunPCA`; `RunUMAP` | 17 PCs; UMAP from PCs 1–17 |
| Graph clustering | `FindNeighbors`; `FindClusters` | PCs 1–17; resolution `0.2` |
| Marker detection | `JoinLayers` on `RNA`; `FindAllMarkers` | `only.pos = TRUE`, `min.pct = 0.25`, `logfc.threshold = 0.25` |
| Annotation | `DotPlot`; `RenameIdents` | Marker-guided identities in `CellType` |
| Exploration | `DimPlot`; `FeaturePlot`; `VlnPlot` | Sample mixing, marker expression, cell-type composition |

### Seurat assay semantics

- **`integrated`**: corrected representation used for joint PCA, UMAP and graph clustering. Integration can attenuate biological differences associated with exposure.
- **`RNA` / `counts`**: original UMI counts; used for pseudobulk aggregation and count-based modeling.
- **`RNA` / `data`**: log-normalized expression; used for gene-level visualization and marker analysis.
- **`JoinLayers(assay = "RNA")`**: joins corresponding sample-specific layers into common `counts`/`data` layers. It does **not** re-integrate, sum nuclei, or overwrite the `integrated` assay.

Cluster indices are algorithmic identifiers, **not invariant biological labels**. The script contains provisional annotations (`undefined`, `pending`) that should be validated against markers before biological interpretation. The `Drd2-MSN` subset used downstream is the population assigned that exact label, **not automatically every Drd2-expressing subtype** (e.g. `Drd2-MSN-2`).

### Hand-off to workflow 02

The second script requires `All_Groups_log` with metadata fields `CellType`, `Sample`, `Sex`, and `Stim`. It can be run in the **same R session** after workflow 01. To persist the object explicitly:

```r
saveRDS(All_Groups_log, "results/All_Groups_log.rds")
# In a new R session:
All_Groups_log <- readRDS("results/All_Groups_log.rds")
```

These `saveRDS`/`readRDS` commands are recommended additions, **not automatic steps in the supplied scripts**.

## Workflow 02 — Drd2-MSN pseudobulk differential expression

**Entry point:** `scripts/pseudobulk analysis.R`  
**Input:** annotated Seurat object `All_Groups_log` from workflow 01.

1. Subset nuclei using `CellType == "Drd2-MSN"`; tabulate nuclei by `Sample`.
2. Set `DefaultAssay(Drd2) <- "RNA"`.
3. Generate pseudobulk counts with `AggregateExpression(assays = "RNA", group.by = "Sample")`: UMI counts are **summed across nuclei within each sample**, producing a genes × samples matrix.
4. Construct the binary treatment factor with `Saline` as reference.
5. Filter low-expression genes with `edgeR::filterByExpr`.
6. Fit `DESeqDataSetFromMatrix(..., design = ~ group)` and `DESeq()`.
7. Extract Wald-test results (`results(dds)`), summarize adjusted p values, and classify genes using **FDR < 0.1** and **absolute log2 fold change > 1**.
8. Visualize with `rlog`, variance-stabilizing transformation, `plotPCA`, `EnhancedVolcano`, and row-z-scored `pheatmap`.

**Contrast interpretation:** with `Saline` as the factor reference, positive `log2FoldChange` corresponds to higher expression in Cocaine and negative values to higher expression in Saline. For clarity and reproducibility, `results(dds, contrast = c("group", "Cocaine", "Saline"))` is an explicit equivalent to the intended comparison.

**Statistical limitation:** the four libraries provide only two samples per treatment and one per sex × treatment combination. Treat nuclei as subsamples, not independent biological replicates. The analysis does not model sex or sex-by-treatment interaction, and its FDR estimates should not be interpreted as robust confirmation of a treatment effect. This dataset is appropriate for demonstrating the computational workflow and generating hypotheses.

## Outputs and interpretation

| Artifact | Generated by | Interpretation |
|---|---|---|
| Cluster marker table (`results/genes markers por cluster.csv`) | `FindAllMarkers` | Positive cluster-enriched markers; not treatment DEGs |
| Annotated UMAP and sample-mixing plots | `DimPlot` | Joint embedding and sample representation |
| Cell-type composition bar plot | `ggplot2` | Fractions of **recovered nuclei**, not direct tissue abundances |
| `results/DEGs Drd2 MSNs coc vs saline DESeq2.csv` | DESeq2 script | Genes passing script-defined FDR and log2FC thresholds |
| `plots/heatmap_DESeq2.pdf` | `pheatmap` | Relative expression patterns for selected genes, if available |
| PCA / volcano plots | `plotPCA`, `EnhancedVolcano` | Sample-level structure and DE summary; displayed interactively unless separately saved |

**Note:** the supplied DESeq2 script assumes at least some selected DEGs when constructing its heatmap. If none pass the thresholds, skip the heatmap or select a clearly labeled descriptive set of genes; do not invent significant results.

## Reproducibility and quality checks

Before interpreting outputs, verify:

```r
# Input metadata and cluster representation
stopifnot(all(c("Sample", "Stim", "Sex", "CellType") %in%
              colnames(All_Groups_log@meta.data)))
table(All_Groups_log$Sample, All_Groups_log$CellType)

# Pseudobulk matrix / metadata correspondence (after constructing them)
stopifnot(identical(colnames(count.matrix), rownames(sample.table)))
stopifnot(all(count.matrix >= 0))
```

Check sufficient nuclei per sample in the selected subtype; inspect missing cell types, genes absent from the reference annotation, and whether integration preserves expected marker patterns. For scientific reuse, document R/package versions, input matrix provenance, QC distributions, and any changes to integration dimensions, clustering resolution, or annotation.

## Data provenance and citation

**Dataset:** rat nucleus accumbens snRNA-seq, cocaine vs saline, four libraries (female/male × treatment).  
**Public accession / primary publication:** *to be verified against the exact matrices used in class before release*.  
**Software:** Seurat, DESeq2, edgeR, EnhancedVolcano and associated R packages; cite the relevant software publications in derivative scientific work.

## Intended use

Hands-on postgraduate instruction in functional genomics and single-cell transcriptomics. The scripts are designed to expose individual processing steps and analytical assumptions, not to serve as a fully validated production pipeline.
