# ==========================================================
# PBMC 3K scRNA-seq Analysis using Seurat
# Dataset: 10x Genomics PBMC 3K
# Author: Muhammad Abrar


# ==========================================================
# Load libraries
# Load required packages
library(dplyr)
library(Seurat)
library(patchwork)


# ==========================================================
# Load Cell Ranger data
pbmc_raw <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19")


# ==========================================================
# Create a Seurat object from the raw count matrix
pbmc <- CreateSeuratObject(counts = pbmc_raw, 
                           project = "PBMC3k", 
                           min.cells = 3,
                           min.features = 200)


# ==========================================================
# Quality Control (QC)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


head(pbmc@meta.data, 6)

VlnPlot(pbmc, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 2500 & 
                        percent.mt < 5)


# ==========================================================
# Normalization
pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 1e4)



# ==========================================================
# Feature Selection (Highly Variable Genes)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)

plot1 <- VariableFeaturePlot(pbmc)
plot2 <-LabelPoints(
  plot = plot1, 
  points = top10, 
  repel = TRUE
  )

plot1 + plot2


# ==========================================================
# Scaling the data
# all_genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, 
                  features = VariableFeatures(pbmc), 
                  vars.to.regress = "percent.mt")

# pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")


# ==========================================================
# Linear Dimensional Reduction (PCA)


# ==========================================================
# Clustering


# ==========================================================
# UMAP Visualization


