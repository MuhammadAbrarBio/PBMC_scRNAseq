# ==========================================================
# Project: PBMC 3K scRNA-seq Analysis using Seurat
# Dataset: 10x Genomics PBMC 3K
# Author: Muhammad Abrar
# Description:
# End-to-end analysis of the PBMC 3K dataset including
# quality control, normalization, clustering,
# differential expression, and cell-type annotation.

# ==========================================================
# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# ==========================================================
# Load Cell Ranger data
pbmc_raw <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19")

# ==========================================================
# Create a Seurat object from the raw count matrix.

# Apply initial filtering based on minimum cell and feature thresholds.
pbmc <- CreateSeuratObject(counts = pbmc_raw, 
                           project = "PBMC3k", 
                           min.cells = 3,
                           min.features = 200
                           )

# ==========================================================
# Quality control (QC)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Calculate the percentage of mitochondrial transcripts
# and store the values in the metadata for quality control.

qc_violin <- VlnPlot(pbmc, 
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     ncol = 3
                     )


ggsave(filename = "output/figures/qc_violin.png", 
       plot = qc_violin,
       height = 4, 
       width = 9,
       dpi = 300
       )

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(pbmc, feature1 = "nFeature_RNA", feature2 = "percent.mt")

plot1 + plot2 + plot3

# Examine relationships between sequencing depth,
# detected genes, and mitochondrial transcript percentage.

pbmc <- subset(pbmc,
               subset = nFeature_RNA > 200 &
                        nFeature_RNA < 2500 &
                        percent.mt < 5
               )

# nFeature_RNA > 200 → remove low-quality cells/empty droplets.
# nFeature_RNA < 2500 → remove potential doublets or multiplets (droplets containing more than one cell).
# percent.mt < 5 → remove cells with high mitochondrial RNA, which are more likely to be stressed or damaged.
# ==========================================================
# Normalization

pbmc <- NormalizeData(pbmc,
                      normalization.method = "LogNormalize",
                      scale.factor = 10000
                      )

# Normalize expression across cells using LogNormalize.
# Each cell is scaled to 10,000 counts prior to log transformation.

# SCTransform() provides an alternative workflow that combines normalization, 
# variable feature selection, and scaling.
# ==========================================================
# Feature Selection (Highly Variable Genes)
pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = "vst",
                             nfeatures = 2000
                             )
# To see the top 10 most differentially/variable expressed genes.
top10_variable_genes <- head(VariableFeatures(pbmc), 10)

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1,
                     points = top10_variable_genes,
                     repel = TRUE
                     )
plot1 + plot2

# ==========================================================
# Scaling the data
# Standardize gene expression prior to PCA.
# Expression values are centered (mean = 0) and scaled (variance = 1).

all_genes <- rownames(pbmc)

pbmc <- ScaleData(pbmc,
                  features = all_genes
                  )

# ==========================================================
# Linear Dimensional Reduction (PCA)
pbmc <- RunPCA(pbmc,
               features = VariableFeatures(pbmc)
               )

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca" )

DimPlot(pbmc, reduction = "pca") + NoLegend()

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

elbow_plot <- ElbowPlot(pbmc, ndims = 25)

ggsave(filename = "output/figures/elbow_plot.png", 
       plot = elbow_plot,
       width = 6,
       height = 4, 
       dpi = 300
)

# ==========================================================
# Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Larger datasets often require a higher clustering resolution.

# ==========================================================
# Non-linear dimensional reduction (UMAP)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)

# ==========================================================
# Differential expression analysis (cluster marker genes)

# Identify marker genes for cluster 2 versus all other clusters.
cluster2_biomarkers <- FindMarkers(pbmc, ident.1 = 2)

# Compare cluster 5 against clusters 0 and 3.
cluster5_biomarkers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5_biomarkers, n = 5)

pbmc_markers <- FindAllMarkers(pbmc, only.pos = TRUE)
head(pbmc_markers, 10)

# cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# Visualize marker expression using raw counts.
VlnPlot(pbmc, 
        features = c("NKG7", "PF4"), 
        layer = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 10 markers (or all markers if less than 10) for each cluster.

top10_markers <- pbmc_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10)  %>%
  ungroup()

DoHeatmap(pbmc, 
          features = top10_markers$gene) + NoLegend() 

# ==========================================================
# Annotate clusters using canonical marker genes.
cluster_new_ids <- c(
  "Naive CD4 T", 
  "CD14+ Mono", 
  "Memory CD4 T", 
  "B", 
  "CD8 T", 
  "FCGR3A+ Mono", 
  "NK", 
  "DC", 
  "Platelet"
  )
names(cluster_new_ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, cluster_new_ids)

DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# ==========================================================
# Save publication-quality UMAP figure
umap_plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))


ggsave(filename = "output/figures/pbmc_umap.png", 
       plot = umap_plot,
       height = 7, 
       width = 12,
       dpi = 300
       )

saveRDS(pbmc, file = "output/pbmc_processed.rds")
# ==========================================================
