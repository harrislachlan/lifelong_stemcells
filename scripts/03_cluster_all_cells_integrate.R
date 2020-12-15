# Purpose: use SCTransform, integration anchors, UMAP to visualise data ----

# Required packages ----
library(reticulate)
library(umap)
library(Matrix)
library(lubridate)
library(glue)
library(RColorBrewer)
library(gridExtra)
library(gdata)
library(slingshot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(gridExtra)
library(cowplot)
library(SummarizedExperiment)
library(Seurat)
library(tidyverse)

# Split object by GEMgroup and SCTransform ----
seurat_object.list <- SplitObject(seurat_object, split.by = "GEMgroup")
for (i in 1:length(seurat_object.list)) {
  seurat_object.list[[i]] <- SCTransform(seurat_object.list[[i]], 
                  verbose = TRUE, return.only.var.genes = FALSE)
}

# Select features for to identify anchors ----
seurat_object.features <- 
  SelectIntegrationFeatures(object.list = seurat_object.list, nfeatures = 5000)
seurat_object.list <- 
  PrepSCTIntegration(object.list = seurat_object.list, anchor.features = seurat_object.features, 
                     verbose = TRUE)

# Identify anchors between datasets to integrate ----
seurat_object.anchors <- FindIntegrationAnchors(object.list = seurat_object.list, normalization.method = "SCT", 
                                               anchor.features = seurat_object.features, verbose = TRUE)
seurat_object.integrated <- IntegrateData(anchorset = seurat_object.anchors, normalization.method = "SCT", 
                                         verbose = TRUE)

# Visualisation/clustering ----
seurat_object.integrated <- RunPCA(seurat_object.integrated, npcs = 40, verbose = TRUE) %>% 
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.3)
ElbowPlot(seurat_object.integrated, ndims = 30, reduction = "pca")
DimPlot(seurat_object.integrated, reduction = "umap", label = TRUE)

# Save data ----
gdata::keep(seurat_object.integrated, sure = TRUE)
save.image(file="all_data_all_clusters1.RData")
rm(list = ls())

# End Script 3 ----

