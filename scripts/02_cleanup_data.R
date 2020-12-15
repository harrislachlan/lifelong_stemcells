# Purpose: remove doublets and remove poor quality cells based on no. of detected genes and %mito ----

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

print(glue("Core dataset contains {length(rownames(seurat_object@meta.data))} cells prior to subsetting"))

# Removedoublets with oligos ----
seurat_object <- subset(seurat_object, subset = Mog > 4 & Aldoc > 4, slot = 'counts', invert = TRUE) %>%
        subset(subset = Mog > 4 & C1qc > 4, slot = 'counts', invert = TRUE) %>%
        subset(subset = Mog > 4 & Cldn5 > 4, slot = 'counts', invert = TRUE) %>%
        subset(subset = Mog > 4 & Stmn2 > 4, slot = 'counts', invert = TRUE) %>%
        # Removedoublets with astro
        subset(subset = Aldoc > 4 & C1qc > 4, slot = 'counts', invert = TRUE) %>%
        subset(subset = Aldoc > 4 & Cldn5 > 4, slot = 'counts', invert = TRUE) %>%
        subset(subset = Aldoc > 4 & Stmn2 > 4, slot = 'counts', invert = TRUE) %>%
        # Removedoublets with microglia
        subset(subset = C1qc > 4 & Cldn5 > 4, slot = 'counts', invert = TRUE)  %>%
        subset(subset = C1qc > 4 & Stmn2 > 4, slot = 'counts', invert = TRUE)  %>%
        subset(subset = Cldn5 > 4 & Stmn2 > 4, slot = 'counts', invert = TRUE)

# Set thresholds for gene and mt content ----
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Chemistry") 
VlnPlot(seurat_object, features = "percentMito", group.by = "GEMgroup", same.y.lims = TRUE, pt.size = 0.1) + ylim(0, 50)
mitoHi <- 10
nGeneLo <- 500
seurat_object <- subset(seurat_object, subset = nFeature_RNA > nGeneLo & percentMito < mitoHi)

# Visually inspect that filtering worked ----
VlnPlot(seurat_object, features = "percentMito", group.by = "GEMgroup", same.y.lims = TRUE, pt.size = 0.1) + ylim(0, 50)
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Chemistry") + ylim(0, 1000)

print(glue("After removal based on doublets/low gene count/high mt content there are {length(rownames(seurat_object@meta.data))} cells remaining)}"))

# End Script 2 ----


