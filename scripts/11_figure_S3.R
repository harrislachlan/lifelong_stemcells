# Supp figure 3
# Load packages and define data ----
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


# Scoring tdtomato expression and adding to metadata ----
load("all_data_all_clusters.RData")
rnagenes <- t(GetAssayData(seurat_object.integrated, assay = "RNA", slot = "counts"))
rnagenes <- as.data.frame(rnagenes[ , c("WPRE","tdTomatoLoxP", "bGHpolyA")])
rnagenes$td <- rnagenes$WPRE + rnagenes$bGHpolyA
seurat_object.integrated <- AddMetaData(object = seurat_object.integrated, metadata = rnagenes)

# Create new objects  ----
Idents(seurat_object.integrated) <- seurat_object.integrated$Batch
seurat_object.integrated_Batch2 <- subset(seurat_object.integrated, idents = "2")


# Visualise first sample wherein GFP+ and td+ cells were sequenced separately ----
FeatureScatter(seurat_object.integrated_Batch2, feature1 = "td", 
               feature2 = "tdTomatoLoxP", pt.size = 0.5, group.by = "cells", 
               cols = c("green", "red")) + geom_line(aes(y = 3.5)) + geom_line(aes(x = 1.5)) + ylim(0, 20) + xlim(0, 10)
ggsave(filename = "plots/tdtomato_cutoffs.png", width = 5, height = 3, scale = 1)


# Enforcing of thresholds ----
# V3 of kit - if count of td 4 or more and tdtomatoloxp less than 8 than cell is tdTomato positive. 
seurat_object.integrated@meta.data$'tdthresh' <- ifelse(test = seurat_object.integrated@meta.data$'td' > 3 & seurat_object.integrated@meta.data$'Chemistry' == "V3", yes = seurat_object.integrated@meta.data$'td', no = 0) 
seurat_object.integrated@meta.data$'tdthresh' <- ifelse(test = seurat_object.integrated@meta.data$'tdTomatoLoxP' > 7 & seurat_object.integrated@meta.data$'Chemistry' == "V3", yes = 0, no = seurat_object.integrated@meta.data$'tdthresh') 
# V2 of kit - if count of td 2 or more and tdtomatoloxp less than 4 than cell is tdTomato positive. 
seurat_object.integrated@meta.data$'tdthresh' <- ifelse(test = seurat_object.integrated@meta.data$'td' > 1 & seurat_object.integrated@meta.data$'Chemistry' == "V2", yes = seurat_object.integrated@meta.data$'td', no = seurat_object.integrated@meta.data$'tdthresh') 
seurat_object.integrated@meta.data$'tdthresh' <- ifelse(test = seurat_object.integrated@meta.data$'tdTomatoLoxP' > 3 & seurat_object.integrated@meta.data$'Chemistry' == "V2", yes = 0, no = seurat_object.integrated@meta.data$'tdthresh')

# Binarisation of tdTomato thresholds, td pos = 1 or td neg = 0 
seurat_object.integrated@meta.data$'tdTomato+' <- ifelse(test = seurat_object.integrated@meta.data$'tdthresh' > 0 , yes = 1, no = 0)


# plots ----
DefaultAssay(seurat_object.integrated) <- "RNA"
g1 <- FeaturePlot(seurat_object.integrated, features = "eGFP", split.by = "cells", max.cutoff = 1,order = TRUE)
g2 <- FeaturePlot(seurat_object.integrated, features = "tdTomato+", split.by = "cells", max.cutoff = 1,order = TRUE)
plot_grid(g1, g2, nrow = 2)
G <- arrangeGrob(g1, g2, nrow = 2)
ggsave(filename = "plots/tdtomato_cell_detection.png", G, width = 10, height = 7, scale = 1)

rm(list = ls())