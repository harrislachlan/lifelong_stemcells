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

# Visualisation of all cellular clusters ----
load("all_data_all_clusters.RData")
DimPlot(seurat_object.integrated, reduction = "umap", label = TRUE)
# create list of all numbered cluster names
current.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10","11","12", "13", "14", "15")
new.cluster.ids <- c("OPCs","MOGs","Neuroblasts","NSCs","Microglia","NSC/IPCs","Astrocytes","Endo","activ. Microglia","Endo","COPs","OPCs","multiplets", "OPCs", "Interneuron", "VLMC")
# current cluster names replaced with new names
seurat_object.integrated@active.ident <- plyr::mapvalues(x = seurat_object.integrated@active.ident, from = current.cluster.ids, to = new.cluster.ids)
#add to metadata 
seurat_object.integrated <- AddMetaData(object = seurat_object.integrated, metadata = seurat_object.integrated@active.ident, col.name = "Cluster") 
A <- DimPlot(seurat_object.integrated, reduction = "umap", label = FALSE)
A <- A + theme_minimal()
A
ggsave(filename = "plots/all_cluters_labelled.png", width = 9, height = 5, scale = 0.618)
dev.off()
rm(list = ls())


# Visualise expression to isolate NSCs & IPCs ----
load("all_data_all_clusters.RData")
DefaultAssay(object = seurat_object.integrated) <- "RNA"
seurat_object.integrated <- NormalizeData(seurat_object.integrated)
g1 <- FeaturePlot(object = seurat_object.integrated, features = "Aqp4",
                  reduction = "umap", cols = c("grey","blue4"), order = TRUE, 
                  pt.size = 0.1, combine = TRUE, ncol = 2, label.size = 0.5, max.cutoff = 5) 
g2 <- FeaturePlot(object = seurat_object.integrated, features = "Hopx",
                  reduction = "umap", cols = c("grey","blue4"), order = TRUE, 
                  pt.size = 0.1, combine = TRUE, ncol = 2, label.size = 0.5, max.cutoff = 5)
g3 <- FeaturePlot(object = seurat_object.integrated, features = "Eomes",
                  reduction = "umap", cols = c("grey","blue4"), order = TRUE, 
                  pt.size = 0.1, combine = TRUE, ncol = 2, label.size = 0.5, max.cutoff = 3)
g4 <- FeaturePlot(object = seurat_object.integrated, features = "Neurod1",
                  reduction = "umap", cols = c("grey","blue4"), order = TRUE, 
                  pt.size = 0.1, combine = TRUE, ncol = 2, label.size = 0.5, max.cutoff = 4)
cowplot::plot_grid(g1, g2, g3, g4)
G <- arrangeGrob(g1, g2, g3, g4, nrow = 2)
ggsave(filename = "plots/nsc_ipc_identification.png", G, width = 10, height = 7, scale = 1)
rm(list = ls())


# Visualise two primary clusters (quiescent and cycling cells (NSCs and IPCs)) ----
load("nsc_ipc.RData")
DimPlot(nsc_ipc.integrated, reduction = "umap", label = TRUE) + theme_minimal()
# create list of all numbered cluster names
current.cluster.ids <- c("0", "1", "2","3","4","5","6","7","8","9")
new.cluster.ids <- c("NSC", "NSC", "IPC/NSC","NSC","IPC/NSC","NSC","IPC/NSC","IPC/NSC","NSC","NSC")
nsc_ipc.integrated@active.ident <- plyr::mapvalues(x = nsc_ipc.integrated@active.ident, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(nsc_ipc.integrated, reduction = "umap", label = FALSE) + theme_minimal()
ggsave(filename = "plots/quiescent_cycling_NSPCs.png", width = 8, height = 5, scale = 0.618)
rm(list = ls())

# Visualise classification of NSCs vs IPCs----
load("nsc_ipc.RData")
DefaultAssay(nsc_ipc.integrated) <- "RNA"
nsc_ipc.integrated <- NormalizeData(nsc_ipc.integrated)
g3 <- FeaturePlot(object = nsc_ipc.integrated, features = "Apoe",
                  reduction = "umap", cols = c("grey","blue4"), order = FALSE, coord.fixed = TRUE,
                  pt.size = 0.5, label.size = 0.5, max.cutoff = 4) + theme_void()
g5 <- FeaturePlot(object = nsc_ipc.integrated, features = "Hopx",
                  reduction = "umap", cols = c("grey","blue4"), order = FALSE, coord.fixed = TRUE,
                  pt.size = 0.5, label.size = 0.5, max.cutoff = 1) + theme_void()
g16 <- FeaturePlot(object = nsc_ipc.integrated, features = "Neurod1",
                   reduction = "umap", cols = c("grey","blue4"), order = FALSE, coord.fixed = TRUE,
                   pt.size = 0.5,  label.size = 0.5, max.cutoff = 1) + theme_void()
g22 <- FeaturePlot(object = nsc_ipc.integrated, features = "Dcx",
                   reduction = "umap", cols = c("grey","blue4"), order = FALSE, coord.fixed = TRUE,
                   pt.size = 0.5, label.size = 0.5, max.cutoff = 1) + theme_void()
cowplot::plot_grid(g3, g5, g16, g22)
g <- arrangeGrob(g3, g5, g16, g22, nrow = 2)
ggsave(file="plots/nsc_classification.png", g, width = 6, scale = 1) 
rm(list = ls())
