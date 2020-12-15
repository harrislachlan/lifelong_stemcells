# Purpose: isolate NSCs/IPCs, recluster, isolate NSCs, recluster ——— 

# Required packages
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

# Load data —---
load("all_data_all_clusters.Rdata")

# Create list of all numbered cluster names, replace old names for the clusters of NSCs and IPCs, check labelling ———
current.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10","11","12", "13", "14", "15")
new.cluster.ids <- c("0","1","2","NSC","4","IPC","6","7","8","9","10","11","12","13", "14", "15")
seurat_object.integrated@active.ident <- plyr::mapvalues(x = seurat_object.integrated@active.ident, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(seurat_object.integrated, reduction = "umap", label = TRUE)

# Remove all cells that are not in NSC and IPC cluster —---
nsc_ipc <- subset(x = seurat_object.integrated, idents = c("NSC", "IPC"))
DimPlot(nsc_ipc, reduction = "umap", label = TRUE)
print(glue("After keeping only the NSC and IPC clusters there are {length(rownames(seurat_object.integrated@meta.data))} cells remaining"))
rm(seurat_object.integrated)

# Split object by Batch and SCTransform. 
#Split by batch as too few cells from GEMgroup 2 (tdTomato+ only cells) are present to split by GEMgroup ———-
nsc_ipc.list <- SplitObject(nsc_ipc, split.by = "Batch")
for (i in 1:length(nsc_ipc.list)) {
  nsc_ipc.list[[i]] <- SCTransform(nsc_ipc.list[[i]], 
return.only.var.genes = FALSE, verbose = TRUE)
}

# Select features for identification of anchors and downstream integration ———-
nsc_ipc.features <- SelectIntegrationFeatures(object.list = nsc_ipc.list, nfeatures = 5000)
nsc_ipc.list <- PrepSCTIntegration(object.list = nsc_ipc.list, anchor.features = nsc_ipc.features, 
                                         verbose = TRUE)

# Identify anchors —---
nsc_ipc.anchors <- FindIntegrationAnchors(object.list = nsc_ipc.list, normalization.method = "SCT", 
                                                anchor.features = nsc_ipc.features, verbose = TRUE)
nsc_ipc.integrated <- IntegrateData(anchorset = nsc_ipc.anchors, normalization.method = "SCT", 
                                          verbose = TRUE)

# Visualisation/clustering ———-
nsc_ipc.integrated <- RunPCA(nsc_ipc.integrated, verbose = TRUE) %>%
        RunUMAP(dims = 1:25) %>%
        FindNeighbors(reduction = "pca", dims = 1:25) %>%
        FindClusters(resolution = 0.7)
ElbowPlot(nsc_ipc.integrated, ndims = 30)
DimPlot(nsc_ipc.integrated, reduction = "umap", label = TRUE)

# Removal of GFP-negative astrocytes that come from GEMgroup 2 (10% tdtomato+ cells + 90% double negative cells) ———-
DefaultAssay(nsc_ipc.integrated) <- "RNA"
nsc_ipc.integrated <- subset(x = nsc_ipc.integrated, subset = S100b > 0 & eGFP == 0, slot = 'counts', invert = TRUE)
DimPlot(nsc_ipc.integrated, reduction = "umap", label = TRUE)

# Removal of small doublet clusters ----
nsc_ipc.integrated <- subset(x = nsc_ipc.integrated, idents = c("10", "11"), invert = TRUE)
DimPlot(nsc_ipc.integrated, reduction = "umap", label = TRUE)

# Save dataset of nscs and ipcs ----
gdata::keep(nsc_ipc.integrated, sure = TRUE)
save.image(file="nsc_ipc.RData")

# Add cell-cycle classification of G2/S/M (G0/G1 distinction not made). 
# Classification appears to work better on the raw counts data in RNA slot —---

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene
s.genes <- tolower(s.genes)
g2m.genes <- tolower(g2m.genes)
s.genes <- R.utils::capitalize(s.genes)
g2m.genes <- R.utils::capitalize(g2m.genes)
DefaultAssay(nsc_ipc.integrated) <- "RNA"
nsc_ipc.integrated <- CellCycleScoring(nsc_ipc.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
nsc_ipc.integrated$CC.Difference <- nsc_ipc.integrated$S.Score - nsc_ipc.integrated$G2M.Score
DimPlot(nsc_ipc.integrated, reduction = "umap", label = FALSE, group.by = "Phase")
gdata::keep(nsc_ipc.integrated, sure = TRUE)

# Remove all cells that are not in obvious NSC clusters ———
DimPlot(nsc_ipc.integrated, reduction = "umap", label = TRUE)
toremove <- subset(nsc_ipc.integrated, idents = c("2", "4", "6", "7"))

# Define which of these you'd definitely like to keep as stem cells 
#that have been forced together with IPCs due to effect of cell-cycle
Idents(toremove) <- toremove$Chemistry
tokeep1 <- subset(toremove, idents = "V2")
tokeep1 <- subset(tokeep1, subset = Neurod2 < 1)
tokeep1 <- subset(tokeep1, subset = Dcx < 1)
tokeep1 <- subset(tokeep1, subset = Hopx|Apoe > 5)
tokeep2 <- subset(toremove, idents = "V3")
tokeep2 <- subset(tokeep2, subset = Neurod2 < 1)
tokeep2 <- subset(tokeep2, subset = Dcx < 1)
tokeep2 <- subset(tokeep2, subset = Hopx|Apoe > 11)

# Visualise subsetting —---
DimPlot(tokeep2, reduction = "umap", label = TRUE)
DimPlot(tokeep1, reduction = "umap", label = TRUE, group.by = "Chemistry")
DimPlot(toremove, reduction = "umap", label = TRUE)

# Isolate NSC clusters, and merge with NSCs that fell in with cycling cells —---
DimPlot(nsc_ipc.integrated, reduction = "umap", label = TRUE)
nsc_ipc.integrated <- subset(nsc_ipc.integrated, idents = c('0', "1", "3", "5", "8", "9"))
nsc.integrated <- merge(nsc_ipc.integrated, tokeep1)
nsc.integrated <- merge(nsc.integrated, tokeep2)

# Split object by Chemistry (main technical effect: splitting by Batch/GEMgroup dampens biological differences) —---
nsc.list <- SplitObject(nsc.integrated , split.by = "Chemistry")
for (i in 1:length(nsc.list)) {
  nsc.list[[i]] <- SCTransform(nsc.list[[i]], return.only.var.genes = FALSE,  verbose = TRUE)
}

# Select features for identification of anchors and downstream integration ———-
nsc.features <- SelectIntegrationFeatures(object.list = nsc.list, nfeatures = 5000)
nsc.list <- PrepSCTIntegration(object.list = nsc.list, anchor.features = nsc.features, 
                                     verbose = TRUE)

# Identify anchors —---
nsc.anchors <- FindIntegrationAnchors(object.list = nsc.list, normalization.method = "SCT", 
                                            anchor.features = nsc.features, verbose = TRUE)
nsc.integrated <- IntegrateData(anchorset = nsc.anchors, normalization.method = "SCT", 
                                      verbose = TRUE)

# Visualisation/clustering ——--
nsc.integrated <- RunPCA(nsc.integrated, ndims.print = c(1:10), verbose = TRUE) %>%
        RunUMAP(dims = 1:25) %>%
        FindNeighbors(reduction = "pca", dims = 1:25) %>%
        FindClusters(resolution = 0.7)
ElbowPlot(nsc.integrated, ndims = 30)
DimPlot(nsc.integrated, reduction = "umap", label = TRUE)

# Confirm these are doublet clusters and remove ——--
nine <- FindMarkers(nsc.integrated, ident.1 = "9", ident.2 = "0", only.pos = TRUE, test.use = "t", min.diff.pct = 0.2)
ten <- FindMarkers(nsc.integrated, ident.1 = "10", ident.2 = "0", only.pos = TRUE, test.use = "t", min.diff.pct = 0.2)
eleven <- FindMarkers(nsc.integrated, ident.1 = "11", ident.2 = "0", only.pos = TRUE, test.use = "t", min.diff.pct = 0.2)
nsc.integrated <- subset(nsc.integrated, idents = c("9", "10", "11"), invert = TRUE)
DimPlot(nsc.integrated, reduction = "umap", label = TRUE)

# Save ——--
gdata::keep(nsc.integrated, sure = TRUE)
save.image(file="nsc.RData")

# End Script 4 —---

