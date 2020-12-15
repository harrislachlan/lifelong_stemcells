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
library(gam)

# Slingshot of all NSCs ----
load("nsc.RData")
nsc.integrated <- NormalizeData(nsc.integrated)
seurat <- nsc.integrated
seurat.sce <- as.SingleCellExperiment(seurat, assay = "RNA")
DimPlot(nsc.integrated)
seurat.sce.sling <- slingshot(seurat.sce, reducedDim = 'UMAP')
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(seurat.sce.sling$slingPseudotime_1, breaks=100)]
plot(reducedDims(seurat.sce.sling)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(seurat.sce.sling), lwd=2, col='black')

# look at the 2000 most variable genes
Y <- log1p(SummarizedExperiment::assays(seurat.sce.sling)$logcounts)
var1K <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:2000]
Y <- Y[var1K,]


# fit a GAM with a loess term for pseudotime
t <- seurat.sce.sling$slingPseudotime_1
gam.pval_guillemot <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})
gam.pval_guillemot <- as.data.frame(gam.pval_guillemot) %>%
  rownames_to_column(var = "genes") %>%
  filter(gam.pval_guillemot < 0.01)

##### Shin et al

input_file <- "GSE71485_Single_TPM.txt"
tpm <- read.delim(input_file)
seurat_object <- CreateSeuratObject(counts = tpm, 
                                    project = "Waterfall", min.cells = 3, min.features = 1)


#Add MT content
seurat_object[["percentMito"]] <- PercentageFeatureSet(seurat_object, pattern = "^Mt")
mitoHi <- 10
nGeneLo <- 500
seurat_object <- subset(x = seurat_object, subset = nFeature_RNA > nGeneLo & percentMito < mitoHi)

seurat_object.features <- FindVariableFeatures(seurat_object, nfeatures = 3000)
seurat_object.features <- ScaleData(seurat_object.features)
seurat_object <- RunPCA(seurat_object.features, verbose = TRUE)
sigPC <- 5
seurat_object <- RunUMAP(seurat_object, dims = 1:sigPC)
seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:sigPC)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
DimPlot(seurat_object, dims = c(1, 2), reduction = "umap", label = TRUE)
seurat_object <- NormalizeData(seurat_object, assay = "RNA")
#FeaturePlot(seurat_object, "Hopx", reduction = "umap")
#FeaturePlot(seurat_object, "Eomes", reduction = "umap")
#FeaturePlot(seurat_object, "Ccnd2", reduction = "umap")
#FeaturePlot(seurat_object, "CFP",reduction = "umap")

seurat_object <- subset(seurat_object, idents = c("0", "2"))
seurat_object <- NormalizeData(seurat_object, assay = "RNA")


#Slingshot - identifying pseudotime genes
seurat <- seurat_object
seurat.sce <- as.SingleCellExperiment(seurat, assay = "RNA")
seurat.sce.sling <- slingshot(seurat.sce, clusterLabels = seurat@active.ident, start.clus = "0", reducedDim = 'UMAP')
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(seurat.sce.sling$slingPseudotime_1, breaks=100)]
plot(reducedDims(seurat.sce.sling)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(seurat.sce.sling), lwd=2, col='black')
t <- seurat.sce.sling$slingPseudotime_1
SummarizedExperiment::assayNames(seurat.sce.sling)


# look at the 2000 most variable genes
Y <- log1p(SummarizedExperiment::assays(seurat.sce.sling)$logcounts)
var1K <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:2000]
Y <- Y[var1K,]

# fit a GAM with a loess term for pseudotime
t <- seurat.sce.sling$slingPseudotime_1
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})
gam.pval_shin <- as.data.frame(gam.pval) %>%
  rownames_to_column(var = "genes") %>%
  filter(gam.pval < 0.01)

intersect(gam.pval_guillemot$genes, gam.pval_shin$genes)

length(intersect(rownames(seurat_object@assays$RNA@counts), rownames(nsc.integrated@assays$RNA@counts)))

#Of the ~14000 gene in common in dataset, 1768 were associated with pseudotime in our
#data (12.6%). By chance would expect 12.6% of 516 genes in Shin dataset to appear in ours,
#whereas the overlap was much larger it was 70.5% (364/(364+152)).
