# Purpose: add cell-cycle and tdTomato metadata ----
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
# Load data ----
load("nsc.RData")

# Identify cells that express detectable levels of genes stably, highly expressed across G1 & S-phase ----

DefaultAssay(nsc.integrated) <- "RNA"
nsc.integrated <- NormalizeData(nsc.integrated, assay = "RNA")
nsc.integrated[["percentMCM2"]]<- PercentageFeatureSet(nsc.integrated, pattern = "Mcm2")
nsc.integrated[["percentMCM3"]] <- PercentageFeatureSet(nsc.integrated, pattern = "Mcm3")
nsc.integrated[["percentMCM4"]] <- PercentageFeatureSet(nsc.integrated, pattern = "Mcm4")
nsc.integrated[["percentMCM5"]] <- PercentageFeatureSet(nsc.integrated, pattern = "Mcm5")
nsc.integrated[["percentMCM6"]] <- PercentageFeatureSet(nsc.integrated, pattern = "Mcm6")
nsc.integrated[["percentMCM7"]] <- PercentageFeatureSet(nsc.integrated, pattern = "Mcm7")
nsc.integrated[["percentCCNE1"]] <- PercentageFeatureSet(nsc.integrated, pattern = "Ccne1")
nsc.integrated[["percentCCNE2"]] <- PercentageFeatureSet(nsc.integrated, pattern = "Ccne2")
nsc.integrated[["percentPCNA"]] <- PercentageFeatureSet(nsc.integrated, pattern = "Pcna")

# Binarise as positive or negative
nsc.integrated@meta.data$'percentMCM2' <- ifelse(test = nsc.integrated@meta.data$'percentMCM2' > 0, yes = 1, no = 0) 
nsc.integrated@meta.data$'percentMCM3' <- ifelse(test = nsc.integrated@meta.data$'percentMCM3' > 0, yes = 1, no = 0) 
nsc.integrated@meta.data$'percentMCM4' <- ifelse(test = nsc.integrated@meta.data$'percentMCM4' > 0, yes = 1, no = 0) 
nsc.integrated@meta.data$'percentMCM5' <- ifelse(test = nsc.integrated@meta.data$'percentMCM5' > 0, yes = 1, no = 0) 
nsc.integrated@meta.data$'percentMCM6' <- ifelse(test = nsc.integrated@meta.data$'percentMCM6' > 0, yes = 1, no = 0) 
nsc.integrated@meta.data$'percentMCM7' <- ifelse(test = nsc.integrated@meta.data$'percentMCM7' > 0, yes = 1, no = 0) 
nsc.integrated@meta.data$'percentCCNE1' <- ifelse(test = nsc.integrated@meta.data$'percentCCNE1' > 0, yes = 1, no = 0) 
nsc.integrated@meta.data$'percentCCNE2' <- ifelse(test = nsc.integrated@meta.data$'percentCCNE2' > 0, yes = 1, no = 0)
nsc.integrated@meta.data$'percentPCNA' <- ifelse(test = nsc.integrated@meta.data$'percentPCNA' > 0, yes = 1, no = 0)

# Add this index together; max value 9, min value is 0 ----
nsc.integrated[["G1_index"]] <- (nsc.integrated@meta.data$'percentMCM2' 
                                + nsc.integrated@meta.data$'percentMCM3'
                                + nsc.integrated@meta.data$'percentMCM4' 
                                + nsc.integrated@meta.data$'percentMCM5'
                                + nsc.integrated@meta.data$'percentMCM6'
                                + nsc.integrated@meta.data$'percentMCM7'
                                +nsc.integrated@meta.data$'percentCCNE1'
                                +nsc.integrated@meta.data$'percentCCNE2'
                                +nsc.integrated@meta.data$'percentPCNA')


# Define cell as active if expresses > 1 index gene sequenced using V2 or >3 index gene if sequenced using V3 ----
nsc.integrated[["G1"]] <- ifelse(nsc.integrated[["G1_index"]] > 1 & nsc.integrated@meta.data$'Chemistry' == "V2", yes = "1", no = "0")
nsc.integrated[["G1"]] <- ifelse(nsc.integrated[["G1_index"]] > 3 & nsc.integrated@meta.data$'Chemistry' == "V3", yes = "1", no = nsc.integrated@meta.data$'G1')
nsc.integrated@meta.data[["G1"]] <- as.numeric(nsc.integrated@meta.data[["G1"]])


# Scoring whether a cell is tdTomato+ or tdTomato negative. Add WPRE and bGHpolyA together (both elements of tdTomato recombined site) for more robustness and sensitivity of td detection
rnagenes <- t(GetAssayData(nsc.integrated, assay = "RNA", slot = "counts"))
rnagenes <- as.data.frame(rnagenes[ , c("WPRE","tdTomatoLoxP", "bGHpolyA")])
rnagenes$td <- rnagenes$WPRE + rnagenes$bGHpolyA
nsc.integrated <- AddMetaData(object = nsc.integrated, metadata = rnagenes)

# Enforcing of thresholds ----
# V3 of kit - if count of td 4 or more and tdtomatoloxp less than 8 than cell is tdTomato positive. 
nsc.integrated@meta.data$'tdthresh' <- ifelse(test = nsc.integrated@meta.data$'td' > 3 & nsc.integrated@meta.data$'Chemistry' == "V3", yes = nsc.integrated@meta.data$'td', no = 0) 
nsc.integrated@meta.data$'tdthresh' <- ifelse(test = nsc.integrated@meta.data$'tdTomatoLoxP' > 7 & nsc.integrated@meta.data$'Chemistry' == "V3", yes = 0, no = nsc.integrated@meta.data$'tdthresh') 
# V2 of kit - if count of td 2 or more and tdtomatoloxp less than 4 than cell is tdTomato positive. 
nsc.integrated@meta.data$'tdthresh' <- ifelse(test = nsc.integrated@meta.data$'td' > 1 & nsc.integrated@meta.data$'Chemistry' == "V2", yes = nsc.integrated@meta.data$'td', no = nsc.integrated@meta.data$'tdthresh') 
nsc.integrated@meta.data$'tdthresh' <- ifelse(test = nsc.integrated@meta.data$'tdTomatoLoxP' > 3 & nsc.integrated@meta.data$'Chemistry' == "V2", yes = 0, no = nsc.integrated@meta.data$'tdthresh')

# Binarisation of tdTomato thresholds, td pos = 1 or td neg = 0 
nsc.integrated@meta.data$'tdTomato+' <- ifelse(test = nsc.integrated@meta.data$'tdthresh' > 0 , yes = 1, no = 0)


# Adding metadata to define dormant (tdtomato- & quiescent), resting (tdtomato+ & quiescent) and active (tdtomtato+ or - & active) ----
nsc.integrated[["quiescence"]] <- (nsc.integrated@meta.data$'G1' + nsc.integrated@meta.data$'tdTomato+')
nsc.integrated[["quiescence"]] <- ifelse(nsc.integrated@meta.data$'tdTomato+' == 0 & nsc.integrated@meta.data$'G1' < 1, 
                                          yes = 'dormant', no = nsc.integrated@meta.data$'quiescence')
nsc.integrated[["quiescence"]] <- ifelse(nsc.integrated@meta.data$'tdTomato+' > 0 & nsc.integrated@meta.data$'G1' < 1, 
                                         yes = 'resting', no = nsc.integrated@meta.data$'quiescence')
nsc.integrated[["quiescence"]] <- ifelse(nsc.integrated@meta.data$'tdTomato+' <= 1 & nsc.integrated@meta.data$'G1' > 0, 
                                         yes = 'active', no = nsc.integrated@meta.data$'quiescence')

# Next two lines ensure no cell scored a G2M or S by CellCycleScoring function in Seurat are accidentally called G0 by my thresholding ----
nsc.integrated[["quiescence"]] <- ifelse(nsc.integrated$Phase == "G2M", yes = 'active',
                                          no = nsc.integrated@meta.data$'quiescence')
nsc.integrated[["quiescence"]] <- ifelse(nsc.integrated$Phase == "S", yes = 'active',
                                          no = nsc.integrated@meta.data$'quiescence')

# Define G1/G0 cells into plot ----
nsc.integrated[["Phase"]] <- ifelse(nsc.integrated$quiescence == "dormant"|nsc.integrated$quiescence == "resting",
                                    yes = 'Go', no = nsc.integrated$Phase)
nsc.integrated$Phase <- as.character(nsc.integrated$Phase)
DimPlot(nsc.integrated, group.by = "Phase")
DimPlot(nsc.integrated, group.by = "quiescence")

# Purpose: remove 3 outlier cells that will artificially elongate pseudotime, in Slingshot analysis below. 
nsc.integrated <- 
  subset(nsc.integrated, cells = c("GGGCCATCAGTCCCGA-7", 
                                   "CTGATCCAGCGCCTAC-8", 
                                   "GTCTCGTCACCACCAG-3"),
         invert = TRUE)

# Save data ----
gdata::keep(nsc.integrated, sure = TRUE)
save.image(file="nsc.RData")

# Slingshot of all NSCs ----
nsc.integrated <- NormalizeData(nsc.integrated)
seurat <- nsc.integrated
seurat.sce <- as.SingleCellExperiment(seurat, assay = "RNA")
DimPlot(nsc.integrated)
seurat.sce.sling <- slingshot(seurat.sce, reducedDim = 'UMAP')
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(seurat.sce.sling$slingPseudotime_1, breaks=100)]
plot(reducedDims(seurat.sce.sling)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(seurat.sce.sling), lwd=2, col='black')
pseudotimevalues <- slingPseudotime(seurat.sce.sling)
Age <- seurat.sce.sling@colData@listData$Age
quiescence <- seurat.sce.sling@colData@listData$quiescence
pseudotimevalues <- cbind(pseudotimevalues, Age)
pseudotimevalues <- cbind(pseudotimevalues, quiescence)

# Save cell ordering
dir.create("spreadsheets")
setwd('spreadsheets/')
write.csv(pseudotimevalues, file = "pseudotime_NSC.csv")
setwd("..")
remove(list = ls())

# End Script 5



