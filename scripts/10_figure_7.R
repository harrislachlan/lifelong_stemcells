# Main text Figure 7
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

# Ascl1/Huwe1/ID4 vs Age (1v6) in violin plots ----
load("nsc.RData")
Idents(nsc.integrated) <- nsc.integrated$Age
nsc.integrated <- subset(nsc.integrated, idents = "2mo", invert = TRUE)
Idents(nsc.integrated) <- nsc.integrated$quiescence
nsc.integrated <- subset(nsc.integrated, idents = "active", invert = TRUE)
Idents(nsc.integrated) <- nsc.integrated$Age
nsc.integrated <- NormalizeData(nsc.integrated, assay = "RNA")
DefaultAssay(object = nsc.integrated) <- "RNA"
V1 <- VlnPlot(nsc.integrated, features = "Ascl1", group.by = "Age",pt.size = 0.2)  + NoLegend() + labs(x = NULL, y = NULL)
V2 <- VlnPlot(nsc.integrated, features = "Id4", group.by = "Age", pt.size = 0.2) + NoLegend() + labs(x = NULL, y = NULL)
V3 <- VlnPlot(nsc.integrated, features = "Huwe1", group.by = "Age", pt.size = 0.2) +ylim(0, 2.5) + NoLegend() + labs(x = NULL, y = NULL)
G <- arrangeGrob(V1, V2, V3, nrow = 1)
plot_grid(V1, V2, V3, ncol = 3)
ggsave(filename = "plots/ascl1_huwe1_age.png", G, width = 5, height = 3, scale = 1)
dev.off()
rm(list = ls())

