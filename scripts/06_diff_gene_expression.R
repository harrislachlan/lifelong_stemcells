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

# Purpose: SCTransform of 1mo and 6mo NSCs to maximise no. of Pearson Residuals returned and DGE ----
# Load data and subset cell of interest, return all Pearson Residuals by re-applying SCTransform
#first analysis: DGE of dormant and resting NSCs in 6 vs 1 month old mice. First remove 2-mo old cells and active cells. 
load("nsc.RData")
Idents(nsc.integrated) <- nsc.integrated$Age
nsc.integrated <- subset(nsc.integrated, idents = "2mo", invert = TRUE)
Idents(nsc.integrated) <- nsc.integrated$quiescence
nsc.integrated <- subset(nsc.integrated, idents = "active", invert = TRUE)
nsc.integrated <- SCTransform(nsc.integrated, return.only.var.genes = FALSE)

# Determine sig genes and FDR adjust p-value
Idents(nsc.integrated) <- nsc.integrated$Age
sig_genes_age <- FindMarkers(nsc.integrated, assay = "SCT", slot = "scale.data", ident.1 = "6mo",
                             ident.2 = "1mo", min.pct = 0.0, 
                             logfc.threshold = 0, test.use = "t")
sig_genes_pvalue <- sig_genes_age$p_val
FDR <- p.adjust(sig_genes_pvalue, "fdr")
sig_genes_age$p_val_adj <- FDR
sig_genes_age <- sig_genes_age %>% filter(pct.1 >= 0.2|pct.2>= 0.2)

# save spreadsgeet
setwd("spreadsheets/")
write.csv(sig_genes_age, "sig_genes_age.csv")
DT::datatable(sig_genes_age, colnames = c('gene' = 2))
setwd("..")
remove(list = ls())


# Purpose :SCTransform of NSCs to maximise no. of Pearson Residuals returned and DGE of resting vs dormant ----
# Load data
load("nsc.RData")
nsc.integrated <- SCTransform(nsc.integrated, return.only.var.genes = FALSE)

# Determine sig genes
Idents(nsc.integrated) <- nsc.integrated$quiescence
sig_genes_resting_dormant <- FindMarkers(nsc.integrated, assay = "SCT", ident.1 = "resting",
                         ident.2 = "dormant", slot = "scale.data", min.pct = 0.0,
                         logfc.threshold = 0, test.use = "t")
sig_genes_pvalue <- sig_genes_resting_dormant$p_val
FDR <- p.adjust(sig_genes_pvalue, "fdr")
sig_genes_resting_dormant$p_val_adj <- FDR
sig_genes_resting_dormant <- sig_genes_resting_dormant %>% filter(pct.1 >= 0.2|pct.2 >= 0.2)

setwd("spreadsheets/")
write.csv(sig_genes_resting_dormant, "sig_genes_resting_dormant.csv")
setwd("..")
remove(list = ls())

# Purpose :SCTransform of NSCs to maximise no. of Pearson Residuals returned and DGE of resting vs active ----

# Load data
load("nsc.RData")
nsc.integrated <- SCTransform(nsc.integrated, return.only.var.genes = FALSE)

# Determine sig genes
Idents(nsc.integrated) <- nsc.integrated$quiescence
sig_genes_resting_active <- FindMarkers(nsc.integrated, assay = "SCT", ident.1 = "resting",
                         ident.2 = "active", slot = "scale.data", min.pct = 0.0, 
                         logfc.threshold = 0, test.use = "t")
sig_genes_pvalue <- sig_genes_resting_active$p_val
FDR <- p.adjust(sig_genes_pvalue, "fdr")
sig_genes_resting_active$p_val_adj <- FDR
sig_genes_resting_active <- sig_genes_resting_active %>% filter(pct.1 >= 0.2| pct.2 >= 0.2)

setwd("spreadsheets/")
write.csv(sig_genes_resting_active , "sig_genes_resting_active.csv")
setwd("..")
remove(list = ls())
