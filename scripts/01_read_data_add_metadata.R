# Purpose: define data, create seurat object, attach metadata associated with experiment/sample ----

# Required packages, save session info ----
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
today=lubridate::today()
sink(glue("sessionInfo_{today}.txt"))
sessionInfo()
sink()
rm(today)

# Read in and create Seurat object,  if downloaded from GEO, delete file prefixes "GSE159768_". ----
input_folder <- "data/"
seurat_object <- Read10X(data.dir = input_folder)
seurat_object <- CreateSeuratObject(counts = seurat_object, project = "Aggregate")

# Add GEMgroup ----
seurat_object@meta.data$GEMgroup[grepl("-1", rownames(seurat_object@meta.data))] <- "1"
seurat_object@meta.data$GEMgroup[grepl("-2", rownames(seurat_object@meta.data))] <- "2"
seurat_object@meta.data$GEMgroup[grepl("-3", rownames(seurat_object@meta.data))] <- "3"
seurat_object@meta.data$GEMgroup[grepl("-4", rownames(seurat_object@meta.data))] <- "4"
seurat_object@meta.data$GEMgroup[grepl("-5", rownames(seurat_object@meta.data))] <- "5"
seurat_object@meta.data$GEMgroup[grepl("-6", rownames(seurat_object@meta.data))] <- "6"
seurat_object@meta.data$GEMgroup[grepl("-7", rownames(seurat_object@meta.data))] <- "7"
seurat_object@meta.data$GEMgroup[grepl("-8", rownames(seurat_object@meta.data))] <- "8"

# Add Batch no - got this from 10X website ----
seurat_object@meta.data$Batch[grepl("-1", rownames(seurat_object@meta.data))] <- "1"
seurat_object@meta.data$Batch[grepl("-2|-3", rownames(seurat_object@meta.data))] <- "2"
seurat_object@meta.data$Batch[grepl("-4", rownames(seurat_object@meta.data))] <- "3"
seurat_object@meta.data$Batch[grepl("-5", rownames(seurat_object@meta.data))] <- "4"
seurat_object@meta.data$Batch[grepl("-6", rownames(seurat_object@meta.data))] <- "5"
seurat_object@meta.data$Batch[grepl("-7", rownames(seurat_object@meta.data))] <- "6"
seurat_object@meta.data$Batch[grepl("-8", rownames(seurat_object@meta.data))] <- "7"

# Add Age of samples ----
seurat_object@meta.data$Age[grepl("-1", rownames(seurat_object@meta.data))] <- "1mo"
seurat_object@meta.data$Age[grepl("-2|-3|-4|-5", rownames(seurat_object@meta.data))] <- "2mo"
seurat_object@meta.data$Age[grepl("-6|-7|-8", rownames(seurat_object@meta.data))] <- "6mo"

# GFP and tdtomato or combined ----
seurat_object@meta.data$cells[grepl("-1|-4|-5|-6|-7|-8", rownames(seurat_object@meta.data))] <- "combined"
seurat_object@meta.data$cells[grepl("-2", rownames(seurat_object@meta.data))] <- "tdTomato"
seurat_object@meta.data$cells[grepl("-3", rownames(seurat_object@meta.data))] <- "Gfp"

# Add 10x chemistry ----
seurat_object@meta.data$Chemistry[grepl("-1|-6|-7|-8", rownames(seurat_object@meta.data))] <- "V3"
seurat_object@meta.data$Chemistry[grepl("-2|-3|-4|-5", rownames(seurat_object@meta.data))] <- "V2"

# Add MT content ----
seurat_object[["percentMito"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")

# Add sex ----
SexMatrix1 <- t(GetAssayData(seurat_object, assay = "RNA", slot = "counts"))
SexMatrix1 <- as.data.frame(SexMatrix1[ , c("Ddx3y","Uty","Eif2s3y","Xist","Tsix") ]) 
seurat_object[["sexScore"]] <- (SexMatrix1$Ddx3y + SexMatrix1$Uty + SexMatrix1$Eif2s3y) - (SexMatrix1$Xist + SexMatrix1$Tsix)  
seurat_object[["Sex"]] <- ifelse(seurat_object[["sexScore"]] < 0, "Female", "Male") 

# End Script 1 ----

