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

#script to determine threshold for calling cells in G0/G1 ----
#load NSCs and IPCs
load("nsc_ipc.RData")

#nsc_ipc.integrated <- NormalizeData(nsc_ipc.integrated)
DefaultAssay(nsc_ipc.integrated) <- "RNA"

#add cc phase and plot
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene
s.genes <- tolower(s.genes)
g2m.genes <- tolower(g2m.genes)
s.genes <- R.utils::capitalize(s.genes)
g2m.genes <- R.utils::capitalize(g2m.genes)
nsc_ipc.integrated <- CellCycleScoring(nsc_ipc.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
nsc_ipc.integrated$CC.Difference <- nsc_ipc.integrated$S.Score - nsc_ipc.integrated$G2M.Score
DimPlot(nsc_ipc.integrated, features = "Phase")
ggsave(filename = "plots/phase.png", width = 6, height = 4, scale = 1)

#extract bona-fide cells in S-phase 
Idents(nsc_ipc.integrated) <- nsc_ipc.integrated$Phase
s_phase <- subset(nsc_ipc.integrated, idents = "S")
DimPlot(s_phase)


#extract bona-fide qNSCs and label these as such in the metadata
Idents(nsc_ipc.integrated) <- nsc_ipc.integrated$integrated_snn_res.0.7
DimPlot(nsc_ipc.integrated, label = TRUE)
nsc_ipc.integrated$Phase <- as.character(nsc_ipc.integrated$Phase)
nsc_ipc.integrated$Phase <- ifelse(nsc_ipc.integrated$integrated_snn_res.0.7 == "1", yes = "G0", no = nsc_ipc.integrated$Phase)
bona_qnsc <- subset(nsc_ipc.integrated, idents = "1")
DimPlot(bona_qnsc)

#merge bona-fide qNSCs and s-phase NSPCs
s_qnsc <- merge(x = s_phase, y = bona_qnsc)

#Add G1 index genes
DefaultAssay(s_qnsc) <- "RNA"
s_qnsc[["percentMCM2"]]<- PercentageFeatureSet(s_qnsc, pattern = "Mcm2")
s_qnsc[["percentMCM3"]] <- PercentageFeatureSet(s_qnsc, pattern = "Mcm3")
s_qnsc[["percentMCM4"]] <- PercentageFeatureSet(s_qnsc, pattern = "Mcm4")
s_qnsc[["percentMCM5"]] <- PercentageFeatureSet(s_qnsc, pattern = "Mcm5")
s_qnsc[["percentMCM6"]] <- PercentageFeatureSet(s_qnsc, pattern = "Mcm6")
s_qnsc[["percentMCM7"]] <- PercentageFeatureSet(s_qnsc, pattern = "Mcm7")
s_qnsc[["percentCCNE1"]] <- PercentageFeatureSet(s_qnsc, pattern = "Ccne1")
s_qnsc[["percentCCNE2"]] <- PercentageFeatureSet(s_qnsc, pattern = "Ccne2")
s_qnsc[["percentPCNA"]] <- PercentageFeatureSet(s_qnsc, pattern = "Pcna")
#Binarise
s_qnsc@meta.data$'percentMCM2' <- ifelse(test = s_qnsc@meta.data$'percentMCM2' > 0, yes = 1, no = 0) 
s_qnsc@meta.data$'percentMCM3' <- ifelse(test = s_qnsc@meta.data$'percentMCM3' > 0, yes = 1, no = 0) 
s_qnsc@meta.data$'percentMCM4' <- ifelse(test = s_qnsc@meta.data$'percentMCM4' > 0, yes = 1, no = 0) 
s_qnsc@meta.data$'percentMCM5' <- ifelse(test = s_qnsc@meta.data$'percentMCM5' > 0, yes = 1, no = 0) 
s_qnsc@meta.data$'percentMCM6' <- ifelse(test = s_qnsc@meta.data$'percentMCM6' > 0, yes = 1, no = 0) 
s_qnsc@meta.data$'percentMCM7' <- ifelse(test = s_qnsc@meta.data$'percentMCM7' > 0, yes = 1, no = 0) 
s_qnsc@meta.data$'percentCCNE1' <- ifelse(test = s_qnsc@meta.data$'percentCCNE1' > 0, yes = 1, no = 0) 
s_qnsc@meta.data$'percentCCNE2' <- ifelse(test = s_qnsc@meta.data$'percentCCNE2' > 0, yes = 1, no = 0)
s_qnsc@meta.data$'percentPCNA' <- ifelse(test = s_qnsc@meta.data$'percentPCNA' > 0, yes = 1, no = 0)

#Add G1-index genes together
s_qnsc[["G1_index"]] <- (s_qnsc@meta.data$'percentMCM2' 
                                 + s_qnsc@meta.data$'percentMCM3'
                                 + s_qnsc@meta.data$'percentMCM4' 
                                 + s_qnsc@meta.data$'percentMCM5'
                                 + s_qnsc@meta.data$'percentMCM6'
                                 + s_qnsc@meta.data$'percentMCM7'
                                 +s_qnsc@meta.data$'percentCCNE1'
                                 +s_qnsc@meta.data$'percentCCNE2'
                                +s_qnsc@meta.data$'percentPCNA' )
                  
Idents(s_qnsc) <- s_qnsc$Chemistry
v3_s_qnsc <- subset(s_qnsc, idents = "V3")
v2_s_qnsc <- subset(s_qnsc, idents = "V2")
VlnPlot(s_qnsc, features = "G1_index")
#here the G1 cells correspond to the bona-fide G0, this just hasn't been updated in the metadata
V1 <- VlnPlot(v3_s_qnsc, features = "G1_index", group.by = "Phase")
V2 <- VlnPlot(v2_s_qnsc, features = "G1_index", group.by = "Phase")
plot_grid(V1, V2, nrow = 2)
G <- arrangeGrob(V1, V2, nrow = 2)
ggsave(filename = "plots/hourglass.png", G, width = 4, height = 6, scale = 1)
rm(list =ls())

#Independent pan-set of genes ----
load("nsc.RData")
Idents(nsc.integrated) <- nsc.integrated$quiescence
DefaultAssay(object = nsc.integrated) <- "RNA"
nsc.integrated <- NormalizeData(nsc.integrated, assay = "RNA")
nsc.integrated$quiescence <- as.factor(nsc.integrated$quiescence)
nsc.integrated$quiescence <- factor(nsc.integrated$quiescence, levels = c("dormant", "resting", "active"))
V7 <- VlnPlot(nsc.integrated, features = "Gins1", group.by = "quiescence", cols = c("blue4", "cyan1", "red")) + NoLegend() + labs(x = NULL, y = NULL)
V8 <- VlnPlot(nsc.integrated, features = "Gins2", group.by = "quiescence", cols = c("blue4", "cyan1", "red")) + NoLegend() + labs(x = NULL, y = NULL)
V9 <- VlnPlot(nsc.integrated, features = "Gins4", group.by = "quiescence", cols = c("blue4", "cyan1", "red")) + NoLegend() + labs(x = NULL, y = NULL)
V10 <- VlnPlot(nsc.integrated, features = "E2f1", group.by = "quiescence", cols = c("blue4", "cyan1", "red")) + NoLegend() + labs(x = NULL, y = NULL)
V11 <- VlnPlot(nsc.integrated, features = "E2f3", group.by = "quiescence", cols = c("blue4", "cyan1", "red")) + NoLegend() + labs(x = NULL, y = NULL)
V12 <- VlnPlot(nsc.integrated, features = "Cdt1", group.by = "quiescence", cols = c("blue4", "cyan1", "red")) + NoLegend() + labs(x = NULL, y = NULL)
plot_grid(V7, V8, V9, V10, V11, V12, nrow = 2)
G <- arrangeGrob(V7, V8, V9, V10, V11, V12, nrow = 2)
ggsave(filename = "plots/independent_pan_G1_S.png", G, width = 12, height = 8, scale = 1)

V13 <- VlnPlot(nsc.integrated, features = "Mcm2", group.by = "quiescence", cols = c("blue4", "cyan1", "red")) + NoLegend() + labs(x = NULL, y = NULL)
V14 <- VlnPlot(nsc.integrated, features = "Mcm6", group.by = "quiescence", cols = c("blue4", "cyan1", "red")) + NoLegend() + labs(x = NULL, y = NULL)
plot_grid(V13, V14, nrow = 1)
G <- arrangeGrob(V13, V14, nrow = 1)
ggsave(filename = "plots/MCM2_MCM6.png", G, width = 8, height = 5, scale = 1)
dev.off()
rm(list = ls())





#No of mRNA and gens per cell----
load("nsc.RData")
Idents(nsc.integrated) <- nsc.integrated$quiescence
DefaultAssay(object = nsc.integrated) <- "RNA"
nsc.integrated <- NormalizeData(nsc.integrated, assay = "RNA")
nsc.integrated$quiescence <- as.factor(nsc.integrated$quiescence)
nsc.integrated$quiescence <- factor(nsc.integrated$quiescence, levels = c("dormant", "resting", "active"))
V1 <- VlnPlot(nsc.integrated, features = "nFeature_RNA", group.by = "quiescence", cols = c("blue4", "cyan1", "red")) + NoLegend() + labs(x = NULL, y = NULL)
V2 <- VlnPlot(nsc.integrated, features = "nCount_RNA", group.by = "quiescence", cols = c("blue4", "cyan1", "red")) + NoLegend() + labs(x = NULL, y = NULL)
plot_grid(V1, V2, nrow = 1)
G <- arrangeGrob(V1, V2, nrow = 1)
ggsave(filename = "plots/RNA_GENES.png", G, width = 8, height = 4, scale = 1)
dev.off()
rm(list = ls())

#No of mRNA and gens per cell - STATISTICS ---- 
load("nsc.RData")
nUMI <- nsc.integrated@meta.data$nCount_RNA
nUMI <- as.tibble(cbind(nsc.integrated@meta.data$nCount_RNA, nsc.integrated@meta.data$quiescence))
dormantUMI <- subset(nUMI, nUMI$V2 == "dormant")
dormantUMI <- as.numeric(dormantUMI$V1)
restingUMI <- subset(nUMI, nUMI$V2 == "resting")
restingUMI <- as.numeric(restingUMI$V1)
activeUMI <- subset(nUMI, nUMI$V2 == "active")
activeUMI <- as.numeric(activeUMI$V1)
wilcox.test(dormantUMI, restingUMI, paired = FALSE, conf.int = TRUE)
wilcox.test(dormantUMI, activeUMI, paired = FALSE, conf.int = TRUE)
wilcox.test(restingUMI, activeUMI, paired = FALSE, conf.int = TRUE)

nGene <- nsc.integrated@meta.data$nFeature_RNA
nGene <- as.tibble(cbind(nsc.integrated@meta.data$nFeature_RNA, nsc.integrated@meta.data$quiescence))
dormantnGene <- subset(nGene, nGene$V2 == "dormant")
dormantnGene <- as.numeric(dormantnGene$V1)
restingnGene <- subset(nGene, nGene$V2 == "resting")
restingnGene <- as.numeric(restingnGene$V1)
activenGene <- subset(nGene, nGene$V2 == "active")
activenGene <- as.numeric(activenGene$V1)
wilcox.test(dormantnGene, restingnGene, paired = FALSE, conf.int = TRUE)
wilcox.test(dormantnGene, activenGene, paired = FALSE, conf.int = TRUE) 
wilcox.test(restingnGene, activenGene, paired = FALSE, conf.int = TRUE)
rm(list = ls())


