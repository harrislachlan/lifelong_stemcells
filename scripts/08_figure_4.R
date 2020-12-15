# Main text Figure 4
#  Load packages and define data ----
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
dir.create("plots")
ggsave(filename = "plots/all_cluters_labelled.png", width = 9, height = 5, scale = 0.618)
dev.off()
rm(list = ls())

#Pseudotime_ordering ----
load("nsc.Rdata")
nsc.integrated <- NormalizeData(nsc.integrated)
seurat <- nsc.integrated
seurat.sce <- as.SingleCellExperiment(seurat, assay = "RNA")
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
cell_embeddings <- as.data.frame(seurat@reductions[["umap"]]@cell.embeddings)
pseudotime <- seurat.sce.sling$slingPseudotime_1
cell_embeddings <- cbind(cell_embeddings, pseudotime)
c = ggplot(cell_embeddings, aes(UMAP_1, UMAP_2))
c1 = c + geom_point(aes(colour = pseudotime))
c2 = c1 + scale_colour_gradient(low="blue4", high="red")
c3 <- c2 + theme_minimal()
c4 <- c3 + NoLegend()
c4
ggsave(filename = "plots/pseudotime_shading_umap.png", width = 10, scale = 0.618)
rm(list = ls())

# Apoe & Ccnd2 expression vs pseudotime in GGPLOT ----
load("nsc.Rdata")
rnagenes <- t(GetAssayData(nsc.integrated, assay = "RNA", slot = "counts"))
rnagenes <- as.data.frame(rnagenes[ , c("Apoe","Ccnd2")])
pseudotime_nsc <- read_csv("spreadsheets/pseudotime_NSC.csv")
pseudotime_nsc$Apoe <- rnagenes$Apoe
pseudotime_nsc$Ccnd2 <- rnagenes$Ccnd2
apoe <- ggplot(pseudotime_nsc, aes(x = curve1, y = Apoe)) %>%
  + geom_point(aes(color = curve1)) + 
  scale_color_gradient(low = "blue4", high = "red") 
apoe + geom_smooth(method = "loess", color = "black") + theme_minimal()
dir.create("plots")
ggsave("plots/apoe.png", width = 5, height = 2.2)
ccnd2 <- ggplot(pseudotime_nsc, aes(x = curve1, y = Ccnd2)) %>%
  + geom_point(aes(color = curve1)) + 
  scale_color_gradient(low = "blue4", high = "red") 
ccnd2 + geom_smooth(method = "loess", color = "black") + theme_minimal()
ggsave("plots/ccnd2.png", width = 5, height = 2.2)
rm(list = ls())

#NSCs_grouped_by_states ----
load("nsc.Rdata")
Idents(nsc.integrated) <- nsc.integrated$quiescence
A <- DimPlot(nsc.integrated, reduction = "umap", group.by = "quiescence", order = c("resting"), cols = c("red", "blue4", "cyan1"), pt.size = 2)
A <- A + theme_minimal()
A <- A + NoLegend() + labs(title = "")
A
ggsave(filename = "plots/nscs_dormant_resting_proliferating.png", width = 10, scale = 0.618)
dev.off()
rm(list = ls())

#NSCs_grouped_by_states_violinplots ----
pseudotime_NSC <- read.csv("spreadsheets/pseudotime_NSC.csv")
pseudotime_NSC$quiescence <- factor(pseudotime_NSC$quiescence, levels = c("dormant", "resting", "active"))
h = ggplot(pseudotime_NSC, aes(quiescence, curve1))
h1 = h + geom_violin()
h2 = h1 + geom_violin(aes(fill = quiescence))
h3 = h2 + geom_jitter(aes(alpha = 0.001, shape = "0.01"))
h4 = h3 + labs(y = "Cell Order", x = "Quiescent Subtype", title = "Pseudotime position") + theme_minimal()
h5 <- h4 + guides(fill = "none", shape = "none", alpha = "none")
h5 <- h5 + scale_fill_manual(values=c("blue4", "cyan", "red"))
h5
ggsave(filename = "plots/NSCs_grouped_states.png", width = 6, scale = 1)

#NSCs_grouped_by_states - statistics ----
pseudotime_NSC <- read.csv("spreadsheets/pseudotime_NSC.csv")
dormant <- pseudotime_NSC %>% 
  filter(quiescence == "dormant") %>%
  pull(curve1)
resting <- pseudotime_NSC %>% 
  filter(quiescence == "resting") %>%
  pull(curve1)
active <- pseudotime_NSC %>% 
  filter(quiescence == "active") %>%
  pull(curve1)
#data is non-parametric perform wilcox.test (mann-whitney because they are not paired)
wilcox.test(dormant, resting, paired = FALSE, conf.int = TRUE)
wilcox.test(dormant, active, paired = FALSE, conf.int = TRUE)
wilcox.test(resting, active, paired = FALSE, conf.int = TRUE)
rm(list = ls())

#Dormant NSCs grouped by age ----
load("nsc.RData")
Idents(nsc.integrated) <- nsc.integrated$quiescence
nsc.integrated <- subset(nsc.integrated, idents = "dormant")
Idents(nsc.integrated) <- nsc.integrated$Age
DimPlot(nsc.integrated, group.by = "Age", pt.size = 1.5, order = c("6mo", "2mo", "1mo"), cols = c("yellow2", "orange2", "purple4")) %>%
  + theme_minimal() %>% 
  + NoLegend() + labs(title = "")
ggsave(filename = "plots/dormantNSCs_grouped_age.png", width = 5, height = 3.2)

#Dormant NSCs grouped by age - Violin plots ----
pseudotime_NSC <- read.csv("spreadsheets/pseudotime_NSC.csv")
dormantonly <- subset.data.frame(pseudotime_NSC, pseudotime_NSC$quiescence == "dormant")
g = ggplot(dormantonly, aes(Age, curve1))
g1 = g + geom_violin()
g2 = g + geom_violin(aes(fill = Age))
g3 = g2 + geom_jitter(aes(alpha = 0.01, shape = "0.01"))
g4 = g3 + labs(y = "Cell Order", title = "Pseudotime") + theme_minimal()
g4 = g4 + labs(legend = NULL) + ylim(0, 12)
g5 <- g4 + guides(fill = "none", shape = "none", alpha = "none")
g5 <- g5 + scale_fill_manual(values=c("yellow2", "orange2", "purple4"))
g5
ggsave(filename = "plots/dormantNSCs_grouped_age_violin.png", width = 5, height = 5)

#Dormant NSCs grouped by age - statistics ----
pseudotime_NSC <- read.csv("spreadsheets/pseudotime_NSC.csv")
dormantonly <- pseudotime_NSC %>%
  filter(quiescence == "dormant")
one_month_dormant <- subset(dormantonly, dormantonly$Age == "1mo")
one_month_dormant <- one_month_dormant$curve1
two_month_dormant <- subset(dormantonly, dormantonly$Age == "2mo")
two_month_dormant <- two_month_dormant$curve1
six_month_dormant <- subset(dormantonly, dormantonly$Age == "6mo")
six_month_dormant <- six_month_dormant$curve
wilcox.test(one_month_dormant, six_month_dormant, paired = FALSE, conf.int = TRUE)
wilcox.test(two_month_dormant, six_month_dormant, paired = FALSE, conf.int = TRUE)
wilcox.test(one_month_dormant, two_month_dormant, paired = FALSE, conf.int = TRUE)
kruskal.test(list(one_month_dormant, two_month_dormant, six_month_dormant))
remove(list = ls())

# End Script 8 —---



























