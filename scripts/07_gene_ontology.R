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

# Gene ontology of all differentially expressed genes ----
# 6 month dormant NSCs vs 1 month dormant NSCs 
setwd("spreadsheets/")
ageGO <- read_csv(file = "sig_genes_age.csv") %>%
  filter(p_val_adj < 0.05) %>%
  pull(X1) %>%
  enrichGO(OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', 
           ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05)
write_csv(ageGO@result, "ageGO.csv")
setwd("..")

# RestingvsDormant
setwd("spreadsheets/")
rest_dorm_GO <- read_csv(file = "sig_genes_resting_dormant.csv") %>%
  filter(p_val_adj < 0.05) %>%
  pull(X1) %>%
  enrichGO(OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', 
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05)
write_csv(rest_dorm_GO@result, "rest_dorm_GO.csv")
setwd("..")

# Restingvsactive
setwd("spreadsheets/")
rest_active_GO <- read_csv(file = "sig_genes_resting_active.csv") %>%
  filter(p_val_adj < 0.05) %>%
  pull(X1) %>%
  enrichGO(OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', 
           ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05)
write_csv(rest_active_GO@result, "rest_active_GO.csv")
setwd("..")

# Gene ontology of all upregulated genes ----

# 1vs6 upregulated
setwd("spreadsheets/")
ageGO <- read_csv(file = "sig_genes_age.csv") %>%
  filter(p_val_adj < 0.05 & avg_diff > 0) %>%
  pull(X1) %>%
  enrichGO(OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', 
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05)
write_csv(ageGO@result, "age_upreg_GO.csv")
setwd("..")

# RestingvsDormant upregulated
setwd("spreadsheets/")
rest_dorm_GO <- read_csv(file = "sig_genes_resting_dormant.csv") %>%
  filter(p_val_adj < 0.05 & avg_diff > 0) %>%
  pull(X1) %>%
  enrichGO(OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', 
           ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05)
write_csv(rest_dorm_GO@result, "rest_dorm_upreg_GO.csv")
setwd("..")

# RestingvsActive upregulated
setwd("spreadsheets/")
rest_active_GO <- read_csv(file = "sig_genes_resting_active.csv")%>%
  filter(p_val_adj < 0.05 & avg_diff > 0) %>%
  pull(X1) %>%
  enrichGO(OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', 
           ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05)
write_csv(rest_active_GO@result, "rest_active_upreg_GO.csv")
setwd("..")

# Gene ontology of all downregulated genes ----

# 1vs6 downregulated
setwd("spreadsheets/")
ageGO <- read_csv(file = "sig_genes_age.csv") %>%
  filter(p_val_adj < 0.05 & avg_diff < 0) %>%
  pull(X1) %>%
  enrichGO(OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', 
           ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05)
write_csv(ageGO@result, "age_downreg_GO.csv")
setwd("..")

# RestingvsDormant downregulated
setwd("spreadsheets/")
rest_dorm_GO <- read_csv(file = "sig_genes_resting_dormant.csv") %>%
  filter(p_val_adj < 0.05 & avg_diff < 0) %>%
  pull(X1) %>%
  enrichGO(OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', 
           ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05)
write_csv(rest_dorm_GO@result, "rest_dorm_downreg_GO.csv")
setwd("..")

# Restingvsactive downregulated
setwd("spreadsheets/")
rest_active_GO <- read_csv(file = "sig_genes_resting_active.csv") %>%
  filter(p_val_adj < 0.05 & avg_diff < 0) %>%
  pull(X1) %>%
  enrichGO(OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', 
           ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05)
write_csv(rest_active_GO@result, "rest_active_downreg_GO.csv")
setwd("..")
remove(list = ls())

# End Script 7 —---