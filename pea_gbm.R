# Pathway Expression Analysis on TCGA-GBM cohort
# Based on Cluster Profile vignette - https://yulab-smu.top/biomedical-knowledge-mining-book/
# João Gabriel - 01/03/2024
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
load("~/gbm/dea_dataframe_idh_status_alpha_0.01_lfc_1.5.rda")
gbm_dea_genes_names <- rownames(dea0.01_new)
ids <- bitr(gbm_dea_genes_names, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#Symbols genes ID
head(gbm_dea_genes_names)

x <- enrichDO(gene          = gbm_dea_genes_names,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = ids,
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
head(x)
