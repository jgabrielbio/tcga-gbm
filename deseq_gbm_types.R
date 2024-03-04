# Differential Expression Analysis by studies design on GBM
# Based on DESeq2 Vignette - https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# João Gabriel - 29/08/2023


###### Installing packages------------------------

#use this function to check if each package is on the local machine
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed from Bioconductor and loaded

packages <- c("TCGAWorkflowData", "DESeq2", "TCGAbiolinks", "DT", "Glimma", "SummarizedExperiment", "tidyverse")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
rm(packages)
###### Downloading GBM Primary Tumor via biolinks--------------------

# Creating searching parameters
query_primary <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  #platform = "Illumina HiSeq", 
  sample.type = c("Primary Tumor"),
  workflow.type = c("STAR - Counts")
  #legacy = TRUE
)

# Variable only for results
query_primary$results[[1]] <-  query_primary$results[[1]]
query_primary_results <- as.data.frame(query_primary[[1]][[1]])

# Removing  duplicates 
query_primary_results <- query_primary_results[!duplicated(query_primary_results$cases),]
query_primary[[1]][[1]] <- query_primary_results

# Downloading effectively
GDCdownload(query_primary)
gbm_primary_exp <- GDCprepare(
  query = query_primary, 
  save = FALSE, 
  summarizedExperiment = TRUE
)

# Excluding NAs
gbm_primary_exp <- gbm_primary_exp[,!is.na(gbm_primary_exp$paper_IDH.status)]
gbm_primary_exp <- gbm_primary_exp[,!is.na(gbm_primary_exp$vital_status)]

# Transforming vital status names in integers
gbm_primary_exp@colData@listData <- as.data.frame(gbm_primary_exp@colData@listData) %>%
  mutate(vital_status = as.factor(vital_status)) %>%
  mutate(vital_status = if_else(vital_status %in% 'Dead', 1, 0))

# Saving object
save(gbm_primary_exp, file  = "GBM_PrimaryTumor_Illumina_HiSeq.rda")
# Load object from TCGAWorkflowData package
load("~/gbm/GBM_PrimaryTumor_Illumina_HiSeq.rda")
## Merged matrix commented ----------------
#Merging normal with normal tumor in a S4 object
#exp_gbm_Prim_Normal <- SummarizedExperiment::cbind(gbm_primary_exp , gbm_normal_exp)
#exp_gbm_Prim_Normal@colData@listData <- as.data.frame(exp_gbm_Prim_Normal@colData@listData) %>%
#mutate(vital_status = as.factor(vital_status))

# Getting expression matrix
data <- assay(gbm_primary_exp)
datatable(
  data = data[1:10,], 
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = TRUE
)

## DE analysis 1 - IDH-status_WT_vs_Mutant ----------------
ddsSE <- DESeqDataSet(gbm_primary_exp, design = ~ paper_IDH.status)
keep <- rowSums(counts(ddsSE)) >= 10
ddsSE <- ddsSE[keep,]
ddsSE <- DESeq(ddsSE)
save(ddsSE, file  = "DEA_IDH_Status.rda")
load("~/gbm/DEA_IDH_Status.rda")
resultsNames(ddsSE)
res <- results(ddsSE, name = "paper_IDH.status_WT_vs_Mutant")
# Chosen one for plots (row below)
res0.01 <- results(ddsSE, name = "paper_IDH.status_WT_vs_Mutant", alpha = 0.01)
res0.01_new <- results(ddsSE, name = "paper_IDH.status_WT_vs_Mutant", lfcThreshold = 1.5, alpha = 0.01)
dea <- as.data.frame(res)
summary(dea)
dea0.01 <- as.data.frame(res0.01)
summary(dea0.01)
dea0.01_new <- as.data.frame(res0.01_new)
save(dea0.01_new, file  = "dea_dataframe_idh_status_alpha_0.01_lfc_1.5.rda")
summary(as.data.frame(res0.01_new))
summary(res0.01_new)
# MA plot
plotMA(res)
plotMA(res0.01)
plotMA(res0.01_new)

## DE analysis 2 - vital_status -----------------
ddsSE2 <- DESeqDataSet(gbm_primary_exp, design = ~vital_status)
keep2 <- rowSums(counts(ddsSE2)) >= 10
ddsSE2 <- ddsSE2[keep2,]
ddsSE2 <- DESeq(ddsSE2)
save(ddsSE2, file  = "DEA_Vital_Status.rda")
load("~/gbm/DEA_Vital_Status.rda")
resultsNames(ddsSE2)
res2 <- results(ddsSE2, name = "vital_status")
dea2 <- as.data.frame(res2)
summary(dea2)
summary(res2)
res2_0.01 <- results(ddsSE2, name = "vital_status", alpha = 0.01)
dea2_0.01 <- as.data.frame(res2_0.01)
summary(dea2_0.01)
summary(res2_0.01)
# Chosen one for plots (row below)
res2_0.01_new <- results(ddsSE2, name = "vital_status", lfcThreshold = 1.5, alpha = 0.01)
dea2_0.01_new <- as.data.frame(res2_0.01_new)
summary(dea2_0.01_new)
summary(res2_0.01_new)
# DESeq MA plot
plotMA(res2)
plotMA(res2_0.01)
plotMA(res2_0.01_new)
## MDS plot -----------------
# For DEA by IDH status
htmlwidgets::saveWidget(glimmaMDS(ddsSE, width = 1200, height = 1200, alpha = 0.01), "mds-plot-idh-status.html")
# For DEA by Vital status
htmlwidgets::saveWidget(glimmaMDS(ddsSE2, width = 1200, height = 1200, lfcThreshold = 1.5, alpha = 0.01), "mds-plot-vital-status.html")
## Glimma MA plot -----------
htmlwidgets::saveWidget(glimmaMA(ddsSE, counts = 12083, width = 1200, height = 1200, lfcThreshold = 4, alpha = 0.001), "ma-plot-idh-status-adjusted-0.001.html")
htmlwidgets::saveWidget(glimmaMA(ddsSE2, counts = 13014, width = 1200, height = 1200, lfcThreshold = 1.5, alpha = 0.01), "ma-plot-vital-status-adjusted.html")
## Glimma MD Plot
# First by IDH Status
dea0.01$log10MeanNormCount <- log10(dea0.01$baseMean + 1)
# Plots will change if run row above
idx <-(rowSums(counts(ddsSE)) > 0) 
dea0.01 <- dea0.01[idx,] 
dea0.01$padj[is.na(dea0.01$padj)] <- 1 
status <- as.numeric(dea0.01$padj < 0.05) 
htmlwidgets::saveWidget(glMDPlot(dea0.01[idx,], 
         xval= "baseMean", 
         yval= "log2FoldChange", 
         counts = counts(ddsSE)[idx,], 
         display.columns = c("GeneID"), 
         anno = data.frame(GeneID=rownames(ddsSE)[idx]), 
         groups = gbm_primary_exp$paper_IDH.status, 
         side.xlab = "Group", 
         side.ylab = "Expression (log2)",          
         samples = gbm_primary_exp$sample, 
         status = status, 
         folder = "MDPlot_GBM.IDH_MT.vs.WT", 
         html = "index"), "glMDPlot_IDH_Status.html")
rm(dea0.01, idx, status)
# By Vital Status
dea2$log10MeanNormCount <- log10(dea2$baseMean + 1) 
# Plots will change if run row above
idx <-(rowSums(counts(ddsSE2)) > 0) 
dea2 <- dea2[idx,] 
dea2$padj[is.na(dea2$padj)] <- 1 
status <- as.numeric(dea2$padj < 0.05) 
htmlwidgets::saveWidget(glMDPlot(dea2[idx,], 
                                 xval= "baseMean", 
                                 yval= "log2FoldChange", 
                                 counts = counts(ddsSE2)[idx,], 
                                 display.columns = c("GeneID"), 
                                 anno = data.frame(GeneID=rownames(ddsSE2)[idx]), 
                                 groups = gbm_primary_exp$paper_IDH.status, 
                                 side.xlab = "Group", 
                                 side.ylab = "Expression (log2)",          
                                 samples = gbm_primary_exp$sample, 
                                 status = status, 
                                 folder = "MDPlot_GBM.Vital_Status_Alive.vs.Dead", 
                                 html = "index"), "glMDPlot_Vital_Status_pvalue_0.05.html")
# New try - P. S: Does not make difference to reduce pvalue for 0.01 
dea2_0.01$log10MeanNormCount <- log10(dea2_0.01$baseMean + 1) 
# Plots will change if run row above
idx <-(rowSums(counts(ddsSE2)) > 0) 
dea2_0.01 <- dea2_0.01[idx,] 
dea2_0.01$padj[is.na(dea2_0.01$padj)] <- 1 
status <- as.numeric(dea2_0.01$padj < 0.05) 
htmlwidgets::saveWidget(glMDPlot(dea2_0.01[idx,], 
                                 xval= "baseMean", 
                                 yval= "log2FoldChange", 
                                 counts = counts(ddsSE2)[idx,], 
                                 display.columns = c("GeneID"), 
                                 anno = data.frame(GeneID=rownames(ddsSE2)[idx]), 
                                 groups = gbm_primary_exp$paper_IDH.status, 
                                 side.xlab = "Group", 
                                 side.ylab = "Expression (log2)",          
                                 samples = gbm_primary_exp$sample, 
                                 status = status, 
                                 folder = "MDPlot_GBM.Vital_Status_Alive.vs.Dead", 
                                 html = "index"), "glMDPlot_Vital_Status_pvalue_0.01.html")
