# Methylation associated with Gene Expression Analysis
# Based on TCGA Workflow Vignette - https://www.bioconductor.org/packages/release/workflows/vignettes/TCGAWorkflow/inst/doc/TCGAWorkflow.html#Epigenetic_analysis
# João Gabriel - 26/09/2023
#----------------------------
# Loading package
#----------------------------
packages <- c("TCGAWorkflowData", "TCGAbiolinks", "DT", "SummarizedExperiment", "tidyverse", "ComplexHeatmap", "circlize", "matlab")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
rm(packages)
#----------------------------
# Obtaining DNA methylation
#----------------------------
# Samples
load("~/gbm/clinical_final_dataset.RData")
gbm.clinic.table <- gbm.clinic.table %>% drop_na(vital_status)
samples_gbm <- slice_sample(gbm.clinic.table, weight_by = vital_status, n = 400)
samples_id <- samples_gbm$patient_id
#-----------------------------------
# 1 - Methylation
# ----------------------------------
query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  data.type = "Methylation Beta Value",
  workflow.type = "SeSAMe Methylation Beta Estimation",
  barcode = samples_id
  )
GDCdownload(query)
met <- GDCprepare(
  query = query, 
  save = FALSE
)
# remove probes with NA 
met <- met[rowSums(is.na(assay(met))) == 0,]
met <- met[,!is.na(met$paper_IDH.status)]
# Saving data
save(met, file  = "GBM_Methylation_Illumina_450.rda")
# Loading data, if necessary
load("~/gbm/GBM_Methylation_Illumina_450.rda")
#------- Searching for differentially methylated CpG sites     ----------
dmc <- TCGAanalyze_DMC(
  data = met,
  groupCol = "paper_IDH.status", # a column in the colData matrix
  group1 = "Mutant", # a type of the disease type column
  group2 = "WT", # a type of the disease column
  p.cut = 0.01,
  diffmean.cut = 0.5,
  save = FALSE,
  legend = "State",
  plot.filename = "GBM_methylation_volcano.png",
  cores = 1 # if set to 1 there will be a progress bar
)
save(dmc, file  = "dmc.rda")
load("~/gbm/dmc.rda")
#--------------------------
# DNA Methylation heatmap
#-------------------------
clinical <- plyr::rbind.fill(
  gbm.clinic.table
)

# get the probes that are Hypermethylated or Hypomethylated
# met is the same object of the section 'DNA methylation analysis'
status.col <- "status"
probes <- rownames(dmc)[grep("hypo|hyper", dmc$status, ignore.case = TRUE)]
sig.met <- met[probes,]
# top annotation, which samples are GBM-IDH mutated and not mutated.
# We will add clinical data as annotation of the samples
# we will sort the clinical data to have the same order of the DNA methylation matrix
clinical.ordered <- clinical[match(substr(colnames(sig.met), 1, 12), clinical$patient_id),]
ta <- HeatmapAnnotation(
  df = clinical.ordered[, c("gender", "vital_status", "race")],
  col = list(
    gender = c("male" = "blue", "female" = "pink"),
    vital_status = c("Dead" = "black", "Alive" = "green"),
    race = c("white" = "white", "not_white" = "brown")
  )
)
heatmap  <- Heatmap(
  matrix = assay(sig.met),
  name = "DNA methylation",
  #col = col_fun
  col = matlab::jet.colors(200),
  show_row_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  bottom_annotation = ta,
  column_title = "DNA Methylation",
  use_raster = TRUE
)
# Save to pdf
png("heatmap_matlab_updated.png",width = 600, height = 400)
draw(heatmap, annotation_legend_side =  "bottom")
dev.off()

qsquare_vital_hyper <- chisq.test(table(clinical.ordered$vital_status, dmc$status))
