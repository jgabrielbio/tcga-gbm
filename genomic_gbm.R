######## Maftools ##########
# Joao Gabriel - 17/03/2023 

# maftools: Summarize, Analyze and Visualize Mutation Anotated Files (MAF) Files
# URL: https://www.bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html
# citation("maftools")
# version 3.16




# path to TCGA gbm MAF file

## Installing packages
packages_bioconductor <- c("TCGAbiolinks","maftools")
packages_cran <- c("tidyverse")
#use this function to check if each package is on the local machine
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed from Bioconductor and loaded
package.check <- lapply(packages_bioconductor, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
package.check <- lapply(packages_cran, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

rm(packages_cran, packages_bioconductor, package.check)
setwd('C:/Users/joaog/Documents/Pesquisa/GBM/Scripts')
load("~/Pesquisa/GBM/Scripts/Genomica_Datas.RData")
## Reading gbm Maf files ---------------------------

# download MAF aligned against hg38
# it saves data inside a GDCdata and project name directory 
query <- GDCquery(
  project = "TCGA-GBM", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE, 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query)

#Remover duplicatas dentro do query 
query_results <- as.data.frame(query[[1]][[1]])
query_results <- query_results[!duplicated(query_results$cases),]
query[[1]][[1]] <- query_results

# GDC
maf <- GDCprepare(query)
sort(colnames(maf))

# MAF object contains main maf file, summarized data and any associated sample annotations
gbm.maf <- read.maf(maf = maf, useAll = T) 

# checking
getSampleSummary(gbm.maf) #  samples
getGeneSummary(gbm.maf) #  genes (hugo)
getClinicalData(gbm.maf) #  samples, no clinical data  
getFields(gbm.maf) #  variables 

# writes an output file
write.mafSummary(maf = gbm.maf, basename = 'gbm.maf')

## Reading clinical indexed data ------------------

# Filtered clinical data
# clinical <- gbm.clinic
# same as in data portal
clinical <- GDCquery_clinic(project = "TCGA-GBM", type = "clinical", save.csv = FALSE)
sort(colnames(clinical))
clinical$Tumor_Sample_Barcode <- clinical$bcr_patient_barcode 
clinical$time <- clinical$days_to_death
clinical$time[is.na(clinical$days_to_death)] <- clinical$days_to_last_follow_up[is.na(clinical$days_to_death)]
clinical <- clinical %>%
  mutate(time = as.numeric(time)) %>% 
  mutate(vital_status = if_else(vital_status %in% 'Dead', 1, 0))

# create object for survival analysis 
gbm.mafclin <- read.maf(maf = maf, clinicalData = clinical, isTCGA = T)


## Visualizing -----------------------------------

# displays variants in each sample and variant types summarized by Variant_Classification
plotmafSummary(maf = gbm.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
# Most mutated genes - amplified
mafbarplot(maf = gbm.maf)
# oncoplot for top ten mutated genes (costumize oncoplots!)
oncoplot(maf = gbm.maf, top = 10)

# Transitions and transvertions 
gbm.titv = titv(maf = gbm.maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = gbm.titv)
# Mudança de aminoácidos
#lollipop plot for TP53
lollipopPlot(
  maf = gbm.maf,
  gene = 'TP53',
  AACol = 'Amino_acids',
  showMutationRate = TRUE,
  labelPos = 882
)
plotProtein(gene = "TP53", refSeqID = "NM_000546")
#Rainfall plot
rainfallPlot(maf = gbm.maf, detectChangePoints = TRUE, pointSize = 0.4)
#Variant Allele Frequencies
plotVaf(maf = gbm.maf, vafCol = 'i_TumorVAF_WU')
## Analyzing -----------------------------------

#exclusive/co-occurance event analysis on top 10 mutated genes. 
somaticInteractions(maf = gbm.maf, top = 20, pvalue = c(0.05, 0.01))

#Detecting cancer driver genes based on positional clustering
# Oncodrive has been superseeded by OncodriveCLUSTL. 
# See http://bg.upf.edu/group/projects/oncodrive-clust.php
gbm.maf.sig <- oncodrive(maf = gbm.maf, AACol = 'Amino_acids', minMut = 5, pvalMethod = 'zscore')
head(gbm.maf.sig)
plotOncodrive(res = gbm.maf.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
# Domains of aminoacids changes
gbm.pfam = pfamDomains(maf = gbm.maf, AACol = 'Amino_acids', top = 10)
gbm.pfam$proteinSummary[,1:7, with = FALSE]
gbm.pfam$domainSummary[,1:3, with = FALSE]
#Survival analysis based on grouping of TP53 mutation status
mafSurvival(maf = gbm.mafclin, genes = 'TP53', time = 'time', Status = "vital_status", isTCGA = TRUE)

#Using top 20 mutated genes to identify a set of genes (of size 2) to predict poor prognostic groups
prog_geneset = survGroup(maf = gbm.mafclin, top = 10, geneSetSize = 2, time = "days_to_last_follow_up", Status = "vital_status", verbose = FALSE)
print(prog_geneset)
#Plot KM
mafSurvGroup(maf = gbm.mafclin, geneSet = c("MUC16", "PIK3R1"), time = "days_to_last_follow_up", Status = "vital_status")

##Comparing two cohorts (MAFs)

#Primary APL MAF
# primary.apl = system.file("extdata", "APL_primary.maf.gz", package = "maftools")
# primary.apl = read.maf(maf = primary.apl)

# #Relapse APL MAF
# relapse.apl = system.file("extdata", "APL_relapse.maf.gz", package = "maftools")
# relapse.apl = read.maf(maf = relapse.apl)

# #Considering only genes which are mutated in at-least in 5 samples in one of the cohort to avoid bias due to genes mutated in single sample.
# pt.vs.rt <- mafCompare(m1 = primary.apl, m2 = relapse.apl, m1Name = 'Primary', m2Name = 'Relapse', minMut = 5)
# print(pt.vs.rt)
# forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1)

#Clinical enrichment analysis
fab.ce.f1 = clinicalEnrichment(maf = gbm.mafclin, clinicalFeature = 'tissue_or_organ_of_origin')
fab.ce.f1$groupwise_comparision[p_value < 0.01]
plotEnrichmentResults(enrich_res = fab.ce.f1, pVal = 0.01, geneFontSize = 0.5, annoFontSize = 0.6)

#Drug-Gene Interactions
dgi = drugInteractions(maf = gbm.maf, fontSize = 0.75)
dnmt3a.dgi = drugInteractions(genes = "TP53", drugs = TRUE)
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]

#Oncogenic Signaling Pathways
OncogenicPathways(maf = gbm.maf)
PlotOncogenicPathways(maf = gbm.maf, pathways = "PI3K")
