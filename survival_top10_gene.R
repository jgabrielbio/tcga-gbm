# Installing and Loading Libraries            
packages_bioconductor = c("TCGAbiolinks", "SummarizedExperiment", "genefilter", "maftools")

#use this function to check if each package is on the local machine
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed from Bioconductor and loaded
package.check <- lapply(packages_bioconductor, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
packages_cran = c("DT", "tidyverse","dplyr", "tibble", "stringr", "data.table", "survminer", "survival")

#use this function to check if each package is on the local machine
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed from CRAN and loaded
package.check <- lapply(packages_cran, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

rm(packages_cran, packages_bioconductor, package.check)
#download archive below on github folder
load("~/Pesquisa/GBM/Scripts/Dados Clínicos/clinical_final_dataset.RData")
## Reading gbm Maf files ---------------------------
#download archive below on github folder
load("~/Pesquisa/GBM/Scripts/Dados Clínicos/genomico/Genomica_Datas.RData")
##### Survival analysis based on grouping of PI3K mutation status ------
gbm.clinic$Tumor_Sample_Barcode <- gbm.clinic$bcr_patient_barcode
gbm.clinic$Overall_Survival_Status <- ifelse(gbm.clinic$vital_status %in% c("Dead"),1,0)
gbm.pi3k.surv = mafSurvival(maf = gbm.mafclin, clinicalData = gbm.clinic, genes = 'TP53', time = 'days_to_last_follow_up', Status = 'Overall_Survival_Status', showConfInt = TRUE)
