#### TCGA-GMB Clinical Data Analysis #######
# Beatriz Stransky - 22/04/2023


# install.packages(c('TCGAbiolinks', 'tidyverse', 'skimr', 'finalfit', 'tableone', dependencies = TRUE))
library(TCGAbiolinks)
library(tidyverse)
library(skimr)
library(finalfit)
library(tableone)
library(rio)


## 1. processamento de dados para salvar arquivo .csv ------

# download clinical data using the GDCquery_xxx() from TCGABiolinks
gbm.clinical.raw <- GDCquery_clinic(project = "TCGA-GBM", type = "Clinical", save.csv = FALSE)

# save it in other object
gbm.clinic <- gbm.clinical.raw

# CHECKING DATASET
class(gbm.clinic)
dim(gbm.clinic)
glimpse(gbm.clinic)
View(gbm.clinic)

# CLEANING DATA
skim(gbm.clinic)

# Removendo observaçoes duplicadas
gbm.clinic <- gbm.clinic %>%
  distinct_at("submitter_id", .keep_all = TRUE)

# Removing Logical variables
logical.vars <- gbm.clinic %>%
  select_if(is.logical) %>% 
  names()
gbm.clinic <- gbm.clinic %>%
  select(-all_of(logical.vars))

# Checando inicio e fim 
head(gbm.clinic)
tail(gbm.clinic, 20)

gbm.clinic <- gbm.clinic %>%
 slice(-(600:617))

# Contando NAs 
gbm.clinic.nas <- gbm.clinic %>%
  summarise_all(~ sum(is.na(.))) %>% # summarise_each
  pivot_longer(cols = everything(), names_to = "var", values_to = "NAs") %>% 
  arrange(desc(NAs))

# gbm.clinic.nas2 <- gbm.clinic %>%
#   purrr::map(~ sum(is.na(.))) 
# gbm.clinic.nas2 <- as.data.frame(gbm.clinic.nas2) 

# Unique/distinct values
gbm.clinic.unique <- gbm.clinic %>%
  summarise_if(is.character, n_distinct, na.rm = FALSE) %>% 
  pivot_longer(cols = everything(), names_to = "var", values_to = "unique") %>% 
  arrange(desc(unique))

table(gbm.clinic$created_datetime)
gbm.clinic <- gbm.clinic %>% 
  select(-created_datetime)

# teste ID
sum(gbm.clinic$submitter_id == gbm.clinic$bcr_patient_barcode)
teste_id <- gbm.clinic %>%
  filter(!submitter_id %in% bcr_patient_barcode)
rm(teste_id)

# Checando variaveis numericas
gbm.clinic %>%
  select(is.numeric) %>% 
  summary()

unique(gbm.clinic$days_to_diagnosis)
sum(gbm.clinic$days_to_diagnosis, na.rm = T)

gbm.clinic <- gbm.clinic %>% 
  select(-days_to_diagnosis)

# checar valores
sum(gbm.clinic$age_at_diagnosis %/% 365 == gbm.clinic$age_at_index, na.rm = T)
sum(-1 * (gbm.clinic$days_to_birth %/% 365) == gbm.clinic$age_at_index, na.rm = T)

# remover variaveis de id
gbm.clinic_id <- gbm.clinic %>%
  select((ends_with("_id")))
gbm.clinic <- gbm.clinic %>%
  select(!(ends_with("_id")))
rm(gbm.clinic_id)

skim(gbm.clinic)

# taming data
gbm.clinic <- gbm.clinic %>%
  mutate_if(is.character, as.factor) %>%
  mutate(bcr_patient_barcode = as.character(bcr_patient_barcode)) %>% 
  select(project, bcr_patient_barcode, everything())

write_csv(gbm.clinic, "gbm.clinic.csv")


## 2. modificando gbm.clinic para imprimir tabela ------

# remover variaveis com unico valor
vars_unique <- gbm.clinic.unique %>%
  filter(unique == 1)
gbm.clinic.table <- gbm.clinic %>%
  select(-(vars_unique$var))

gbm.clinic.table <- gbm.clinic.table %>%
  select(c(bcr_patient_barcode, age_at_index, race, gender, ethnicity, vital_status, prior_treatment, treatments_pharmaceutical_treatment_or_therapy, treatments_radiation_treatment_or_therapy))


# Changing variables names
names(gbm.clinic.table)
names(gbm.clinic.table) <- names(gbm.clinic.table) %>%
  str_remove("treatments_")

gbm.clinic.table <- gbm.clinic.table %>%
  rename(
    age = "age_at_index",
    pharmaceutical_treatment = "pharmaceutical_treatment_or_therapy",
    radiation_treatment = "radiation_treatment_or_therapy",
    patient_id = "bcr_patient_barcode"
  )
names(gbm.clinic.table)

skim(gbm.clinic.table)

# taming levels*
gbm.clinic.fct <- gbm.clinic.table %>%
  select_if(is.factor) %>% 
  summary()
gbm.clinic.fct

# agregating levels*
gbm.clinic.table <- gbm.clinic.table %>%
  mutate(
    race = fct_collapse(race, not_white = c("asian", "black or african american")),
    prior_treatment = fct_recode(prior_treatment, NULL = "Not Reported"),
    race = fct_recode(race, NULL = "not reported"),
    gender = fct_recode(gender, NULL = "not reported"),
    ethnicity = fct_recode(ethnicity, NULL = "not reported"),
    vital_status = fct_recode(vital_status, NULL = "Not Reported"),
    pharmaceutical_treatment = fct_recode(pharmaceutical_treatment, NULL = "not reported"),
    radiation_treatment = fct_recode(radiation_treatment, NULL = "not reported")
  )


## 3. Tabelas para dados clinicos ------

# TABLE ONE
gbm.clinic.tb1 <- gbm.clinic.table %>%
  select(-patient_id)

CreateTableOne(data = gbm.clinic.tb1)
myVars <- names(gbm.clinic.tb1)
catVars <- names(gbm.clinic.tb1 %>%
                   select_if(is.factor))

tab1 <- CreateTableOne(vars = myVars, data = gbm.clinic.tb1, factorVars = catVars)
print(tab1, showAllLevels = TRUE, formatOptions = list(big.mark = ","))


# FINAL FIT
explanatory <- gbm.clinic.table %>%
    select(-c(vital_status, patient_id)) %>%
    names
#dependent <- 'vital_status'
 
clinic_descript <- gbm.clinic.table %>%
  #summary_factorlist(dependent, explanatory,
  summary_factorlist(explanatory,
    na_to_p = FALSE,
    na_include = FALSE, add_dependent_label = TRUE
  )

# visualization
clinic_descript




