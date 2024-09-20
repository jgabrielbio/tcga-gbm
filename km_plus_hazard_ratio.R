### Análise de Sobrevivência do TCGA - GBM ####

#Base: Datacamp (https://www.datacamp.com/tutorial/survival-analysis-R)

### Anotações ----
### Equivalências dataset survival com GBM
### futime = days_to_last_follow_up
### fustat = vital_status
### patient_age = age



### Carregando pacotes ----
library(survival)
library(dplyr)
# library(tidyverse)

if(!require('survminer')) {
     install.packages('survminer')
     library('survminer')
}

###Curva de Kaplan-Meier----
#Loading final dataset from clinical_data.R on console
#load("~/Pesquisa/GBM/Scripts/Dados Clínicos/clinical_final_dataset.RData")
load("~/gbm/clinical_final_dataset.RData")
gbm_clinic_surv <- gbm.clinic.table
gbm_clinic_surv <- gbm_clinic_surv %>% 
  mutate(age_group = ifelse(age >=50, "old", "young")) %>% 
  relocate(age_group, .after = age) %>% 
  mutate(days_to_last_follow_up = gbm.clinic$days_to_last_follow_up ) %>%
  rename(
    gender = "sex")
gbm_clinic_surv <- gbm_clinic_surv  %>%
      mutate_if(is.character, as.factor)

gbm_clinic_surv <- gbm_clinic_surv %>% 
     mutate(vital_status = if_else(vital_status %in% 'Dead', 1, 0))

surv_object <- Surv(time = gbm_clinic_surv$days_to_last_follow_up, event = gbm_clinic_surv$vital_status)
surv_object
fit0 <- survfit(Surv(days_to_last_follow_up, vital_status) ~ 1, data = gbm_clinic_surv)
#to know median
ggsurvplot(fit0, data = gbm_clinic_surv, surv.median.line = "hv", pval = TRUE)
print(fit0)
fit1 <- survfit(surv_object ~ sex, data = gbm_clinic_surv)
ggsurvplot(fit1, data = gbm_clinic_surv, pval = TRUE)

fit2 <- survfit(surv_object ~ race, data = gbm_clinic_surv)
ggsurvplot(fit2, data = gbm_clinic_surv, pval = TRUE)

fit3 <- survfit(surv_object ~ age_group, data = gbm_clinic_surv)
ggsurvplot(fit3, data = gbm_clinic_surv, pval = TRUE)

fit4 <- survfit(surv_object ~ pharmaceutical_treatment, data = gbm_clinic_surv)
ggsurvplot(fit4, data = gbm_clinic_surv, pval = TRUE)

fit5 <- survfit(surv_object ~ prior_treatment, data = gbm_clinic_surv)
ggsurvplot(fit5, data = gbm_clinic_surv, pval = TRUE)

fit6 <- survfit(surv_object ~ radiation_treatment, data = gbm_clinic_surv)
ggsurvplot(fit6, data = gbm_clinic_surv, pval = TRUE)

fit7 <- survfit(surv_object~1, data = gbm_clinic_surv)
ggsurvplot(fit7, data = gbm_clinic_surv, pval = TRUE)


### Modelo de risco ----
fit_coxph <- coxph(surv_object ~ race + sex + age_group + pharmaceutical_treatment + prior_treatment + radiation_treatment, data = gbm_clinic_surv)

# Augmented risk HR>1, reducted risk HR<1
ggforest(fit_coxph, data = gbm_clinic_surv)
