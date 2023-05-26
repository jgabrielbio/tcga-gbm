### TCGA - GBM Survival Analysis ####
#Based on it Datacamp (https://www.datacamp.com/tutorial/survival-analysis-R)
# João Gabriel - 26/05/23

### Loading packages ----
library(survival)
library(survminer)
library(dplyr)
if(!require('survminer')) {
     install.packages('survminer')
     library('survminer')
}

###Kaplan-Meier----
#Loading final dataset from clinical_data.R on console, first.
gbm_clinic_surv <- gbm.clinic.table
gbm_clinic_surv <- gbm_clinic_surv %>% mutate(age_group = ifelse(age >=50, "old", "young")) %>% relocate(age_group, .after = age) %>% mutate(days_to_last_follow_up = gbm.clinic$days_to_last_follow_up )
gbm_clinic_surv <- within(gbm_clinic_surv, {
     race <- factor(race, labels = c('white', 'not_white'))
     gender <- factor(gender, labels = c("female", "male"))
     age_group <- factor(age_group, labels = c("old", "young"))
     pharmaceutical_treatment <- factor(pharmaceutical_treatment, labels = c("yes", "no"))
     prior_treatment <- factor(prior_treatment, labels = c("Yes", "No"))
     radiation_treatment <- factor(radiation_treatment, labels = c("yes", "no"))
})
gbm_clinic_surv <- gbm_clinic_surv  %>%
      mutate_if(is.character, as.factor)
gbm_clinic_surv <- gbm_clinic_surv %>% 
     mutate(vital_status = if_else(vital_status %in% 'Dead', 1, 0))
surv_object <- Surv(time = gbm_clinic_surv$days_to_last_follow_up, event = gbm_clinic_surv$vital_status)
surv_object  

fit1 <- survfit(surv_object ~ gender, data = gbm_clinic_surv)
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


### Hazard Ratio  ----

fit_coxph <- coxph(surv_object ~ race + gender + age_group + pharmaceutical_treatment + prior_treatment + radiation_treatment, data = gbm_clinic_surv)

# Increased Risk (HR>1), Decreased HR<1
ggforest(fit_coxph, data = gbm_clinic_surv)
