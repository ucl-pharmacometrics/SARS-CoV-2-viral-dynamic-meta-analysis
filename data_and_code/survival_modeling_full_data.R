# Code to underake survival analysis on SARS-CoV-2
library(tidyverse)
library(survival)
library(ggthemes)
library(ggpubr)
library(survminer) 
library(ggfortify)
library(forestmodel)
#library(MASS)
library(ggplot2)
####### read in master data ####
viral_dynamics <- read_csv(here::here("data_and_code", "SARSCoV2_viral_dynamics_200620.csv"))
# length(unique(viral_dynamics$id_anon))
# 645
# unique(viral_dynamics$drug_quality)
# This analysis will be limited to upper or lower respiratory tract or stool, and
#  drug quality zero or 1 (i.e. knowledge of whether or not patient received
#     antivirals but not necessarily time varying)
#  Rationale: multiple imputation or other methods may not be
#      reliable for drug data, particularly since such heterogeneous drugs
#   blood, urine, breast milk ocular may be negative for virus even when 
#     resp tract is positive, so if these are taken at later time point may
#     falsely indicate viral cleaance.
#  
viral_dynamics <- viral_dynamics %>%
  filter(drug_quality <= 2) %>% # removed quality 3 
  filter(sample_site_group != "blood/plasma") %>% 
  filter(sample_site_group != "breast_milk" ) %>% 
  filter(sample_site_group != "ocular") %>%  
  filter(sample_site_group != "urine")  
# unique(viral_dynamics$sample_site_group) 
# length(unique(viral_dynamics$id_anon))
# unique(viral_dynamics$drug_quality)
# remove unspecified
viral_dynamics <- viral_dynamics[viral_dynamics$drug_name != "antiviralsbutnotspecified", ]
viral_dynamics <- viral_dynamics[viral_dynamics$drug_name != "convalescentplasma+antivirals", ]
viral_dynamics <- viral_dynamics[viral_dynamics$drug_name != "convplasma", ]
# conv plasma alone is just 1 patient
##  
# length(unique(viral_dynamics$id_anon))
#  354
# unique(viral_dynamics$drug_quality)
# unique(viral_dynamics$drug_name[viral_dynamics$drug_quality == 0])
# This analysis will not look at time varying drug effects, rather seeks 
#  signals for the pharmacodynamic model-based approach.  therefore put each antiviral as
#  non-time varying
viral_dynamics <- viral_dynamics %>%
  group_by(id_anon) %>%
  mutate(drug_azit_ever = case_when(max(drug_azit) == 1 ~ "1_azi", 
                                    TRUE ~ "0_azi")) %>% 
  mutate(drug_cqhcq_ever = case_when(max(drug_cqhcq) == 1 ~ "1_cqhcq", 
                                     TRUE ~ "0_cqhcq")) %>% 
  mutate(drug_ifn_ever = case_when(max(drug_ifn) == 1 ~ "1_ifn", 
                                   TRUE ~ "0_ifn")) %>%
  mutate(drug_ifn_alpha_ever = case_when(max(drug_ifn_alpha) == 1 ~ "1_ifn_a",
                                         TRUE ~ "0_ifn_a")) %>%
  mutate(drug_ifn_beta_ever = case_when(max(drug_ifn_beta) == 1 ~ "1_ifn_b",
                                        TRUE ~ "0_ifn_b")) %>%
  mutate(drug_lpvr_ever = case_when(max(drug_lpvr) == 1 ~ "1_lpvr", 
                                    TRUE ~ "0_lpvr")) %>%  
  mutate(drug_riba_ever = case_when(max(drug_riba) == 1 ~ "1_riba", 
                                    TRUE ~ "0_riba")) %>%  
  mutate(drug_remd_ever = case_when(max(drug_remd) == 1 ~ "1_remd", 
                                    TRUE ~ "0_remd")) %>%  
  mutate(drug_thym_ever = case_when(max(drug_thym) == 1 ~ "1_thym", 
                                    TRUE ~ "0_thym")) %>% 
  mutate(drug_umif_ever = case_when(max(drug_umif) == 1 ~ "1_umif", 
                                    TRUE ~ "0_umif")) %>% 
  mutate(drug_any_ever = case_when(max(drug_azit, drug_cqhcq, drug_ifn,  
                                       drug_lpvr, drug_riba, drug_thym,
                                       drug_umif) == 1 ~ "1_any", 
                                   TRUE ~ "0_none"))
# filter out patients for whom we don't know what drugs they got
# table(viral_dynamics$drug_name)
# antiviralsbutnotspecified 
# 119 length(unique(viral_dynamics$id_anon))
# filter out patients for whom we don't know what drugs they got
# table(viral_dynamics$drug_name)
# antiviralsbutnotspecified 
# 119 length(unique(viral_dynamics$id_anon))
# label for specific combinations
viral_dynamics <- viral_dynamics %>%
  # monotherapy
  mutate(drug_lpvr_mono_ever = case_when(drug_lpvr_ever == "1_lpvr" & 
                                           drug_riba_ever == "0_riba" &
                                           drug_ifn_beta_ever == "0_ifn_b" ~ "1_lpvr",
                                         TRUE ~ "0_lpvr")) %>%  
  mutate(drug_riba_mono_ever = case_when(drug_lpvr_ever == "0_lpvr" & 
                                           drug_riba_ever == "1_riba" &
                                           drug_ifn_beta_ever == "0_ifn_b" ~ "1_riba",
                                         TRUE ~ "0_riba")) %>%  
  # lpvr + riba
  mutate(drug_lpvr_riba_ever = case_when(drug_lpvr_ever == "1_lpvr" & 
                                           drug_riba_ever == "1_riba" &
                                           drug_ifn_beta_ever == "0_ifn_b" ~ "1_lpvr_riba",
                                         TRUE ~ "0_lpvr_riba")) %>%
  # lpvr + riba + ifn b
  mutate(drug_lpvr_riba_ifn_beta_ever = case_when(drug_lpvr_ever == "1_lpvr" & 
                                                    drug_riba_ever == "1_riba" &
                                                    drug_ifn_beta_ever == "1_ifn_b" ~ "1_lpvr_ifnb_riba",
                                                  TRUE ~ "0_lpvr_ifnb_riba")) %>% 
  # umif mono
  mutate(drug_umif_mono_ever = case_when(drug_umif_ever == "1_umif" & 
                                           drug_ifn_alpha_ever == "0_ifn_a" ~ "1_umif",
                                         TRUE ~ "0_umif")) %>% 
  mutate(drug_ifn_alpha_mono_ever = case_when(drug_umif_ever == "0_umif" & 
                                                drug_ifn_alpha_ever ==  "0_ifn_a" ~ "1_ifn_a",
                                              TRUE ~ "0_ifn_a")) %>%  
  mutate(drug_umif_ifn_alpha_ever = case_when(drug_umif_ever == "1_umif" & 
                                                drug_ifn_alpha_ever ==  "1_ifn_a" ~ "1_ifna_umif",
                                              TRUE ~ "0_ifna_umif")) 
# Now select the "first-last" below limit of detection sample
#  i.e. if there are a series of negative samples prior to last take the first of these
# find, for each patient, the event variable = first of the last string of consecutive
#  1s in blod_1
## examples   0 0 0 0 0 0 1 0 0 0 0 X 1 1 1 1 1  event in x (1)
##            0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 X  censored in x (0)
##            0 0 0 0 1 0 0 0 0 0 1 1 1 1 0 0 0 0 X 1 1 event in x (1)
# Firstly remove all the patients with BLD = 0 on last observed day
viral_dynamics <- viral_dynamics %>%
  group_by(id_anon) %>%
  # if any of the final values are above lod, set censored = 1
  mutate(detectable = case_when(min(blod_1[t_symp_days == max(t_symp_days)]) == 1 ~ 0,
                                TRUE ~ 1))
vd_cens <- viral_dynamics %>%
  filter(detectable == 1) %>%
  arrange(id_anon, -t_symp_days, lod_value_log10_cpml) %>%
  # take last time point and lowest detection limit assay if there were more than one
  group_by(id_anon) %>%
  slice(1)
# length(unique(vd_cens$id_anon)) nrow(vd_cens)
vd_event <- viral_dynamics %>%
  filter(detectable == 0) %>%
  mutate(all_bld = case_when(min(blod_1) == 1 ~ 1,
                             TRUE ~ 0))
vd_event_all_bld <- vd_event %>%
  filter(all_bld == 1) %>%
  # select lowest time and lowest lod value
  arrange(id_anon, t_symp_days, lod_value_log10_cpml) %>%
  group_by(id_anon) %>%
  slice(1)
# length(unique(vd_event_all_bld$id_anon)) nrow(vd_event_all_bld)
vd_event <- vd_event %>%
  filter(all_bld == 0) 
#length(unique(vd_event$id_anon)) + length(unique(vd_cens$id_anon))  + length(unique(vd_event_all_bld$id_anon)) 
# Now for event pick first in a series of bld vlues
vd_event <- vd_event %>%
  # sort by reverse time and assay sensitivity so if 2 assays on same day we pick most sensitive
  arrange(id_anon, -t_symp_days, blod_1, ) %>%
  group_by(id_anon) %>%
  # select BLQ values where next one is not BLQ
  mutate(bld_select = case_when(blod_1 != lead(blod_1) ~ 1,
                                TRUE ~ 0)) %>%
  # Find latest time of that value
  mutate(time_event = max(t_symp_days[bld_select == 1])) %>%
  filter(t_symp_days == time_event) %>%  
  # Find lowest Blod value
  mutate(lod_event = min(lod_value_log10_cpml)) %>%
  # sum(vd_event$bld_select)
  filter(bld_select == 1)
# length(unique(vd_event$id_anon)) 
# nrow(vd_event) + nrow(vd_event_all_bld)  + nrow(vd_cens)
#now bind back together
vd_event <- vd_event[ , -which(names(vd_event) %in% c("detectable", "all_bld", 
                                                      "bld_select", "time_event",           
                                                      "lod_event"))]
vd_event_all_bld <- vd_event_all_bld[ , -which(names(vd_event_all_bld) %in% c("detectable", "all_bld"))]
vd_cens <- vd_cens[ , -which(names(vd_cens) %in% c("detectable"))]
# viral dynamics now has 1 row per patient
viral_dynamics <- rbind(vd_event, vd_event_all_bld, vd_cens)
viral_dynamics <- viral_dynamics %>%
  mutate(time_to_event = t_symp_days) %>%
  mutate(event = blod_1)
# length(unique(viral_dynamics $id_anon)) nrow(viral_dynamics) 

viral_dynamics <- viral_dynamics %>%
  mutate(sample_site_group_cat = case_when(sample_site_group == "lower_resp_tract" ~ 
                                             "2=lower_resp_tract",
                                           sample_site_group == "upper_resp_tract" ~ 
                                             "1=upper_resp_tract",
                                           sample_site_group == "stool/rectal" ~ 
                                             "3=stool/rectal"
  )) %>%
  mutate(vl_quality_cat = case_when(vl_quality == 1 ~ "1",
                                    vl_quality == 2 ~ "2",
                                    vl_quality == 3 ~ "2",)) %>%
  mutate(age_cat = case_when(age_y < 40 ~ "0-39",
                             age_y >= 40 & age_y < 60 ~ "40-59",
                             age_y >= 60 & age_y < 80 ~ "60-79",
                             age_y >= 80 ~ "80+")) %>%
  mutate(male_rounded = case_when(round(male) == 1 ~ "male",
                                  round(male) == 0 ~ "female")) %>%
  mutate(disease_status_cat = case_when(disease_status == 0 ~ "0=Asymptomatic",
                                        disease_status == 1 ~ "1=Mild",
                                        disease_status == 2 ~ "2=Moderate",
                                        disease_status == 3 ~ "3=Severe"))



######## survival modelling ##########
# 1. Check data
# length(unique(viral_dynamics$id_anon))  
# 354 individuals
# unique(viral_dynamics$drug_name)  
#  unique(viral_dynamics$drug_quality)
#  table(viral_dynamics$male)
#  sum(is.na(viral_dynamics$male))
#  table(viral_dynamics$disease_status)
#  sum(is.na(viral_dynamics$disease_status))
#  table(viral_dynamics$age_y)
#  summary(viral_dynamics$age_y)
#  unique(viral_dynamics$paper[is.na(viral_dynamics$age_y)])
# No missing data, for sex need to consider how to handle those with proportions
#  option 1: multiple imputation, but this does not use prior knowledge on
#      proportions in studies
#  option 2: bootstrap through proportion data to randomly assign
#  option 3: set to most likely i.e. prop <0.5 = female
#  conclusion: use option 2 and compare coefficients with option 3
#   depending on similarity, take forward option 2 or 3 to final model
######### 2. Model building #######
######### Model 1, base model ########
km.mod1 <- survfit(Surv(time_to_event, event) ~ 1, 
                   data = viral_dynamics)
summary(km.mod1, times = seq(0, 60, by = 10))
surv.plot.1 <- ggsurvplot(km.mod1, conf.int = TRUE, pval = FALSE, 
                          risk.table = TRUE, risk.table.y.text.col = TRUE,
                          xlab = "Time (days)", ylab = "p(negativity)")
surv.plot.1
# Model building:
#  Start with data structure covariates: site and lod value
#  Then demographics
#  Then antivirals yes/no
######## Model 2, sample site #########
# Rationale: shedding in e.g. stool and LRT seems longer than
#  URT so important to correct for this a priori
# table(viral_dynamics$sample_site_group)
km.mod2 <- survfit(Surv(time_to_event, event) ~ sample_site_group_cat, 
                   data = viral_dynamics)
table(viral_dynamics$sample_site_group_cat)
cox.mod2 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat, 
                  data = viral_dynamics)
summary(cox.mod2)
surv.plot.2 <- ggsurvplot(km.mod2, conf.int = FALSE, pval = TRUE, 
                          risk.table = TRUE, risk.table.y.text.col = TRUE,
                          xlab = "Time (days)", ylab = "p(negativity)") +
  NULL
surv.plot.2
######## Model 3, viral load quality #########
# Rationale: less sensitive assays will be more likely to call negative?
table(viral_dynamics$vl_quality_cat)
km.mod3 <- survfit(Surv(time_to_event, event) ~ vl_quality_cat, 
                   data = viral_dynamics)
cox.mod3 <- coxph(Surv(time_to_event, event) ~ vl_quality_cat, 
                  data = viral_dynamics)
summary(cox.mod3)
surv.plot.3 <- ggsurvplot(km.mod3, conf.int = FALSE, pval = FALSE, 
                          risk.table = TRUE, risk.table.y.text.col = TRUE,
                          xlab = "Time (days)", ylab = "p(negativity)") + 
  NULL
surv.plot.3
# where viral load calculated lod usually lod and lod not separately repored 
# with(viral_dynamics, table(vl_quality_cat, lod_value_log10_cpml))
# with(viral_dynamics, table(vl_quality_cat, lod_value_log10_cpml))
cox.mod2.3 <- coxph(Surv(time_to_event, event) ~ vl_quality_cat + sample_site_group_cat, 
                    data = viral_dynamics)
summary(cox.mod2.3)
cox.mod2.3_inter <- coxph(Surv(time_to_event, event) ~ vl_quality_cat * sample_site_group_cat, 
                          data = viral_dynamics)
summary(cox.mod2.3_inter)
anova(cox.mod2.3_inter, cox.mod2.3, test = "Chisq")
# no improvement with interaction term
######## Model 4, viral load lod #########
km.mod4 <- survfit(Surv(time_to_event, event) ~ lod_value_log10_cpml, 
                   data = viral_dynamics)
cox.mod4 <- coxph(Surv(time_to_event, event) ~ lod_value_log10_cpml, 
                  data = viral_dynamics)
summary(cox.mod4)
cox.mod2.4 <- coxph(Surv(time_to_event, event) ~ lod_value_log10_cpml + sample_site_group_cat, 
                    data = viral_dynamics)
summary(cox.mod2.4)
cox.mod2.4_inter <- coxph(Surv(time_to_event, event) ~ lod_value_log10_cpml * sample_site_group_cat, 
                          data = viral_dynamics)
summary(cox.mod2.4_inter)
anova(cox.mod2.4_inter, cox.mod2.4, test = "Chisq")
# no improvement with interaction term
cox.mod3.4 <- coxph(Surv(time_to_event, event) ~ lod_value_log10_cpml + vl_quality_cat, 
                    data = viral_dynamics)
summary(cox.mod2.4)
cox.mod3.4_inter <- coxph(Surv(time_to_event, event) ~ lod_value_log10_cpml * vl_quality_cat, 
                          data = viral_dynamics)
summary(cox.mod3.4_inter)
anova(cox.mod3.4_inter, cox.mod3.4, test = "Chisq")
# interaction better
######## Model 5, age #########
# for plotting use age categories but treat age as continuous in coxph, check again at end 
table(viral_dynamics$age_cat)
#
km.mod5 <- survfit(Surv(time_to_event, event) ~ age_cat, data = viral_dynamics)
cox.mod5 <- coxph(Surv(time_to_event, event) ~ age_y, data = viral_dynamics)
summary(cox.mod5)
# Investigate nonlinearity
cox.mod5.sp2<-  coxph(Surv(time_to_event, event) ~ splines::ns(age_y, df=2), data = viral_dynamics) 
cox.mod5.sp3<-  coxph(Surv(time_to_event, event) ~ splines::ns(age_y, df=3), data = viral_dynamics) 

BIC(cox.mod5, cox.mod5.sp2, cox.mod5.sp3) 
## cox.mod5 is preferable (use BIC, not anova as the models aren't nested)

surv.plot.5 <- ggsurvplot(km.mod5, conf.int = FALSE, pval = FALSE, 
                          risk.table = TRUE, risk.table.y.text.col = TRUE,
                          xlab = "Time (days)", ylab = "p(negativity)") +
  NULL
surv.plot.5
# compare combination and interaction
cox.mod2.5 <- coxph(Surv(time_to_event, event) ~ age_y + sample_site_group_cat, 
                    data = viral_dynamics)
cox.mod2.5_inter <- coxph(Surv(time_to_event, event) ~ age_y * sample_site_group_cat, 
                          data = viral_dynamics)
anova(cox.mod2.5_inter, cox.mod2.5, test = "Chisq")
# no improvement with interaction term
cox.mod3.5 <- coxph(Surv(time_to_event, event) ~ age_y + vl_quality_cat*lod_value_log10_cpml, 
                    data = viral_dynamics)
cox.mod3.5_inter <- coxph(Surv(time_to_event, event) ~ age_y * vl_quality_cat * lod_value_log10_cpml, 
                          data = viral_dynamics)
anova(cox.mod3.5_inter, cox.mod3.5, test = "Chisq")
# no improvement with interaction term
######## Model 6, sex #########
# We will bootstrap 10000 male values using the true ones when known and 
#  randomly generated values based on reported proportions when unknown, and 
#  compare coefficients to when the value is rounded based on proportions
#  i.e. 

km.mod6 <- survfit(Surv(time_to_event, event) ~ male_rounded, data = viral_dynamics)
cox.mod6 <- coxph(Surv(time_to_event, event) ~ male_rounded, data = viral_dynamics)
summary(cox.mod6)
surv.plot.6 <- ggsurvplot(km.mod6, conf.int = FALSE, pval = TRUE, 
                          risk.table = TRUE, risk.table.y.text.col = TRUE,
                          xlab = "Time (days)", ylab = "p(negativity)") +
  NULL
surv.plot.6
boot_out <- data.frame(matrix(nrow = 10000, ncol = 5))
colnames(boot_out) <- colnames(summary(cox.mod5)$coef)
for(i in c(1:10000)){
  set.seed(1 + i)
  male_boot <- rbinom(n = length(viral_dynamics$male), 
                      size = 1, prob = viral_dynamics$male)
  male_boot[male_boot == 1] <- "male"
  male_boot[male_boot == 0] <- "female"
  viral_dynamics$male_boot  <- male_boot
  cox.mod.boot <- coxph(Surv(time_to_event, event) ~ male_boot,
                        data = viral_dynamics)
  boot_out[i,] <- summary(cox.mod.boot)$coeff
}
summary(boot_out$coef)
summary(boot_out$`Pr(>|z|)`)
summary(cox.mod6)
# conclusion - on average the p-value is ~0.5, ie suggests
#  similar that sex alone not significant
# Plan - move forward with male_rounded 
cox.mod2.6 <- coxph(Surv(time_to_event, event) ~ male_rounded + sample_site_group_cat, 
                    data = viral_dynamics)
cox.mod2.6_inter <- coxph(Surv(time_to_event, event) ~ male_rounded * sample_site_group_cat, 
                          data = viral_dynamics)
anova(cox.mod2.6_inter, cox.mod2.6, test = "Chisq")
# no improvement with interaction term
cox.mod3.6 <- coxph(Surv(time_to_event, event) ~ male_rounded + vl_quality_cat * lod_value_log10_cpml, 
                    data = viral_dynamics)
cox.mod3.6_inter <- coxph(Surv(time_to_event, event) ~ male_rounded * 
                            vl_quality_cat * lod_value_log10_cpml, 
                          data = viral_dynamics)
anova(cox.mod3.6_inter, cox.mod3.6, test = "Chisq")
# no improvement with interaction term
cox.mod5.6 <- coxph(Surv(time_to_event, event) ~ male_rounded + age_y, 
                    data = viral_dynamics)
summary(cox.mod5.6)
cox.mod5.6_inter <- coxph(Surv(time_to_event, event) ~ male_rounded * age_y, 
                          data = viral_dynamics)
summary(cox.mod5.6_inter)
anova(cox.mod5.6_inter, cox.mod5.6, test = "Chisq")
# no improvement with interaction term
######## Model 7, disease severity #########
km.mod7 <- survfit(Surv(time_to_event, event) ~ disease_status_cat, 
                   data = viral_dynamics)
cox.mod7 <- coxph(Surv(time_to_event, event) ~ disease_status_cat, 
                  data = viral_dynamics)
summary(cox.mod7)
surv.plot.7 <- ggsurvplot(km.mod7, conf.int = FALSE, pval = TRUE, 
                          risk.table = TRUE, risk.table.y.text.col = TRUE,
                          xlab = "Time (days)", ylab = "p(negativity)")
surv.plot.7
# compare combination and interaction
cox.mod2.7 <- coxph(Surv(time_to_event, event) ~ disease_status_cat + sample_site_group_cat, 
                    data = viral_dynamics)
cox.mod2.7_inter <- coxph(Surv(time_to_event, event) ~ disease_status_cat * sample_site_group_cat, 
                          data = viral_dynamics)
anova(cox.mod2.7_inter, cox.mod2.7, test = "Chisq")
# no improvement with interaction term
cox.mod3.7 <- coxph(Surv(time_to_event, event) ~ disease_status_cat + 
                      lod_value_log10_cpml + vl_quality_cat, 
                    data = viral_dynamics)
cox.mod3.7_inter <- coxph(Surv(time_to_event, event) ~ disease_status_cat * 
                            (vl_quality_cat + lod_value_log10_cpml), 
                          data = viral_dynamics)
anova(cox.mod3.7_inter, cox.mod3.7, test = "Chisq")
# borderline improvement with interaction term, but this is a complex interaction to explain
cox.mod5.7 <- coxph(Surv(time_to_event, event) ~ disease_status_cat + age_y, 
                    data = viral_dynamics)
cox.mod5.7_inter <- coxph(Surv(time_to_event, event) ~ disease_status_cat * age_y, 
                          data = viral_dynamics)
anova(cox.mod5.7_inter, cox.mod5.7, test = "Chisq")
# no improvement with interaction term
cox.mod6.7 <- coxph(Surv(time_to_event, event) ~ disease_status_cat + male_rounded, 
                    data = viral_dynamics)
cox.mod6.7_inter <- coxph(Surv(time_to_event, event) ~ disease_status_cat * male_rounded, 
                          data = viral_dynamics)
anova(cox.mod6.7, cox.mod6.7_inter, test = "Chisq")
# no improvement with interaction term
######## Model 8, multivariable on demographics  #########
# do not use male/disease severity in interaction:
#cox.mod8 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
#                    vl_quality_cat * lod_value_log10_cpml + age_y + 
#                    disease_status_cat * male_rounded + 
#                    drug_any_ever,
#                  data = viral_dynamics)
#summary(cox.mod8)
cox.mod8 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                    vl_quality_cat * lod_value_log10_cpml + age_y + 
                    male_rounded + disease_status_cat,  
                  data = viral_dynamics)
summary(cox.mod8)
cox.mod8_nointer <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                            vl_quality_cat + lod_value_log10_cpml + age_y + 
                            male_rounded + disease_status_cat,  
                          data = viral_dynamics)
summary(cox.mod8_nointer)
anova(cox.mod8, cox.mod8_nointer, test = "Chisq")
# interaction still better on mv
########## model 8 diagnostics ##########
zph.mod8 <- cox.zph(cox.mod8)
print(zph.mod8) ### proportional hazard assumptions tests

myfun = function(p){return(log(-log(p)))}
plot(km.mod7, col=c("black", "green", "blue", "red"), fun="cloglog")
legend(x = 1, y = 1, legend = c("Asymptomatic", "Mild Symptoms",
                                "Moderate Symptoms","Severe Symptoms"),
       lty=rep(1, 4),col=c("black", "green", "blue", "red"), cex = 3/4)
# ph assumption not held for disease severity, likely since asymptomatic very different
# re-run model excluding asymptomatics:
cox.mod8a <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                     vl_quality_cat * lod_value_log10_cpml + age_y + 
                     male_rounded + disease_status_cat,  
                   data = viral_dynamics[viral_dynamics$disease_status_cat != "0=Asymptomatic", ])
summary(cox.mod8a)
zph.mod8a <- cox.zph(cox.mod8a)
print(zph.mod8a) ### proportional hazard assumptions tests
# Does not completely fix as severe also quite different
# remove disease severity
cox.mod8b <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                     vl_quality_cat * lod_value_log10_cpml + age_y + 
                     male_rounded,  
                   data = viral_dynamics)
summary(cox.mod8b)
zph.mod8b <- cox.zph(cox.mod8b)
print(zph.mod8b) ### proportional hazard assumptions tests
anova(cox.mod8, cox.mod8b, test = "Chisq")
# Disease severity improving fit 
#  Conclusion: distinction between severe (ventilated) and moderate (O2 only)
#   considered too great to merge and disease severity offering some imporvement in fit, 
#   so we will accept ph assumption not being met for this covariate and retain 
#   it in the model
plot(zph.mod8)
ggcoxdiagnostics(cox.mod8) ### seems OK
# We now take the multivariable model forward to go drug by drug, looking for
#  further interaction terms only at the drug level
######## model 9 lopinavir #########
with(viral_dynamics, table(drug_lpvr_ever))
with(viral_dynamics, table(drug_riba_ever))
with(viral_dynamics, table(drug_ifn_ever))
with(viral_dynamics, table(drug_ifn_alpha_ever))
with(viral_dynamics, table(drug_ifn_beta_ever))
with(viral_dynamics, table(drug_cqhcq_ever))
with(viral_dynamics, table(drug_umif_ever))
with(viral_dynamics, table(drug_azit_ever))
with(viral_dynamics, table(drug_thym_ever))
with(viral_dynamics, table(drug_remd_ever))
#
cox.mod9 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                    vl_quality_cat * lod_value_log10_cpml + age_y + 
                    male_rounded + disease_status_cat + 
                    drug_lpvr_ever,
                  data = viral_dynamics)
summary(cox.mod9)
forest_model(cox.mod9)
# LPV/r not significant
######## model 10 ribavirin ########
with(viral_dynamics, table(drug_riba_ever, drug_lpvr_ever))
# ribavirin only given without lpvr in 3 patients
cox.mod10 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                     vl_quality_cat * lod_value_log10_cpml + age_y + 
                     male_rounded + disease_status_cat + 
                     drug_riba_ever,
                   data = viral_dynamics)
summary(cox.mod10)
forest_model(cox.mod10)
# ribavirin significant BUT almost never given without LPV/r
cox.mod9.10 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                       vl_quality_cat * lod_value_log10_cpml + age_y + 
                       male_rounded + disease_status_cat + 
                       drug_lpvr_ever + drug_riba_ever,
                     data = viral_dynamics)
summary(cox.mod9.10)
# since hardly any riba given without lpvr almost an interaction term?
cox.mod9.10_inter <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                             vl_quality_cat * lod_value_log10_cpml + age_y + 
                             male_rounded + disease_status_cat + 
                             drug_lpvr_ever * drug_riba_ever,
                           data = viral_dynamics)
summary(cox.mod9.10_inter)
anova(cox.mod9.10_inter, cox.mod9.10, test = "Chisq")
# no improvement with interaction term, but additive is interaction?
######## model 11 interferon ########
with(viral_dynamics, table(drug_ifn_ever, drug_riba_ever, drug_lpvr_ever))
with(viral_dynamics, table(drug_ifn_alpha_ever, drug_riba_ever, drug_lpvr_ever))
with(viral_dynamics, table(drug_ifn_beta_ever, drug_riba_ever, drug_lpvr_ever))
# any interferon
cox.mod11 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                     vl_quality_cat * lod_value_log10_cpml + age_y + 
                     male_rounded + disease_status_cat + 
                     drug_ifn_ever,
                   data = viral_dynamics)
summary(cox.mod11)
cox.mod11a <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                      vl_quality_cat * lod_value_log10_cpml + age_y + 
                      male_rounded + disease_status_cat + 
                      drug_ifn_alpha_ever,
                    data = viral_dynamics)
summary(cox.mod11a)
cox.mod11b <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                      vl_quality_cat * lod_value_log10_cpml + age_y + 
                      male_rounded + disease_status_cat + 
                      drug_ifn_beta_ever,
                    data = viral_dynamics)
summary(cox.mod11b)
# but ifn beta only ever given in triple therapy with riba and lpvr
cox.mod9.11a <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                        vl_quality_cat * lod_value_log10_cpml + age_y + 
                        male_rounded + disease_status_cat + 
                        drug_lpvr_ever + drug_ifn_alpha_ever,
                      data = viral_dynamics)
summary(cox.mod9.11a)
forest_model(cox.mod9.11a)
cox.mod9.11a_inter <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                              vl_quality_cat * lod_value_log10_cpml + age_y + 
                              male_rounded + disease_status_cat + 
                              drug_lpvr_ever * drug_ifn_alpha_ever,
                            data = viral_dynamics)
anova(cox.mod9.11a_inter, cox.mod9.11a, test = "Chisq")
# no interaction
#ifn beta mainly given in triple combination with riba and lpv/r, but both are similar so keep
#  merged interferons 
# 
######## model 12 chloroquine/hcq #########
with(viral_dynamics, table(drug_cqhcq_ever, drug_lpvr_ever))
with(viral_dynamics, table(drug_cqhcq_ever, drug_riba_ever))
with(viral_dynamics, table(drug_cqhcq_ever, drug_ifn_alpha_ever))
# only 1 overlap so no point checking interactions with these
cox.mod12 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                     vl_quality_cat * lod_value_log10_cpml + age_y + 
                     male_rounded + disease_status_cat + 
                     drug_cqhcq_ever,
                   data = viral_dynamics)
summary(cox.mod12)
forest_model(cox.mod12)
# appears to be in wrong direction, cq/hcq worsening time to viral clearance
######## model 13 azithromycin #########
with(viral_dynamics, table(drug_azit_ever, drug_lpvr_ever))
with(viral_dynamics, table(drug_azit_ever, drug_riba_ever))
with(viral_dynamics, table(drug_azit_ever, drug_ifn_alpha_ever))
with(viral_dynamics, table(drug_cqhcq_ever, drug_azit_ever))
# only reasonable overlap so check interactions with these
cox.mod13 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                     vl_quality_cat * lod_value_log10_cpml + age_y + 
                     male_rounded + disease_status_cat + 
                     drug_azit_ever,
                   data = viral_dynamics)
summary(cox.mod13)
# no effect
cox.mod9.13 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                       vl_quality_cat * lod_value_log10_cpml + age_y + 
                       male_rounded + disease_status_cat + 
                       drug_lpvr_ever + drug_azit_ever,
                     data = viral_dynamics)
cox.mod9.13_inter <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                             vl_quality_cat * lod_value_log10_cpml + age_y + 
                             male_rounded + disease_status_cat + 
                             drug_lpvr_ever * drug_azit_ever,
                           data = viral_dynamics)
anova(cox.mod9.13_inter, cox.mod9.13, test = "Chisq")
# no interaction 
cox.mod10.13 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                        vl_quality_cat * lod_value_log10_cpml + age_y + 
                        male_rounded + disease_status_cat + 
                        drug_riba_ever + drug_azit_ever,
                      data = viral_dynamics)
cox.mod10.13_inter <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                              vl_quality_cat * lod_value_log10_cpml + age_y + 
                              male_rounded + disease_status_cat + 
                              drug_riba_ever * drug_azit_ever,
                            data = viral_dynamics)
anova(cox.mod10.13_inter, cox.mod10.13, test = "Chisq")
# no interaction
cox.mod12.13 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                        vl_quality_cat * lod_value_log10_cpml + age_y + 
                        male_rounded + disease_status_cat + 
                        drug_cqhcq_ever + drug_azit_ever,
                      data = viral_dynamics)
cox.mod12.13_inter <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                              vl_quality_cat * lod_value_log10_cpml + age_y + 
                              male_rounded + disease_status_cat + 
                              drug_cqhcq_ever * drug_azit_ever,
                            data = viral_dynamics)
anova(cox.mod12.13_inter, cox.mod12.13, test = "Chisq")
# no interaction
######## model 14 umifenovir (Arbidol) #########
with(viral_dynamics, table(drug_umif_ever))
with(viral_dynamics, table(drug_umif_ever, drug_lpvr_ever))
with(viral_dynamics, table(drug_umif_ever, drug_riba_ever))
with(viral_dynamics, table(drug_umif_ever, drug_ifn_alpha_ever))
with(viral_dynamics, table(drug_umif_ever, drug_ifn_beta_ever))
with(viral_dynamics, table(drug_cqhcq_ever, drug_umif_ever))
with(viral_dynamics, table(drug_azit_ever, drug_umif_ever))
# only test interactions with ifn, lpvr and riba
cox.mod14 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                     vl_quality_cat * lod_value_log10_cpml + age_y + 
                     male_rounded + disease_status_cat + 
                     drug_umif_ever,
                   data = viral_dynamics)
summary(cox.mod14)
# no effect 
cox.mod9.14 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                       vl_quality_cat * lod_value_log10_cpml + age_y + 
                       male_rounded + disease_status_cat + 
                       drug_lpvr_ever + drug_umif_ever,
                     data = viral_dynamics)
cox.mod9.14_inter <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                             vl_quality_cat * lod_value_log10_cpml + age_y + 
                             male_rounded + disease_status_cat + 
                             drug_lpvr_ever * drug_umif_ever,
                           data = viral_dynamics)
anova(cox.mod9.14_inter, cox.mod9.14, test = "Chisq")
# no interaction 
cox.mod10.14 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                        vl_quality_cat * lod_value_log10_cpml + age_y + 
                        male_rounded + disease_status_cat + 
                        drug_riba_ever + drug_umif_ever,
                      data = viral_dynamics)
cox.mod10.14_inter <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                              vl_quality_cat * lod_value_log10_cpml + age_y + 
                              male_rounded + disease_status_cat + 
                              drug_riba_ever * drug_umif_ever,
                            data = viral_dynamics)
anova(cox.mod10.14_inter, cox.mod10.14, test = "Chisq")
# no interaction ? too small numbers anyhow
cox.mod11.14 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                        vl_quality_cat * lod_value_log10_cpml + age_y + 
                        male_rounded + disease_status_cat + 
                        drug_ifn_alpha_ever + drug_umif_ever,
                      data = viral_dynamics)
cox.mod11.14_inter <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                              vl_quality_cat * lod_value_log10_cpml + age_y + 
                              male_rounded + disease_status_cat + 
                              drug_ifn_alpha_ever * drug_umif_ever,
                            data = viral_dynamics)
anova(cox.mod11.14_inter, cox.mod11.14, test = "Chisq")
# ?weak interaction 
######## model 15 thymalfasin #########
with(viral_dynamics, table(drug_thym_ever))
with(viral_dynamics, table(drug_thym_ever, drug_lpvr_ever))
with(viral_dynamics, table(drug_thym_ever, drug_riba_ever))
with(viral_dynamics, table(drug_thym_ever, drug_ifn_alpha_ever))
with(viral_dynamics, table(drug_thym_ever, drug_umif_ever))
with(viral_dynamics, table(drug_azit_ever, drug_thym_ever))
# only test interactions with ifn, umif
cox.mod15 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                     vl_quality_cat * lod_value_log10_cpml + age_y + 
                     male_rounded + disease_status_cat + 
                     drug_thym_ever,
                   data = viral_dynamics)
summary(cox.mod15)
# no effect 
cox.mod14.15 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                        vl_quality_cat * lod_value_log10_cpml + age_y + 
                        male_rounded + disease_status_cat + 
                        drug_umif_ever + drug_thym_ever,
                      data = viral_dynamics)
cox.mod14.15_inter <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                              vl_quality_cat * lod_value_log10_cpml + age_y + 
                              male_rounded + disease_status_cat + 
                              drug_umif_ever * drug_thym_ever,
                            data = viral_dynamics)
anova(cox.mod14.15_inter, cox.mod14.15, test = "Chisq")
# no interaction
######### model 16 Remdesevir ########
with(viral_dynamics, table(drug_remd_ever))
with(viral_dynamics, table(drug_remd_ever, drug_lpvr_ever))
with(viral_dynamics, table(drug_remd_ever, drug_riba_ever))
with(viral_dynamics, table(drug_remd_ever, drug_ifn_alpha_ever))
with(viral_dynamics, table(drug_remd_ever, drug_umif_ever))
with(viral_dynamics, table(drug_remd_ever, drug_thym_ever))
with(viral_dynamics, table(drug_remd_ever, drug_azit_ever))
# always given alone
cox.mod16 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                     vl_quality_cat * lod_value_log10_cpml + age_y + 
                     male_rounded + disease_status_cat + 
                     drug_remd_ever,
                   data = viral_dynamics)
summary(cox.mod16)
# not significant 
########### model 17 multivariable model ##########
# options for interactions:
# 1. no interactions
cox.mod17 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                     vl_quality_cat + lod_value_log10_cpml + age_y + 
                     male_rounded + disease_status_cat + 
                     drug_lpvr_ever + drug_riba_ever + drug_ifn_ever +
                     drug_umif_ever + drug_cqhcq_ever + drug_azit_ever + 
                     drug_thym_ever + 
                     drug_remd_ever,
                   data = viral_dynamics)
summary(cox.mod17)
forest_model(cox.mod17)
cox.mod17ifn <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                     vl_quality_cat + lod_value_log10_cpml + age_y + 
                     male_rounded + disease_status_cat + 
                     drug_lpvr_ever + drug_riba_ever +
                     drug_ifn_alpha_ever + drug_ifn_beta_ever +
                     drug_umif_ever + drug_cqhcq_ever + drug_azit_ever + 
                     drug_thym_ever + 
                     drug_remd_ever,
                   data = viral_dynamics)
summary(cox.mod17ifn)
forest_model(cox.mod17ifn)
with(viral_dynamics, table(drug_riba_ever, drug_ifn_alpha_ever))
with(viral_dynamics, table(drug_riba_ever, drug_ifn_beta_ever))

cox.mod17_inter_both <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                                vl_quality_cat * lod_value_log10_cpml + age_y +
                                male_rounded + disease_status_cat + 
                                (drug_lpvr_ever + drug_riba_ever) * drug_ifn_ever +
                                drug_umif_ever +
                                drug_cqhcq_ever + drug_azit_ever +  drug_thym_ever +
                                drug_remd_ever,
                              data = viral_dynamics)
summary(cox.mod17_inter_both)
cox.mod17_inter_rib <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                               vl_quality_cat * lod_value_log10_cpml + age_y +
                               male_rounded + disease_status_cat + 
                               drug_lpvr_ever + drug_riba_ever * drug_ifn_ever +
                               drug_umif_ever +
                               drug_cqhcq_ever + drug_azit_ever +  drug_thym_ever +
                               drug_remd_ever,
                             data = viral_dynamics)
summary(cox.mod17_inter_rib)
cox.mod17_inter_lpv <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                               vl_quality_cat * lod_value_log10_cpml + age_y +
                               male_rounded + disease_status_cat + 
                               drug_lpvr_ever * drug_riba_ever + drug_ifn_ever +
                               drug_umif_ever +
                               drug_cqhcq_ever + drug_azit_ever +  drug_thym_ever +
                               drug_remd_ever,
                             data = viral_dynamics)
summary(cox.mod17_inter_lpv)
cox.mod17_inter_all <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                               vl_quality_cat * lod_value_log10_cpml + age_y +
                               male_rounded + disease_status_cat + 
                               drug_ifn_ever * drug_lpvr_ever +
                               drug_ifn_ever * drug_riba_ever +
                               drug_lpvr_ever * drug_riba_ever +
                               drug_umif_ever +
                               drug_cqhcq_ever + drug_azit_ever +  drug_thym_ever +
                               drug_remd_ever,
                             data = viral_dynamics)
summary(cox.mod17_inter_all)
forest_model(cox.mod17_inter_all)
anova(cox.mod17, cox.mod17_inter_lpv, test = "Chisq")
# slight improvement
anova(cox.mod17, cox.mod17_inter_rib, test = "Chisq")
# big improvement
anova(cox.mod17_inter_both, cox.mod17_inter_rib, test = "Chisq")
# no improvement
anova(cox.mod17_inter_all, cox.mod17_inter_rib, test = "Chisq")
#  Present double interactions in final model to show relative effects of all combinations
########### model 18 multivariable model refinement ##########
cox.mod18 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat + 
                     vl_quality_cat * lod_value_log10_cpml + age_y + 
                     male_rounded + disease_status_cat + 
                     drug_ifn_ever * drug_lpvr_ever +
                     drug_ifn_ever * drug_riba_ever +
                     drug_lpvr_ever * drug_riba_ever +
                     drug_umif_ever +
                     drug_cqhcq_ever + drug_azit_ever +  drug_thym_ever +
                     drug_remd_ever,
                   data = viral_dynamics)
summary(cox.mod18)
forest_model(cox.mod18)
cox.mod18_agecat <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat +
                            vl_quality_cat * lod_value_log10_cpml + age_cat + 
                            male_rounded + disease_status_cat + 
                            drug_azit_ever + drug_cqhcq_ever + drug_remd_ever +
                            drug_umif_ever +  drug_thym_ever +
                            drug_ifn_ever * drug_lpvr_ever +
                            drug_ifn_ever * drug_riba_ever +
                            drug_lpvr_ever * drug_riba_ever,
                          data = viral_dynamics)
summary(cox.mod18_agecat)
anova(cox.mod18, cox.mod18_agecat, test = "Chisq")
# not much worse fit with categories and this makes visualising effect easier
cox.mod18_lod <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat +
                         vl_quality_cat + lod_value_log10_cpml + 
                         age_cat + male_rounded + disease_status_cat + 
                         drug_azit_ever + drug_cqhcq_ever + drug_remd_ever +
                         drug_umif_ever +  drug_thym_ever +
                         drug_ifn_ever * drug_lpvr_ever +
                         drug_ifn_ever * drug_riba_ever +
                         drug_lpvr_ever * drug_riba_ever,
                       data = viral_dynamics)
summary(cox.mod18_lod)
anova(cox.mod18_lod, cox.mod18_agecat, test = "Chisq")
# Removing lod/quality interation term does not change results in terms of significant effects 
#  and easier to visualise
########### model 19 final multivariable model ##########
cox.mod19 <- coxph(Surv(time_to_event, event) ~ sample_site_group_cat +
                     vl_quality_cat + lod_value_log10_cpml + 
                     age_cat + male_rounded + disease_status_cat + 
                     drug_azit_ever + drug_cqhcq_ever + drug_remd_ever +
                     drug_umif_ever +  drug_thym_ever +
                     drug_ifn_ever * drug_lpvr_ever +
                     drug_ifn_ever * drug_riba_ever +
                     drug_lpvr_ever * drug_riba_ever,
                   data = viral_dynamics)
summary(cox.mod19)
forest_model(cox.mod19)
