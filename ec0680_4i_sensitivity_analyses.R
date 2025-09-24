## EC0680
## Sensitivity analyses
## A. Sherris

source("code/4_analysis/ec0680_4a_prep_for_analysis.R")
source("code/4_analysis/ec0680_4h_prep_for_sensitivity_analyses.R")

load("processed data/2_dataset_creation/MICE/wf_imp_survival.RData")
load("processed data/2_dataset_creation/MICE/wf_imp.RData")
load("processed data/2_dataset_creation/wf_analytic_dataset_survival_deid.RData")
load("processed data/2_dataset_creation/wf_analytic_dataset_deid.RData")


## TRIMESTER-SPECIFIC MODELS (GLM)

res_preterm_trimesters_mice_coadjust <- exposure_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_trimesters_mice",
                         outcome =  "birth_ga_preterm",
                         exposure_period = "Trimester",
                         exposure = x,
                         "covars_primary", wf_imp)) %>% 
  mutate(exposure_period = rep(c("T1","T2","T3"), 14))

save(res_preterm_trimesters_mice_coadjust, file = "output/2_preterm_birth/res_preterm_trimesters_mice.RData")


# Logistic regression / EXPOSURE PERIOD: 32 WEEKS (GLM) 

res_preterm_wk32 <- exposure_list %>%
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_unpooled_mice",
                         outcome = "birth_ga_preterm",
                         exposure_period = "wk32",
                         exposure = x,
                         covars = "covars_primary",
                         data = wf_imp))
save(res_preterm_wk32, file = "output/8_sensitivity/res_preterm_wk32.RData")

load("output/2_preterm_birth/res_survival_preterm_mice.RData")
res_compare <- bind_rows(
  mutate(res_preterm_wk32, model = "Logistic regression (0-32 weeks)"),
  mutate(res_survival_preterm_mice, model = "Pooled logistic regression (survival))"))

plot_compare_models(res_compare)

# TEMPERATURE

res_survival_preterm_temp <- exposure_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice",
                         outcome = "preterm",
                         exposure_period = "surv", 
                         exposure = x,
                         covars = "covars_temp", 
                         data = wf_imp_survival)) %>% 
  relabel_exposure()
save(res_survival_preterm_temp, file = "output/8_sensitivity/preterm_temp.RData")

# TOTAL PM  

res_survival_preterm_totalPM <- exposure_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice",
                         outcome = "preterm",
                         exposure_period = "surv", 
                         exposure = x,
                         covars = "covars_totalPM", 
                         data = wf_imp_survival)) %>% 
  relabel_exposure()
save(res_survival_preterm_totalPM, file = "output/8_sensitivity/preterm_totalPM.RData")


# NO MICE 

res_survival_preterm_nomice <- exposure_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_no_mice",
                         outcome = "preterm",
                         exposure_period = "surv", 
                         exposure = x,
                         covars = "covars_primary", 
                         data = wf_analytic_survival)) %>% 
  relabel_exposure()

save(res_survival_preterm_nomice, file = "output/8_sensitivity/res_survival_preterm_nomice.rData")

# GLM models (cohort fixed effect)

res_survival_preterm_glm_mice <- exposure_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glm_pooled_fixed_effect",
                         outcome = "preterm",
                         exposure_period = "surv", 
                         exposure = x,
                         covars = "covars_cohort", 
                         data = wf_imp_survival))

save(res_survival_preterm_glm_mice, file = "output/8_sensitivity/res_survival_preterm_glm_mice.RData")

# cohort random effect

res_survival_preterm_randomeffect <- exposure_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_random_effect",
                         outcome = "preterm",
                         exposure_period = "surv", 
                         exposure = x,
                         covars = "covars_primary", 
                         data = wf_imp_survival))

save(res_survival_preterm_randomeffect, file = "output/8_sensitivity/res_survival_preterm_randomeffect.RData")

# end
