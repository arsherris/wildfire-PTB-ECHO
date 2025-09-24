## EC0680
## Stratified analyses
## Start date: 10/1/23
## A. Sherris

# prep for analysis
source("code/4_analysis/ec0680_4a_prep_for_analysis.R")

# load imputed analytic datasets
load("processed data/2_dataset_creation/MICE/wf_imp_survival.RData")

# WESTERN US PTB SURVIVAL MODELS -------------------------------------------------------------------

wf_imp_west_surv <- mice::filter(wf_imp_survival, census_region == "West")

# primary covariates
res_survival_preterm_west <- exposure_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice", 
                         outcome = "preterm", 
                         exposure_period = "surv", 
                         exposure = x, 
                         covars = "covars_primary", 
                         data = wf_imp_west_surv))

save(res_survival_preterm_west, file = "output/3_preterm_birth_western/res_survival_preterm_west.RData")


# extended covariates
res_survival_preterm_west_ext <- exposure_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice", 
                         outcome = "preterm", 
                         exposure_period = "surv", 
                         exposure = x, 
                         covars = "covars_extended", 
                         data = wf_imp_west_surv))

res_survival_preterm_west_ext
save(res_survival_preterm_west_ext, file = "output/3_preterm_birth_western/res_survival_preterm_west_ext.RData")

# primary covariates - restricted sample
extended_ids <- wf_analytic$xParticipantID[wf_analytic$extended_pop==1]
wf_imp_west_restricted <- mice::filter(wf_imp_west_surv, xParticipantID %in% extended_ids)

res_survival_preterm_west_restricted <- exposure_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice",
                         outcome = "preterm",
                         exposure_period = "surv", 
                         exposure = x,
                         covars = "covars_primary", 
                         data = wf_imp_west_restricted)) 

save(res_survival_preterm_west_restricted, file = "output/3_preterm_birth_westernUS/res_survival_preterm_west_restricted.Rdata")

# end
