## EC0680
# ANALYSIS: WILDFIRE EXPOSURE AND PRETERM BIRTH
## A. Sherris

# prep for analysis
source("code/4_analysis/ec0680_4a_prep_for_analysis.R")

# load imputed dataset for survival analysis
load("processed data/2_dataset_creation/MICE/wf_imp_survival.RData")

# NATIONWIDE PTB SURVIVAL MODELS -------------------------------------------------------------------

# Primary covariates
Sys.time()
res_survival_preterm_mice <- exposure_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice",
                         outcome = "preterm",
                         exposure_period = "surv", 
                         exposure = x,
                         covars = "covars_primary", 
                         data = wf_imp_survival)) 

Sys.time()
save(res_survival_preterm_mice, file = "output/2_preterm_birth_nationwide/res_survival_preterm_mice.RData")
plot_res_single_period(res_survival_preterm_mice)


# Extended covariates 
Sys.time()
res_survival_preterm_extended <- exposure_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice",
                         outcome = "preterm",
                         exposure_period = "surv", 
                         exposure = x,
                         covars = "covars_extended", 
                         data = wf_imp_survival))
Sys.time()
save(res_survival_preterm_extended, file = "output/2_preterm_birth_nationwide/res_survival_preterm_extended.RData")

plot_res_single_period(res_survival_preterm_extended)

# SENSITIVITY ANALYSIS: primary covars - restricted population 

extended_ids <- wf_analytic$xParticipantID[wf_analytic$extended_pop==1]
wf_imp_survival_restricted <- mice::filter(wf_imp_survival, xParticipantID %in% extended_ids)

res_survival_preterm_restricted <- exposure_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice",
                         outcome = "preterm",
                         exposure_period = "surv", 
                         exposure = x,
                         covars = "covars_primary", 
                         data = wf_imp_survival_restricted)) 

save(res_survival_preterm_restricted, file = "output/2_preterm_birth_nationwide/res_survival_preterm_restricted.RData")
plot_res_single_period(res_survival_preterm_restricted)

# compare models
res_survival_preterm_all <- bind_rows(
  mutate(res_survival_preterm_glmer_mice, model = "Primary"),
  mutate(res_survival_preterm_extended, model = "Extended"),
  mutate(res_survival_preterm_restricted, model = "Primary/restricted")
)

plot_compare_models(res_survival_preterm_all)

# end
