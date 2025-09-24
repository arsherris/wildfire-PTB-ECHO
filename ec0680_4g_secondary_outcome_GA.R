## EC0680
## Secondary outcomes
## A. Sherris

## OUTCOME: Gestational age ---------------------------------------------------------------

load("processed data/2_dataset_creation/MICE/wf_imp.RData")

# overall pregnancy (32 weeks)
res_ga_wk32_mice <- exposure_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "lmer_primary_mice", 
                         outcome = "birth_ga", 
                         exposure_period = "wk32",
                         exposure = x, 
                         covars = "covars_primary", 
                         data = wf_imp))

save(res_ga_wk32_mice, file = "output/6_secondary_outcomes/res_ga_wk32_mice.RData")
plot_res_single_period(res_ga_wk32_mice)
Sys.time()

