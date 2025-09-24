## EC0680
## Weekly cumulative exposure metrics
## A. Sherris

# load analytic dataset
load("processed data/2_dataset_creation/MICE/wf_imp.RData")
load("processed data/2_dataset_creation/wf_analytic_dataset_deid.RData")
load("processed data/1_exposure_assessment/wf_pm25_exposure_weekly.RData")

# remove weekly exposure after birth week
wf_exposure_weekly_join <- wf_exposure_weekly %>% 
  left_join(select(wf_analytic, xParticipantID, birth_ga)) %>% 
  filter(gest_week < birth_ga) %>% 
  select(-birth_ga)
  
# add weekly exposure to MICE dataset --------------------------------------

wf_imp_complete <- mice::complete(wf_imp, action = "long", include = T)

wf_imp_weekly <- wf_imp_complete %>% 
  select(-(preg_smokePM_mean:wk32_smoke_wave4_over20)) %>% 
  left_join(wf_exposure_weekly_join, relationship = "many-to-many")  %>% 
  mutate(newid = row_number()) %>% 
  as.mids(.id = "newid")

## PTB WEEKLY ASSOCIATIONS  ---------------------------------------------------------------

week_list <- c(0:35)

# smoke days >0 ug/m3
res_weekly_preterm_smoke_day <- week_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_unpooled_mice", 
                         outcome = "birth_ga_preterm", 
                         exposure_period = "weekly",
                         exposure = "smoke_days", 
                         covars = "covars_primary", 
                         data = mice::filter(wf_imp_weekly, gest_week == x))) %>% 
  mutate(gest_week = week_list)

# smoke days >2.5 ug/m3
res_weekly_preterm_smoke_day_over25 <- week_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_unpooled_mice", 
                         outcome = "birth_ga_preterm", 
                         exposure_period = "weekly",
                         exposure = "smoke_days_over25", 
                         covars = "covars_primary", 
                         data = mice::filter(wf_imp_weekly, gest_week == x))) %>% 
  mutate(gest_week = week_list)

# smoke days >5 ug/m3
res_weekly_preterm_smoke_day_over5 <- week_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_unpooled_mice", 
                         outcome = "birth_ga_preterm", 
                         exposure_period = "weekly",
                         exposure = "smoke_days_over5", 
                         covars = "covars_primary", 
                         data = mice::filter(wf_imp_weekly, gest_week == x))) %>% 
  mutate(gest_week = week_list)

# smoke days >10 ug/m3
res_weekly_preterm_smoke_day_over10 <- week_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_unpooled_mice", 
                         outcome = "birth_ga_preterm", 
                         exposure_period = "weekly",
                         exposure = "smoke_days_over10", 
                         covars = "covars_primary", 
                         data = mice::filter(wf_imp_weekly, gest_week == x))) %>% 
  mutate(gest_week = week_list)


# combine results
res_weekly_preterm <- bind_rows(res_weekly_preterm_smoke_day,
                             res_weekly_preterm_smoke_day_over25,
                             res_weekly_preterm_smoke_day_over5,
                             res_weekly_preterm_smoke_day_over10)

save(res_weekly_preterm, file = "output/4_weekly_exposure/res_weekly_preterm.RData")

# end
