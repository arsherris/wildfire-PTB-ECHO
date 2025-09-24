## EC0680
## Stratified analyses
## Start date: 10/1/23
## A. Sherris

# prep for analysis
source("code/4_analysis/ec0680_4a_prep_for_analysis.R")

# load imputed analytic datasets
load("processed data/2_dataset_creation/MICE/wf_imp.RData")
load("processed data/2_dataset_creation/MICE/wf_imp_survival.RData")

exposure_list_em <- c("smokePM_mean", "smoke_days")

## CENSUS REGIONS ------------------------------------------------

# stratify imputed dataset
wf_imp_midwest_surv <- mice::filter(wf_imp_survival, census_region == "Midwest")
wf_imp_northeast_surv <- mice::filter(wf_imp_survival, census_region == "Northeast")
wf_imp_south_surv <- mice::filter(wf_imp_survival, census_region == "South")

res_survival_preterm_midwest <- exposure_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice", outcome = "preterm", exposure_period = "surv", 
                         exposure = x, covars = "covars_primary", 
                         data = wf_imp_midwest_surv))

res_survival_preterm_northeast <- exposure_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice", outcome = "preterm", exposure_period = "surv", 
                         exposure = x, covars = "covars_primary", 
                         data = wf_imp_northeast_surv))

res_survival_preterm_south <- exposure_list %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice", outcome = "preterm", exposure_period = "surv", 
                         exposure = x, covars = "covars_primary", 
                         data = wf_imp_south_surv))

load("output/3_preterm_birth_western/res_survival_preterm_west.RData")

res_preterm_strat_region <- bind_rows(
  mutate(res_survival_preterm_west, region = "West"),
  mutate(res_survival_preterm_midwest, region = "Midwest"),
  mutate(res_survival_preterm_northeast, region = "Northeast"),
  mutate(res_survival_preterm_south, region = "South")
)

save(res_preterm_strat_region, file = "output/5_stratification/res_preterm_strat_region.RData")


## CHILD SEX --------------------------------------------------

# stratify imputed dataset 
wf_imp_survival_female <- mice::filter(wf_imp_survival, child_sex == 0)
wf_imp_survival_male <- mice::filter(wf_imp_survival, child_sex == 1)

res_survival_preterm_female <- exposure_list_em %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice",
                         outcome = "preterm",
                         exposure_period = "surv", 
                         exposure = x,
                         covars = "covars_primary", 
                         data = wf_imp_survival_female))

res_survival_preterm_male <- exposure_list_em %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice",
                         outcome = "preterm",
                         exposure_period = "surv", 
                         exposure = x,
                         covars = "covars_primary", 
                         data = wf_imp_survival_male))
Sys.time()

res_preterm_strat_sex <- bind_rows(
  mutate(res_survival_preterm_female, sex = "Female"),
  mutate(res_survival_preterm_male, sex = "Male")
)

save(res_preterm_strat_sex, file = "output/5_stratification/res_preterm_strat_sex.RData")

## POVERTY --------------------------------------------------------------------------

# stratify imputed dataset
wf_imp_pov1 <- mice::filter(wf_imp_survival, pov==1)
wf_imp_pov2 <- mice::filter(wf_imp_survival, pov==2)
wf_imp_pov3 <- mice::filter(wf_imp_survival, pov==3)

res_survival_preterm_pov1 <- exposure_list_em %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice", outcome = "preterm", exposure_period = "surv", 
                         exposure = x, covars = "covars_primary", 
                         data = wf_imp_pov1))

res_survival_preterm_pov2<- exposure_list_em %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice", outcome = "preterm", exposure_period = "surv", 
                         exposure = x, covars = "covars_primary", 
                         data = wf_imp_pov2))

res_survival_preterm_pov3 <- exposure_list_em %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice", outcome = "preterm", exposure_period = "surv", 
                         exposure = x, covars = "covars_primary", 
                         data = wf_imp_pov3))

res_preterm_strat_pov <- bind_rows(
  mutate(res_survival_preterm_pov1, strat_var = "First"),
  mutate(res_survival_preterm_pov2, strat_var = "Second"),
  mutate(res_survival_preterm_pov3, strat_var = "Third")
) %>% 
  mutate(modifier = "Poverty tertile")

save(res_preterm_strat_pov, file = "output/5_stratification/res_preterm_strat_pov.RData")
plot_stratified_overall(res_preterm_strat_pov, "strat_var")

## MATERNAL RACE - --------------------------------------------------------------------------

# specify covariate list omitting race
covars_no_race  <- " + ns(mat_age, df = 3) + mat_hispanic + child_sex + ep_poverty + ns(as.numeric(birth_year), df = 4) + conception_season + X1+X2+X3+X4+X5+X6+X7+X8+X9"

# stratify imputed dataset
wf_imp_white <- mice::filter(wf_imp_survival, mat_race_cat == "1 White")
wf_imp_black <- mice::filter(wf_imp_survival, mat_race_cat == "2 Black")
wf_imp_asian <- mice::filter(wf_imp_survival, mat_race_cat == "3 Asian / Pacific Islander")

res_survival_preterm_white <- exposure_list_em %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice", outcome = "preterm", exposure_period = "surv", 
                         exposure = x, covars = "covars_no_race", 
                         data = wf_imp_white))

res_survival_preterm_black<- exposure_list_em %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice", outcome = "preterm", exposure_period = "surv", 
                         exposure = x, covars = "covars_no_race", 
                         data = wf_imp_black))

res_survival_preterm_asian <- exposure_list_em %>% 
  map_df(\(x) extract_mod_output(mod_fxn = "glmer_pooled_mice", outcome = "preterm", exposure_period = "surv", 
                         exposure = x, covars = "covars_no_race", 
                         data = wf_imp_asian))

res_preterm_strat_race <- bind_rows(
  mutate(res_survival_preterm_white, strat_var = "White"),
  mutate(res_survival_preterm_black, strat_var = "Black"),
  mutate(res_survival_preterm_asian, strat_var = "Asian/Pacific Islander")
) %>% 
  mutate(modifier = "Maternal race")

save(res_preterm_strat_race, file = "output/5_stratification/res_preterm_strat_race.RData")
plot_stratified_overall(res_preterm_strat_race, "strat_var")

