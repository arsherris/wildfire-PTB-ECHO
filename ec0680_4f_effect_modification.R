## EC0680
## Effect modification / interaction analyses
## Start date: 10/1/23
## A. Sherris

# prep for analysis
source("//echofile.rti.ns/EC0680/USERS/asherris/ANALYSIS/code/0_EC0680_config.R")
source("code/4_analysis/4a_prep_for_analysis.R")

# load imputed analytic datasets
load("processed data/2_dataset_creation/MICE/wf_imp.RData")
load("processed data/2_dataset_creation/MICE/wf_imp_survival.RData")

# define exposure list for EM
exposure_list_em <- c("smokePM_mean", "smoke_days")

# define model for effect modification analysis
glmer_em_mice_pooled <- function(modifier, ref_level, outcome, exposure_period, exposure, covars, data){
  
  
  data_complete <- mice::complete(data, action = "long", include = T)
  data_complete[[modifier]] <- relevel(as.factor(data_complete[[modifier]]), ref = ref_level)
  data_revised <- mice::as.mids(data_complete)
  
  with(data_revised, glmer(formula(paste0(outcome, " ~ ", exposure_period, "_", exposure, "*", modifier, 
                                  " + factor(start_wk)", covars, "+ (1|CohortID)")),
                   family = "binomial", 
                   control = glmer_control,
                   nAGQ = 0))
  
}

# function to extract and clean model output for smoke metric -------------------------
# mod_fxn, outcome, exposure period, exposure, covars as string
em_extract_mod_output <- function(modifier, ref_level, mod_fxn, outcome, exposure_period, exposure, covars, data){
  
  mod_fxn_assigned <- get(mod_fxn)
  covars_assigned <- get(covars)
  
  mod <- mod_fxn_assigned(modifier, ref_level, outcome, exposure_period, exposure, covars_assigned, data)
  
  # pool results of imputed models
  if (data.class(mod) == "mira") {mod_tidy <- pool(mod)} else {mod_tidy <- mod}
  
  # generate model output
  mod_out <- broom.mixed::tidy(mod_tidy) %>% 
    filter(grepl("smoke|:", term)) %>% 
    mutate(outcome = outcome,
           ref_level= ref_level,
           interaction = grepl(":", term),
           exposure_period = exposure_period,
           exposure = str_remove(term, paste0(exposure_period, "_")),
           covars = covars,
           result = estimate, 
           lower = estimate + qnorm(0.025)*std.error,
           upper = estimate + qnorm(0.975)*std.error,
           fxn = mod_fxn,
           mice = data.class(mod) == "mira") %>%
    select(outcome:mice, p_value = p.value) 
  
  # exponentiate results of binomial models
  if (grepl("glm", mod_fxn)) 
    return(mutate_at(mod_out, vars(result:upper), exp))
  else return(mod_out)
  
}


# categorical interactions - wald test for p-value
em_extract_mod_output_categorical <- function(modifier, ref_level, mod_fxn,
                                      outcome, exposure_period, exposure, covars, data){
  
  mod_fxn_assigned <- get(mod_fxn)
  covars_assigned <- get(covars)
  
  mod <- mod_fxn_assigned(modifier, ref_level, outcome, exposure_period, exposure, covars_assigned, data)
  
  mod_pooled <- pool(mod) 

  # get VCOV matrix from pooled model
  m <- mod_pooled$m
  ubar <- Reduce("+", lapply(mod$analyses, vcov)) / (m)
  b <- mod_pooled$pooled$b
  
  t <- ubar + (1 + 1 / (m)) * b
  
  all.equal(as.numeric(diag(t)), mod_pooled$pooled$t)
  
  # get index of interaction terms
  int_index <- which(grepl(":", mod_pooled$pooled$term))
  
  # wald test for interaction global p-value
  wald_out <- wald.test(Sigma = t,
                        b = mod_pooled$pooled$estimate, 
                        Terms = int_index)
  
  # generate model output
  mod_out <- broom.mixed::tidy(mod_pooled) %>% 
    filter(grepl("smoke|:", term)) %>% 
    mutate(outcome = outcome,
           ref_level= ref_level,
           interaction = grepl(":", term),
           exposure_period = exposure_period,
           exposure = str_remove(term, paste0(exposure_period, "_")),
           covars = covars,
           result = estimate, 
           lower = estimate + qnorm(0.025)*std.error,
           upper = estimate + qnorm(0.975)*std.error,
           fxn = mod_fxn,
           mice = data.class(mod) == "mira",
           wald_pval = wald_out$result$chi2[[3]]) %>%
    select(outcome:wald_pval, p_value = p.value) 
  
  # exponentiate results of binomial models
  if (grepl("glm", mod_fxn)) 
    return(mutate_at(mod_out, vars(result:upper), exp))
  else return(mod_out)

}


# PTB -------------------------------------------------------------------------------

# EM by region
em_region_preterm_northeast <- exposure_list_em %>% 
  map_df(\(x) em_extract_mod_output_categorical(modifier = "census_region",
                            ref_level = "Northeast",
                            mod_fxn = "glmer_em_mice_pooled",
                             outcome = "preterm",
                             exposure_period = "surv", 
                             exposure = x,
                             covars = "covars_primary", 
                             data = wf_imp_survival)) 

save(em_region_preterm_northeast, file = "output/6_effect_mod/em_region_preterm_northeast.RData")


# EM by sex
em_sex_preterm <- exposure_list_em %>% 
  map_df(\(x) em_extract_mod_output(modifier = "child_sex",
                            mod_fxn = "glmer_em_mice_pooled",
                            ref_level = "0",
                            outcome = "preterm",
                            exposure_period = "surv", 
                            exposure = x,
                            covars = "covars_primary", 
                            data = wf_imp_survival)) 
save(em_sex_preterm, file = "output/6_effect_mod/em_sex_preterm.RData")

# em by race
table(wf_imp_survival$data$mat_hispanic)

em_race_preterm <- exposure_list_em%>% 
  map_df(\(x) em_extract_mod_output_categorical(modifier = "mat_race_cat",
                            ref_level = "1 White",
                            mod_fxn = "glmer_em_mice_pooled",
                            outcome = "preterm",
                            exposure_period = "surv", 
                            exposure = x,
                            covars = "covars_primary", 
                            data = wf_imp_survival)) 
save(em_race_preterm, file = "output/6_effect_mod/em_race_preterm.RData")

# em by poverty
em_pov_preterm <- exposure_list_em %>% 
  map_df(\(x) em_extract_mod_output_categorical(modifier = "pov",
                            ref_level = "1",
                            mod_fxn = "glmer_em_mice_pooled",
                            outcome = "preterm",
                            exposure_period = "surv", 
                            exposure = x,
                            covars = "covars_primary", 
                            data = wf_imp_survival)) 

save(em_pov_preterm, file = "output/6_effect_mod/em_pov_preterm.RData")





# test f-test code and function
# 
# library(lmtest)
# library(sandwich)
# library(aod)
# 
# mice_mod <- lmer_em_mice("census_region", "Northeast", "birth_bw",
#                            "term", "smokePM_mean", covars_primary, wf_imp_term)
# mice_mod <- mod
# mod_pooled <- mice::pool(mice_mod)
# mod_pooled$pooled$estimate
# 
# 
# m <- mod_pooled$m
# ubar <- Reduce("+", lapply(mice_mod$analyses, vcov)) / (m)
# b <- mod_pooled$pooled$b
# 
# t <- ubar + (1 + 1 / (m)) * b
# 
# all.equal(as.numeric(diag(t)), mod_pooled$pooled$t)
# 
# wald_out <- wald.test(Sigma = t,
#                      b = mod_pooled$pooled$estimate,
#                      Terms = 32:34)
# 
# # test the function to do the same thing
# test = em_extract_mod_output_categorical(modifier = "census_region", ref_level = "Northeast",
#                                  mod_fxn =  "lmer_em_mice",outcome =  "birth_bw",
#                                  exposure_period = "term", exposure =  "smokePM_mean",
#                                  covars = "covars_primary", data =wf_imp_term)

