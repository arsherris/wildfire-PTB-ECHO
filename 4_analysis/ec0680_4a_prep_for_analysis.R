## EC0680
## Analysis preparation and functions
## Start date: 8/27/23
## A. Sherris

load("processed data/2_dataset_creation/wf_analytic_dataset_deid.RData")

library(mice)
library(lme4)
library(splines)
library(broom.mixed)

# load exposure list
exposure_list <-  c("smokePM_mean",
                    "smoke_days",
                    "smoke_days_over25", "smoke_days_over5", "smoke_days_over10",
                    "smoke_wave2_over25", "smoke_wave3_over25", "smoke_wave4_over25",
                    "smoke_wave2_over5", "smoke_wave3_over5", "smoke_wave4_over5",
                    "smoke_wave2_over10", "smoke_wave3_over10", "smoke_wave4_over10")


# define covariates ----------------------------------------------------------

covars_primary  <- " + ns(mat_age, df = 3) + mat_race_cat + mat_hispanic + child_sex + ep_poverty + ns(as.numeric(birth_year), df = 4) + conception_season + X1+X2+X3+X4+X5+X6+X7+X8+X9"
covars_extended <- paste0(covars_primary, "+ parity_cat + sesIMP_educ5*sesIMP_source + birth_ga_src_cat + anthr_bmi_pp + preg_tobacco + preg_alcohol") 

# define models ------------------------------------------------------------

# model controls to allow faster run time
glmer_control <- glmerControl(optimizer = "nloptwrap", calc.derivs = F)
lmer_control  <- lmerControl(optimizer = "nloptwrap", calc.derivs = F)

# primary PLR MICE model

glmer_pooled_mice <- function(outcome, exposure_period, exposure, covars, data){
  
  with(data, glmer(formula(paste0(outcome, " ~ ", exposure_period, "_", exposure, " + factor(start_wk)", covars, "+ (1|CohortID)")),
                   family = "binomial", 
                   control = glmer_control,
                   nAGQ = 0)
  )
}


# unpooled logistic regression (e.g. for weekly analyses)
glmer_unpooled_mice <- function(outcome, exposure_period, exposure, covars, data){
  
  with(data, glmer(formula(paste0(outcome, " ~ ", exposure_period, "_", exposure, covars, "+ (1|CohortID)")),
                   family = "binomial", 
                   control = glmer_control,
                   nAGQ = 0)
  )
}


# LME models - for continuous outcomes (i.e. gestational age) 

lmer_primary_mice <- function(outcome, exposure_period, exposure, covars, data){
  
  with(data, 
       lmer(formula(paste0(outcome, " ~ ", exposure_period, "_", exposure, covars, "+ (1|CohortID)")),
            control = lmer_control))
  
}


# clean results output for plotting -------------------------------------------------

relabel_exposure <- function(data) {
  
  data %>% 
    mutate(exposure = fct_inorder(exposure),
           exposure_group =
             case_when(grepl("wave", exposure) ~ "Waves",
                       grepl("days", exposure) ~ "Days",
                       grepl("mean", exposure) ~ "Overall",
                       T ~ NA_character_),
           duration = 
             case_when(grepl("wave2", exposure) ~ "2 days",
                       grepl("wave3", exposure) ~ "3 days",
                       grepl("wave4", exposure) ~ "4 days",
                       T ~ NA_character_),
           intensity = fct_inorder(case_when(
             grepl("over25", exposure)  ~ "\u2265 2.5",
             grepl("over5", exposure) ~ "\u2265 5.0",
             grepl("over10", exposure) ~ "\u2265 10",
             exposure == "smoke_days" ~ "Any",
             T ~ NA_character_)))%>%  
    return()
  
}


# function to extract and clean model output for smoke metric -------------------------
# mod_fxn, outcome, exposure period, exposure, covars as string

extract_mod_output <- function(mod_fxn, outcome, exposure_period, exposure, covars, data){
  
  # assign function and covariates from string
  
  mod_fxn_assigned <- get(mod_fxn)
  covars_assigned <- get(covars)
  
  mod <- mod_fxn_assigned(outcome, exposure_period, exposure, covars_assigned, data)
  
  # pool results of imputed models
  if (data.class(mod) == "mira") {mod_tidy <- pool(mod)} else {mod_tidy <- mod}
  
  # generate model output
  mod_out <- broom.mixed::tidy(mod_tidy) %>% 
    filter(grepl("smoke", term)) %>% 
    mutate(outcome = outcome,
           exposure_period = exposure_period,
           exposure = exposure,
           covars = covars,
           result = estimate, 
           lower = estimate + qnorm(0.025)*std.error,
           upper = estimate + qnorm(0.975)*std.error,
           fxn = mod_fxn,
           mice = data.class(mod) == "mira") %>%
    select(term,outcome:mice) %>% 
    relabel_exposure
  
  # exponentiate results of binomial models
  if (grepl("glm", mod_fxn)) 
    return(mutate_at(mod_out, vars(result:upper), exp))
  else return(mod_out)
  
}


# plot results -----------------------------------------------------------------

library(patchwork)

plot_res_single_period <- function(data, plot_title) {
  
  y_axis_label <- case_when(data$outcome[1] %in% c("birth_ga_preterm", "preterm", 
                                                 "lbw", "term_lbw", "early_term") ~ "OR and 95% CI",
                          data$outcome[1] %in% c("birth_ga", "birth_bw", "birth_bwz_sex") ~ "Coefficient and 95% CI",
                          T ~ "Unknown")
  
  yint <- if_else(data$outcome[1] %in% c("birth_bw", "birth_ga", 
                                         "birth_sex_bwz"), 0, 1) 
  
  plot_overall <- ggplot(dplyr::filter(data, exposure_group == "Overall"), 
                         aes(x = exposure, y = result))+
    geom_hline(yintercept = yint, col = "grey")+
    geom_point()+
    geom_linerange(aes(ymin = lower, ymax = upper))+
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          text = element_text(size = 12),
          title = element_text(size = 11))+
    # scale_y_continuous(transform = "log")+
    labs(title = expression(Mean ~ smoke ~ PM[2.5]),
         x = "", y = y_axis_label)
  
  plot_days <- ggplot(dplyr::filter(data, exposure_group == "Days"), 
                      aes(x = intensity, y = result))+
    geom_hline(yintercept = yint, col = "grey")+
    geom_point(position = position_dodge(width = 0.7))+
    geom_linerange(aes(ymin = lower, ymax = upper), 
                  position = position_dodge(width = 0.7))+
    theme_classic() +
    theme(text = element_text(size = 12),
          title = element_text(size = 11), 
          axis.title.y = element_blank())+
    # scale_y_continuous(transform = "log")+
    labs(x = expression("Concentration threshold ("  * mu * g ~ m^{-3} *")"), 
         title = "Smoke days")
  
  
  plot_waves <- ggplot(filter(data, exposure_group == "Waves"), 
                       aes(x = intensity, col = duration, y = result))+
    geom_hline(yintercept = yint, col = "grey")+
    geom_point(position = position_dodge(width = 0.7))+
    geom_linerange(aes(ymin = lower, ymax = upper), 
                  position = position_dodge(width = 0.7))+
    scale_color_viridis_d(option = "A", begin = 0.7, end = 0.2)+
    theme_classic() +
    theme(text = element_text(size = 12),
          title = element_text(size = 11), 
          axis.title.y = element_blank())+
    # scale_y_continuous(transform = "log")+
    labs(x = "", #expression(Over ~ 5 ~ mu * g ~ m^{-3}), 
         col = "Wave duration",
         title = "Smoke waves")
  
  plot_all <- plot_overall + plot_days + plot_waves + 
    plot_layout(ncol = 3, widths = c(2.7, 4, 4)) +
    plot_annotation(title = plot_title, 
                    theme = theme(plot.title = element_text(face = "bold", size = 16)))
  return(plot_all)

}


# plot stratified results

plot_stratified <- function(data, strat_variable) {
  
  y_axis_label <- case_when(data$outcome[1] %in% c("birth_ga_preterm", "preterm", 
                                                   "lbw", "term_lbw", "early_term") ~ "OR and 95% CI",
                            data$outcome[1] %in% c("birth_ga", "birth_bw", "birth_bwz_sex") ~ "Coefficient and 95% CI",
                            T ~ "Unknown")
  
  yint <- if_else(grepl("preterm|lbw|early_term", data$outcome[1]), 1, 0) 
  
  data$strat_var <- data[[strat_variable]]
  
  plot_overall <- ggplot(dplyr::filter(data, exposure_group == "Overall"), 
                         aes(x = exposure, y = result))+
    facet_wrap(~strat_var, nrow = 1)+
    geom_hline(yintercept = yint, col = "grey")+
    geom_point()+
    geom_linerange(aes(ymin = lower, ymax = upper))+
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          text = element_text(size = 12),
          title = element_text(size = 11))+
    labs(title = expression(A. ~ Mean ~ smoke ~ PM[2.5]),
         x = "", y = y_axis_label)
  
  plot_days <- ggplot(dplyr::filter(data, exposure_group == "Days"), 
                      aes(x = intensity, y = result))+
    facet_wrap(~strat_var, nrow = 1)+
    geom_hline(yintercept = yint, col = "grey")+
    geom_point(position = position_dodge(width = 0.7))+
    geom_linerange(aes(ymin = lower, ymax = upper), 
                  position = position_dodge(width = 0.7))+
    theme_classic() +
    theme(text = element_text(size = 12),
          title = element_text(size = 11))+
    labs(x = "",
         y = "",
         title = "B. Smoke days")
  
  
  plot_waves <- ggplot(filter(data, exposure_group == "Waves"),
                       aes(x = intensity, col = duration, y = result))+
    facet_wrap(~strat_var, nrow = 1)+
    geom_hline(yintercept = yint, col = "grey")+
    geom_point(position = position_dodge(width = 0.7))+
    geom_linerange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.7))+
    scale_color_viridis_d(option = "A", begin = 0.7, end = 0.2)+
    theme_classic() +
    theme(text = element_text(size = 12),
          title = element_text(size = 11))+
    labs(x =expression("Concentration threshold ("  * mu * g ~ m^{-3} *")"), 
         y = "", col = "Wave duration",
         title = "C. Smoke waves")
  
  plot_all <- plot_overall + plot_days +plot_waves + plot_layout(ncol = 1)
  return(plot_all)
  
}


# plot to compare models -----------------------------------------------------

plot_compare_models <- function(data) {
  
  y_axis_label <- case_when(data$outcome[1] %in% c("birth_ga_preterm", "preterm", 
                                                   "lbw", "term_lbw") ~ "OR and 95% CI",
                            data$outcome[1] %in% c("birth_ga", "birth_bw", "birth_bwz_sex") ~ "Coefficient and 95% CI",
                            T ~ "Unknown")
  
  yint <- if_else(grepl("preterm|lbw", data$outcome[1]), 1, 0) 
  
  
  plot_overall <- ggplot(dplyr::filter(data, exposure_group == "Overall"), 
                         aes(x = exposure, y = result, shape = model))+
    geom_hline(yintercept = yint, col = "grey")+
    geom_point(position = position_dodge(width = 0.3))+
    geom_linerange(aes(ymin = lower, ymax = upper), width = 0.2,
                   position = position_dodge(width = 0.3))+
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          text = element_text(size = 12),
          title = element_text(size = 11),
          legend.position = "none")+
    labs(title = expression(Mean ~ smoke ~ PM[2.5]),
         x = "", y = y_axis_label)+
    scale_shape_manual(values = c(16, 17, 15, 8, 1, 2, 0))
  
  plot_days <- ggplot(dplyr::filter(data, exposure_group == "Days"), 
                      aes(x = intensity, y = result, shape = model))+
    geom_hline(yintercept = yint, col = "grey")+
    geom_point(position = position_dodge(width = 0.7))+
    geom_linerange(aes(ymin = lower, ymax = upper), 
                   position = position_dodge(width = 0.7))+
    theme_classic() +
    theme(text = element_text(size = 12),
          title = element_text(size = 11),
          legend.position = "none", 
          axis.title.y = element_blank())+
    labs(x =expression("Concentration threshold ("  * mu * g ~ m^{-3} *")"), 
         title = "Smoke days")+
    scale_shape_manual(values = c(16, 17, 15, 8, 1, 2, 0))
  
  
  plot_waves <- ggplot(filter(data, exposure_group == "Waves"),
                       aes(x = intensity, col = duration, y = result, shape = model))+
    geom_hline(yintercept = yint, col = "grey")+
    geom_point(position = position_dodge(width = 0.9), size = 1)+
    geom_linerange(aes(ymin = lower, ymax = upper),
                   position = position_dodge(width = 0.9))+
    scale_color_viridis_d(option = "A", begin = 0.7, end = 0.2)+
    theme_classic() +
    theme(text = element_text(size = 12),
          title = element_text(size = 10), 
          axis.title.y = element_blank())+
    labs(x = "",
         col = "Wave duration",
         title = "Smoke waves",
         shape = "Model")+
    scale_shape_manual(values = c(16, 17, 15, 8, 1, 2, 0))
  
  plot_all <- plot_overall + plot_days + plot_waves + plot_layout(ncol = 3, widths = c(2, 3.5, 4.5))
  return(plot_all)
  
}

# end
