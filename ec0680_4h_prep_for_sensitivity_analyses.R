## EC0680
## Prep for sensitivity analyses
## A. Sherris

# define covariates for sensitivity analyses

covars_temp     <- paste0(covars_primary, "+ preg_mean_temp")
covars_totalPM  <- paste0(covars_primary, "+ preg_totalPM_mean")
covars_cohort  <- paste0(covars_primary, "+ CohortID")

# Sensitivity analyses models -----------------------------------------------------------

# primary unpooled trimester-adjusted model

glmer_trimesters_mice <- function(outcome, exposure_period, exposure, covars, data){
  
  with(data, glmer(formula(paste0(outcome, " ~ ", 
                                  " T1_", exposure,
                                  " + T2_", exposure,
                                  " + T3_", exposure, 
                                  covars, "+ (1|CohortID)")),
                   family = "binomial",
                   control = glmer_control,
                   nAGQ = 0))
}

# no MICE
glmer_pooled_no_mice <- function(outcome, exposure_period, exposure, covars, data){
  
  glmer(formula(paste0(outcome, " ~ ", exposure_period, "_", exposure, " + factor(start_wk)", covars, "+ (1|CohortID)")),
        family = "binomial", 
        control = glmer_control,
        nAGQ = 0,
        data = data)
  
}


# no MICE, default parameters
glmer_pooled_no_mice_default <- function(outcome, exposure_period, exposure, covars, data){
  
  glmer(formula(paste0(outcome, " ~ ", exposure_period, "_", exposure, " + factor(start_wk)", covars, "+ (1|CohortID)")),
        family = "binomial", 
        data = data)
  
}

# sensitivity analysis: cohort fixed effect
glm_pooled_fixed_effect <- function(outcome, exposure_period, exposure, covars, data){
  
  with(data, glm(formula(paste0(outcome, " ~ ", exposure_period, "_", exposure, " + factor(start_wk)", covars)), 
                 family = "binomial"))
  
}


# sensitivity analysis: random effects
glmer_pooled_random_effect <- function(outcome, exposure_period, exposure, covars, data){
  
  with(data, glmer(formula(paste0(outcome, " ~ ", exposure_period, "_", 
                                  exposure, " + factor(start_wk)", covars,
                                  "+ (", exposure_period, "_", exposure, "|CohortID)")),
                   family = "binomial", 
                   control = glmer_control,
                   nAGQ = 0)
  )
}


# PLR MICE model with default nAGQ

glmer_pooled_mice_defaultAGQ <- function(outcome, exposure_period, exposure, covars, data){
  
  with(data, glmer(formula(paste0(outcome, " ~ ", exposure_period, "_", exposure, " + factor(start_wk)", covars, "+ (1|CohortID)")),
                   family = "binomial", 
                   control = glmer_control)
  )
}


# stratified results - only show overall exposures (smokePM_mean and smoke_days)

plot_stratified_overall <- function(data, strat_variable) {
  
  y_axis_label <- case_when(data$outcome[1] %in% c("birth_ga_preterm", "preterm", 
                                                   "lbw", "term_lbw", "early_term") ~ "OR and 95% CI",
                            data$outcome[1] %in% c("birth_ga", "birth_bw", "birth_bwz_sex") ~ "Coefficient and 95% CI",
                            T ~ "Unknown")
  
  yint <- if_else(grepl("preterm|lbw|early_term", data$outcome[1]), 1, 0) 
  
  data$strat_var <- data[[strat_variable]]
  
  plot_overall <- ggplot(data, 
                         aes(x = modifier, y = result,  col = strat_var))+
    facet_wrap(~exposure, scales = "free") +
    geom_hline(yintercept = yint, col = "grey")+
    geom_point(position = position_dodge(width = 0.3))+
    geom_linerange(aes(ymin = lower, ymax = upper), 
                   position = position_dodge(width = 0.3))+
    theme_classic() +
    theme(
      # axis.text.x = element_blank(),
      # axis.ticks.x = element_blank(),
      text = element_text(size = 12),
      title = element_text(size = 11))+
    labs(x = "",
         col = data$modifier[1],
         y = y_axis_label)
  
  return(plot_overall)
  
}


# plot results from all trimesters

plot_res_trimesters <- function(data) {
  
  y_axis_label <- case_when(data$outcome[1] %in% c("birth_ga_preterm", "preterm", 
                                                   "lbw", "term_lbw", "early_term") ~ "OR and 95% CI",
                            data$outcome[1] %in% c("birth_ga", "birth_bw", "birth_bwz_sex") ~ "Coefficient and 95% CI",
                            T ~ "Unknown")
  
  yint <- if_else(grepl("preterm|lbw|early_term|sga|lga", data$outcome[1]), 1, 0) 
  
  plot_overall <- ggplot(dplyr::filter(data, exposure == "smokePM_mean"), 
                         aes(x = exposure, y = result))+
    facet_wrap(~exposure_period, nrow = 1)+
    geom_hline(yintercept = yint, col = "grey")+
    geom_point()+
    geom_linerange(aes(ymin = lower, ymax = upper))+
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          text = element_text(size = 12),
          title = element_text(size = 11))+
    # scale_y_continuous(transform = "log")+
    labs(title = expression(D. ~ Mean ~ smoke ~ PM[2.5]),
         x = "", y = y_axis_label)
  
  plot_days = ggplot(dplyr::filter(data, exposure_group == "Days"), 
                     aes(x = intensity, y = result))+
    facet_wrap(~exposure_period, nrow = 1)+
    geom_hline(yintercept = yint, col = "grey")+
    geom_point(position = position_dodge(width = 0.7))+
    geom_linerange(aes(ymin = lower, ymax = upper), 
                   position = position_dodge(width = 0.7))+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 12),
          title = element_text(size = 11), 
          axis.title.y = element_blank())+
    # scale_y_continuous(transform = "log", breaks = c(0.98, 1.02, 1.04 ,1.06, 1.08))+
    labs(x = expression("Concentration threshold ("  * mu * g ~ m^{-3} *")"), 
         title = "E. Smoke days")
  
  plot_waves = ggplot(filter(data, exposure_group == "Waves"), 
                      aes(x = intensity, col = duration, y = result))+
    facet_wrap(~exposure_period, nrow = 1)+
    geom_hline(yintercept = yint, col = "grey")+
    geom_point(position = position_dodge(width = 0.7), size = 1)+
    geom_linerange(aes(ymin = lower, ymax = upper), linewidth = 0.5, 
                   position = position_dodge(width = 0.7))+
    scale_color_viridis_d(option = "A", begin = 0.7, end = 0.2)+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 12),
          title = element_text(size = 11), 
          axis.title.y = element_blank())+
    # scale_y_continuous(transform = "log", breaks = c(0.5, 1, 1.5,  2))+
    labs(x = "",
         col = "Wave duration",
         title = "F. Smoke waves")
  
  plot_trimesters <- plot_overall + plot_days + plot_waves + plot_layout(ncol = 3, widths = c(2.5, 4, 4))
  plot_trimesters 
  
}
