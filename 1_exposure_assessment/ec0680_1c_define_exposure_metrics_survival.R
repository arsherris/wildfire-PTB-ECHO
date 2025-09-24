## EC0680
## Exposure assessment: Assign exposure metrics for longitudinal survival analysis
## Start date: 8/7/23
## A. Sherris

load("processed data/1_exposure_assessment/wf_pm25_exposure_deidentified.RData")

# function to evaluate cumulative WFS exposure to the start of each gestational week
# for use in long form data / survival analysis

summarize_cumulative_wf_exposure <- function(data, prefix) {
  
  data %>% 
    # calculate exposure metrics by day
    mutate(smoke_day = if_else(smokePM_pred > 0, 1, 0),
           PM_over25 = if_else(smokePM_pred >= 2.5, 1, 0),
           PM_over5  = if_else(smokePM_pred >= 5, 1, 0),
           PM_over10  = if_else(smokePM_pred >= 10, 1, 0)) %>% 
    arrange(xParticipantID, gest_day) %>% 
    group_by(xParticipantID) %>% 
    mutate(gest_week = as.integer(gest_day/7),
           # cumulative cumulative smoke mean and days
           smoke_cummean = cummean(smokePM_pred),
           smoke_cumdays = cumsum(smoke_day),
           cumdays_over25 = cumsum(PM_over25),
           cumdays_over5 = cumsum(PM_over5),
           cumdays_over10 = cumsum(PM_over10),
           # identify smoke waves
           wave2_over25 = if_else(lag(PM_over25, 1) == 1 & PM_over25 == 1, 1, 0),
           wave3_over25 = if_else(lag(PM_over25, 2) == 1 &
                                    lag(PM_over25, 1) == 1 & PM_over25 == 1, 1, 0),
           wave4_over25 = if_else(lag(PM_over25, 3) == 1 &
                                    lag(PM_over25, 2) == 1 &
                                    lag(PM_over25, 1) == 1 & PM_over25 == 1, 1, 0),
           wave2_over5 = if_else(lag(PM_over5, 1) == 1 & PM_over5 == 1, 1, 0),
           wave3_over5 = if_else(lag(PM_over5, 2) == 1 &
                                   lag(PM_over5, 1) == 1 & PM_over5 == 1, 1, 0),
           wave4_over5 = if_else(lag(PM_over5, 3) == 1 &
                                   lag(PM_over5, 2) == 1 &
                                   lag(PM_over5, 1) == 1 & PM_over5 == 1, 1, 0),
           wave2_over10 = if_else(lag(PM_over10, 1) == 1 & PM_over10 == 1, 1, 0),
           wave3_over10 = if_else(lag(PM_over10, 2) == 1 &
                                    lag(PM_over10, 1) == 1 & PM_over10 == 1, 1, 0),
           wave4_over10 = if_else(lag(PM_over10, 3) == 1 &
                                    lag(PM_over10, 2) == 1 &
                                    lag(PM_over10, 1) == 1 & PM_over10 == 1, 1, 0)) %>% 
    # count each wave as one wave (regardless of length of wave)
    mutate_at(vars(starts_with("wave")), ~if_else(lag(.) == 1, 0, .)) %>% 
    # calculate cumulative sum of waves
    mutate_at(vars(starts_with("wave")), ~cumsum(replace_na(., 0))) %>% 
    group_by(xParticipantID, gest_week) %>% 
    # summarize variables for each gestational week
    summarize(birth_ga = first(birth_ga),
              start_day = first(gest_day),
              end_day = last(gest_day),
              weeklymean_smokePM = mean(smokePM_pred),
              smokePM_mean = mean(smoke_cummean),
              smoke_days = max(smoke_cumdays),
              smoke_days_over25 = max(cumdays_over25),
              smoke_days_over5 = max(cumdays_over5),
              smoke_days_over10 = max(cumdays_over10),
              smoke_wave2_over25 = max(wave2_over25),
              smoke_wave3_over25 = max(wave3_over25),
              smoke_wave4_over25 = max(wave4_over25),
              smoke_wave2_over5 = max(wave2_over5),
              smoke_wave3_over5 = max(wave3_over5),
              smoke_wave4_over5 = max(wave4_over5),
              smoke_wave2_over10 = max(wave2_over10),
              smoke_wave3_over10 = max(wave3_over10),
              smoke_wave4_over10 = max(wave4_over10)) %>% 
    ungroup() %>% 
    mutate(start_wk = gest_week + 1,  # start of gestational week - exposure reflects previous completed week of exposure
           end_wk = gest_week + 2) %>% 
    filter(gest_week<37) %>% 
    select(xParticipantID, start_wk, end_wk, start_day, end_day, birth_ga, starts_with("smoke")) %>% 
      rename_at(vars(starts_with("smoke")), ~paste0(prefix, "_", .)) %>% 
      return()
  
}


# entire pregnancy cumulative exposure
wf_exposure_cumulative <- wf_pm25_exposure_deid %>% 
  mutate(final_gest_day = birth_ga*7,
       gest_day = gest_day+1) %>%  # fix mistake in 1a_exposure_assessment_prep
  filter(gest_day >= 0, 
         gest_day <= final_gest_day) %>%
  # calculate exposure metrics by day
  summarize_cumulative_wf_exposure(., "surv") 


save(wf_exposure_cumulative, 
     file = "processed data/1_exposure_assessment/wf_pm25_exposure_cumulative_deid.RData")

