## EC0680
## Exposure assessment: Assign exposure metrics
## Start date: 6/7/23
## A. Sherris

# load data
load("processed data/1_exposure_assessment/wf_pm25_exposure_deidentified.RData")

wf_pm25_exposure <- wf_pm25_exposure_deid %>% 
  mutate(final_gest_day = (birth_ga*7),
         gest_day = gest_day+1)  # fix mistake in 1a_link_exposuer_data version sent to DAC
  

# write function to assign wf exposure

wf_exposure_weekly <- wf_pm25_exposure %>% 
  filter(gest_day >= 0) %>%
  arrange(gest_day) %>% 
  mutate(gest_week = floor(gest_day/7),
         PM_day = if_else(smokePM_pred > 0, 1, 0),
         PM_over25 = if_else(smokePM_pred >= 2.5, 1, 0),
         PM_over5  = if_else(smokePM_pred >= 5, 1, 0),
         PM_over10  = if_else(smokePM_pred >= 10, 1, 0)
         ) %>% 
   # summarize exposure metrics by participant and week
  group_by(xParticipantID, gest_week) %>% 
  summarize(weekly_smoke_days = sum(PM_day),
            weekly_smoke_days_over25 = sum(PM_over25),
            weekly_smoke_days_over5   = sum(PM_over5),
            weekly_smoke_days_over10  = sum(PM_over10)
            ) %>% 
    ungroup()

# save output
save(wf_exposure_weekly,
     file = "processed data/1_exposure_assessment/wf_pm25_exposure_weekly.RData")

# end
