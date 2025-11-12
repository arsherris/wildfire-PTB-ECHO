## EC0680
## Exposure assessment: Assign exposure metrics
## Start date: 6/7/23
## A. Sherris

# load data
load("processed data/1_exposure_assessment/wf_pm25_exposure_deidentified.RData")

# write function to assign wildfire exposure based on daily wildfire-specific PM2.5 
# estimates during pregnancy from Childs et al. 2022

summarize_wf_exposure <- function(data, prefix) {
  
 data %>% 
    arrange(gest_day) %>% 
    group_by(xParticipantID) %>% 
    # calculate exposure metrics by day
    mutate(PM_over25 = if_else(smokePM_pred >= 2.5, 1, 0),
           PM_over5  = if_else(smokePM_pred >= 5, 1, 0),
           PM_over10  = if_else(smokePM_pred >= 10, 1, 0),
            # weekly rolling average of smoke PM2.5
           PM_week_mean = rollmean(smokePM_pred, k = 7, fill = NA, align = "right"),
           # smoke waves
           wave2_over25 = if_else(lag(PM_over25, 1) == 1 & PM_over25 == 1, 1, 0),
           wave3_over25 = if_else(lag(PM_over25, 2) == 1 & 
                                    lag(PM_over25, 1) == 1 & PM_over25 == 1, 1, 0),
           wave4_over25 = if_else(lag(PM_over25, 3) == 1 & 
                                    lag(PM_over25, 2) == 1 & 
                                    lag(PM_over25, 1) == 1 & PM_over25 == 1, 1, 0),
           wave2_over5  = if_else(lag(PM_over5, 1) == 1 & PM_over5 == 1, 1, 0),
           wave3_over5  = if_else(lag(PM_over5, 2) == 1 & 
                                    lag(PM_over5, 1) == 1 & PM_over5 == 1, 1, 0),
           wave4_over5  = if_else(lag(PM_over5, 3) == 1 & 
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
    # summarize exposure metrics by participant
    summarize(gest_day_count = n(),
              smokePM_mean = mean(smokePM_pred),
              smokePM_weekly_max = max(PM_week_mean, na.rm=T),
              smoke_days = sum(smokePM_pred > 0),
              smoke_days_over25 = sum(PM_over25),
              smoke_days_over5   = sum(PM_over5),
              smoke_days_over10  = sum(PM_over10),
              smoke_wave2_over25 = sum(wave2_over25, na.rm = T),
              smoke_wave3_over25 = sum(wave3_over25, na.rm = T),
              smoke_wave4_over25 = sum(wave4_over25, na.rm = T),
              smoke_wave2_over5 = sum(wave2_over5, na.rm = T),
              smoke_wave3_over5 = sum(wave3_over5, na.rm = T),
              smoke_wave4_over5 = sum(wave4_over5, na.rm = T),
              smoke_wave2_over10 = sum(wave2_over10, na.rm = T),
              smoke_wave3_over10 = sum(wave3_over10, na.rm = T),
              smoke_wave4_over10 = sum(wave4_over10, na.rm = T)) %>% 
    ungroup %>% 
    rename_at(vars(-xParticipantID), ~paste0(prefix, "_", .)) %>% 
    return
  
}

# define exposure during pregnancy 

# exposure period doc to dob
wf_pm25_pregnancy <- wf_pm25_exposure %>%
  filter(gest_day >= 0,
         gest_day <= final_gest_day) %>% 
  summarize_wf_exposure(., "preg") 

# exposure period doc to term
wf_pm25_term <- wf_pm25_exposure %>%
  filter(birth_ga >= 37,
         gest_day >= 0,
         gest_day <= 259) %>% 
  summarize_wf_exposure(., "term") 

# exposure period doc to 32 wk
wf_pm25_32wk <- wf_pm25_exposure %>%
  filter(birth_ga >= 32,
         gest_day >= 0,
         gest_day <= 224) %>% 
  summarize_wf_exposure(., "wk32") 


# define exposure during trimesters

wf_pm25_trimesters <- wf_pm25_exposure %>% 
  # remove exposure before conception
  filter(gest_day >= 0) %>%
  group_by(xParticipantID) %>% 
  # Note - check start/end dates for trimesters (eg, T1 goes to end of 13 weeks)
  mutate(T1 = if_else(gest_day %in% 0:97, 1, 0), #0-13 weeks after doc
         T2 = if_else(gest_day %in% 98:188, 1, 0), # 14-26 wks
         T3 = case_when(birth_ga < 31 ~ NA_integer_,
                        gest_day %in% (first(final_gest_day) - 27):first(final_gest_day) ~ 1,
                        T ~ 0)
         ) %>% 
  ungroup 


wf_pm25_T1 <- wf_pm25_trimesters %>%
  filter(T1 == 1) %>%
  summarize_wf_exposure(., "T1")

wf_pm25_T2 <- wf_pm25_trimesters %>%
  filter(T2 == 1) %>%
  summarize_wf_exposure(., "T2")

wf_pm25_T3 <- wf_pm25_trimesters %>%
  filter(T3 == 1) %>%
  summarize_wf_exposure(., "T3")

# combine exposure estimates
wf_pm25_exposure_summary <- wf_pm25_pregnancy %>% 
  left_join(wf_pm25_T1) %>% 
  left_join(wf_pm25_T2) %>% 
  left_join(wf_pm25_T3) %>% 
  left_join(wf_pm25_term) %>% 
  left_join(wf_pm25_32wk)

# save output
save(wf_pm25_exposure_summary,
     file = "processed data/1_exposure_assessment/wf_pm25_exposure_summary_deid.RData")


# end 