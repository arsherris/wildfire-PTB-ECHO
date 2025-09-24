## EC0680
## Exposure assessment: Evaluate total PM2.5 and temperature during pregnancy
## Start date: 6/7/23
## A. Sherris

# load data
load("processed data/1_exposure_assessment/wf_exposure_deidentified.RData")

wf_pm25_exposure <- wf_pm25_exposure_deid %>% 
  mutate(final_gest_day = (birth_ga*7),
         gest_day = gest_day+1)  # fix mistake in 1a_link_exposuer_data version sent to DAC

# write function to assign temperature exposure

summarize_totalPM_temp <- function(data, prefix) {
  
 data %>% 
    arrange(gest_day) %>% 
    group_by(xParticipantID) %>% 
    # summarize exposure metrics by participant
    summarize(percNA_totalPM = sum(is.na(total_pm25)) / n(),  #  missingness is the same for totalPM, mean temp, max temp
              totalPM_mean = mean(total_pm25, na.rm=T),
              mean_temp = mean(mean_temp, na.rm=T),
              mean_maxtemp = mean(max_temp, na.rm=T)
              ) %>%
    mutate_at(vars(totalPM_mean:mean_maxtemp), ~ if_else(percNA_totalPM > 0.5, NA_real_, .)) %>% 
    rename_at(vars(-xParticipantID), ~paste0(prefix, "_", .)) %>% 
    return
  
}


# define exposure during pregnancy 

# exposure period doc to dob
totalPM_pregnancy <- wf_pm25_exposure %>%
  filter(gest_day >= 0,
         gest_day <= final_gest_day) %>% 
  summarize_totalPM_temp(., "preg") 

save(totalPM_pregnancy, file = "processed data/1_exposure_assessment/totalPM_pregnancy.RData")

# end
