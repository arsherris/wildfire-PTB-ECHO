## EC0680
## Project setup
## Start date: 6/7/23
## A. Sherris

# set working directory

setwd("\\\\echofile.rti.ns/EC0680/USERS/asherris/ANALYSIS_PTB")

# load packages

library(tidyverse)
library(lubridate)
library(zoo)
library(mice)

# load exposure list

exposure_list <-  c("smokePM_mean",
                    "smoke_days",
                    "smoke_days_over25", "smoke_days_over5", "smoke_days_over10",
                    "smoke_wave2_over25", "smoke_wave3_over25", "smoke_wave4_over25",
                    "smoke_wave2_over5", "smoke_wave3_over5", "smoke_wave4_over5",
                    "smoke_wave2_over10", "smoke_wave3_over10", "smoke_wave4_over10")

