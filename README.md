## Overview

This repo contains code supporting Sherris et al. (2025) Wildfire-specific fine particulate matter 
and preterm birth: a US ECHO Cohort analysis. The Lancet Planetary Health 101324. 

Scripts to assess exposure and define exposure metrics are in the folder "1_exposure assessment".
Analytic code to define models and generate results are in the folder "4_analysis".

For ECHO analysts, code download and process raw ECHO data and link with wildfire estimates
are available on the ECHO Analysis Workbench under the folder "EC0680/USERS/asherris/ANALYSIS/code".

Please direct queries to Allison Sherris, asherris@uw.edu or arsherris@gmail.com

### Exposure Assessment

This analysis uses daily estimates of census tract wildfire-specific PM2.5 from 
Child et al (2022) Daily local-level estimates of ambient wildfire smoke PM2.5 for the contiguous US.
Output from this analysis are available at https://github.com/echolab-stanford/daily-10km-smokePM.

This analysis also uses estimates of total ambient PM2.5 and ambient temperature described here:

Just AC, Arfer KB, Rush J, Lyapustin A, Kloog I. XIS-PM2.5: A daily spatiotemporal machine-learning model
for PM2.5 in the contiguous United States. Environmental Research. 2025 Apr 15;271:120948. 

Carrión D, Arfer KB, Rush J, Dorman M, Rowland ST, Kioumourtzoglou MA, Kloog I, Just AC. 
A 1-km hourly air-temperature model for 13 northeastern U.S. states using remotely sensed and 
ground-based measurements. Environ Res. 2021 Sept;200:111477. PMCID: PMC8403657


To apply the functions in in 1_exposure_assessment, raw data should be linked with
the above spatiotemporal models and processed in a data table with the following variables:

- xParticipant ID: Unique participant ID
- gest_day: Day of gestational (starting at conception = 0, ending at birth)
- smokePM_pred: Wildfire PM2.5 on the census tract of residence on the day of gestation
- total_pm25: Total ambient PM2.5 on the day of gestation
- mean_temp: Mean ambient temperature on the day of gestation

The output of these functions will be summary estimates for average smoke, smoke days,
and smoke waves across various exposure periods during gestation as follows:

*ec0680_1b_define_exposure_metrics.R*: Evaluate exposure metrics across fixed exposure periods:
pregnancy, 0-32 weeks, and by trimester.

*ec0680_1c_defined_exposure_metrics_survival.R*: Evaluate cumulative exposure metrics from
22 weeks through delivery for survival analysis.

*ec0680_1d_defined_ambientPM25_temp_metrics.R*: Evaluate summary metrics for ambient 
temperature and total PM2.5.

*ec0680_1e_define_weekly_exposure_metrics.R*: Evaluate weekly smoke days across pregnancy.


## Statistical analysis

The script *ec0680_4a_prep_for_analysis.R* defines functions and covariates to execute analytic code.