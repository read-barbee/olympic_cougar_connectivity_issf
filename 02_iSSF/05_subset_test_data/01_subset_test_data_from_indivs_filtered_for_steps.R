#### Create Dataframe of test data from indivs screened from step building process ####

# Author: Read Barbee

# Date:2024-02-15
#Last updated: 2024-02-20

################################ Libraries #################################
library(tidyverse)
library(amt)
library(sf)
library(janitor)

################################ Functions #################################

#Create track for each animal
multi_track <- function(d){
  make_track(d, lon_utm, lat_utm, date_time_utc, check_duplicates = TRUE, all_cols = TRUE, crs = 5070)
}

########################################################################
##
## 1. Subset data for indivs not included in step analysis and export
##
##########################################################################

#Training data before converting to steps
og_locs <- locs_screened <- read_csv("data/gps_data/master_locations/screened/gps_locs_dop_disp_screened_02-13-2024.csv")

#Training data after converting to steps (fewer indivs)
filtered_steps <- read_csv("data/gps_data/steps/2h_steps_unscaled_no_imp_annual_cov_02-13-2024.csv")

#extract animal_ids from full data
og_ids <- og_locs %>% 
  distinct(animal_id) %>% 
  pull(animal_id)

#extract animal_ids from step data
filtered_ids <- filtered_steps %>% 
  distinct(animal_id) %>% 
  pull(animal_id)

#get list of individuals in the og data that aren't in the step data
filtered_indivs <- setdiff(og_ids, filtered_ids)


#subset locations from the original data belonging to those individuals
filtered_indiv_locs <- og_locs %>% 
  filter(animal_id %in% filtered_indivs)


# write_csv(filtered_indiv_locs, "test_dat_step_reject_indivs_pre2024_2-15-24.csv")

