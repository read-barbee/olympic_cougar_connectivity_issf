#### GPS Data Screening ####

# Author: Read Barbee

# Date:2023-06-01
# Last updated:2024-02-08

# Purpose: Screen GPS data for erroneous locations and capture effects

# Inputs:
#   1. Raw GPS locations

# Outputs:
#   1. Screened GPS locations

# Steps:
# 1. Calculate global fix rate 
# 2. Check for missing data
# 3. Remove first 24hr of fixes for each animal to avoid capture effects
# 4. Remove 2D fixes with dop score > 5
# 5. Remove locations outside extent of study area raster 


################################ Libraries #################################
library(tidyverse)
library(terra)
library(sf)
library(lubridate)
library(amt)
library(janitor)
library(flextable)
library(DataExplorer)
library(mapview)

################################ Functions #################################
#calculate fix rate
calc_fix_rate <- function(data){
  
  success <- data %>% 
    filter(is.na(lat_wgs84)==FALSE) %>% 
    count() %>% pull()
  
  total = nrow(data)
  
  fix_rate = success/total
  
  return(fix_rate)
  
}

################################ User-Defined Parameters #################################

#buffer distance around study area in meters
study_area_buffer_dist <- 1000  #500m allows for walking out on beaches at low tide

#########################################################################
##
## 1. Import and format data
##
##########################################################################

#raw locations (April 2023)
locs_raw <- read_csv("data/gps_data/master_locations/raw/gps_locs_master_10-02-2023.csv") #,  col_types = list(fix_type = col_character())

#collar deployments (April 2023)
deployments <- read_csv("data/gps_data/metadata/deployments/collar_deployments_master_7-11-2023.csv")

# Dispersal inventory from Teams (3-13-2023)
dispersals <- read_csv("data/gps_data/metadata/dispersals/project_records/dispersals_master_10-02-2023.csv")

#Import dispersal dates (from nsd analysis)
dispersal_dates <- read_csv("data/gps_data/metadata/dispersals/nsd/dispersals_master_reduced_10-02-2023.csv") %>% 
 mutate(disp_date_nsd = mdy(disp_date_nsd)) %>%
  select(-sex)

#check for duplicates in the dispersal dates file
get_dupes(dispersal_dates, animal_id)

#import study area boundary
op_poly <- st_read("data/spatial/ocp_study_area_poly_wa_only_10-30-23.gpkg")

#buffer study area by specified distance 
op_poly_buffer <- st_buffer(op_poly, study_area_buffer_dist)


#########################################################################
##
## 2. Data Diagnostics
##
##########################################################################

#calculate global gps fix rate
fix_success_rate <- calc_fix_rate(locs_raw)

#calculate gps fix rate for each individual 
indiv_fix_rate <- locs_raw %>% 
  nest_by(animal_id) %>% 
  summarize(fix_success = calc_fix_rate(data))

#plot distribution of individual fix rates
hist(indiv_fix_rate$fix_success)

#Make sure no essential fields are missing data
summary(locs_raw %>% filter(!is.na(lat_wgs84)))
plot_missing(locs_raw)

#make sure no locations with coordinates are missing dop scores. Get the names of the individuals that are if any
locs_raw %>% 
  filter(!is.na(lat_wgs84) & is.na(dop)) %>% 
  distinct(animal_id)


#########################################################################
##
## 3. Add demographic metadata
##
##########################################################################

#Filter deploment list to get one row per individual
first_deps <- deployments %>% 
  clean_names() %>% 
  distinct(name, .keep_all = TRUE) %>% 
  mutate(name = trimws(name))

#Get vector of disperser names
disperser_names <- dispersals %>% 
  clean_names() %>% 
  filter(!is.na(animal_id)) %>% 
  distinct(animal_id) %>% 
  pull(animal_id)


#Add dispersal status to deployment list
dem_cats <- first_deps %>% 
  mutate(dispersal_status = case_when(name %in% disperser_names == TRUE ~ "disperser",
                                      name %in% disperser_names == FALSE ~ "resident")) %>% 
  select(name,
         sex,
         dispersal_status) %>% 
  rename(animal_id=name)


#Nest locations by animal_id and add columns for sex and dispersal status
locs_nested <- locs_raw %>% 
  nest_by(animal_id) %>% 
  left_join(dem_cats, by = join_by(animal_id)) %>% 
  left_join(.,dispersal_dates, by = join_by(animal_id)) %>% 
  select(animal_id, sex, dispersal_status:disp_qual, data)


locs_dem <- locs_nested %>% unnest(cols = data) %>% ungroup()

get_dupes(locs_dem, deployment_id, date_time_utc)

#########################################################################
##
## 4. Create column indicating which steps were during active dispersals
##
##########################################################################

locs_dem <- locs_dem %>% 
  mutate(dispersing = case_when(date_time_local >= disp_date_nsd ~ TRUE,
                                date_time_local < disp_date_nsd ~ FALSE,
                                is.na(disp_date_nsd) == TRUE ~ NA), .after = disp_qual)



#########################################################################
##
## 5. Remove missing fixes and locations with missing DOP scores
##
##########################################################################

#remove missing locations
locs_raw_no_na <- locs_dem %>% filter(!is.na(lat_wgs84))

#remove locations that are missing dop scores
locs_raw_filt <- locs_raw_no_na %>% 
  filter(!is.na(lat_wgs84) & !is.na(dop))

#########################################################################
##
## 6. Remove capture effects
##
##########################################################################


# Remove first 24 hours for each animal to avoid capture effects (removes 1,781 points). Remove 2D fixes with dop > 5 (removes 27,397 points)

cap_eff <- locs_raw_filt %>%  ## Removes Belle who only had 16 locations #
  group_by(animal_id) %>% 
  filter(date_time_local >= (min(date_time_local) + hours(24))) %>% 
  ungroup() %>% 
  filter(!is.na(lat_wgs84))


#compare to minimum date times from original data frame to make sure it worked
locs_raw %>%
  group_by(animal_id) %>%
  summarize(min=min(date_time_local))

cap_eff %>%
  group_by(animal_id) %>%
  summarize(min=min(date_time_local))

#########################################################################
##
## 7. Remove low quality locations (DOP filter)
##
##########################################################################  
  
dop_filt <- cap_eff %>%  filter(dop<=5 | fix_type == "3D") %>%
  mutate(unique_id= 1:nrow(.), .before = deployment_id) 


#########################################################################
##
## 9. Remove locations outside of study area 
##
##########################################################################

#remove_locations outside of study area if not already done in screening stage (should remove 2 pts)
geo_filt <- dop_filt %>% 
  st_as_sf(coords = c("lon_utm", "lat_utm"), crs = 5070, remove=FALSE) %>% 
  mutate(intersects_study_area = lengths(st_intersects(., op_poly_buffer)) > 0) %>% 
  filter(intersects_study_area==TRUE)

#convert geo-sreened sf object back to dataframe for comparison with df before geo-screenieng
geo_screened_df <- geo_filt %>% 
  as.data.frame() %>% 
  select(unique_id, deployment_id, animal_id, collar_id, everything(), -c(geometry, intersects_study_area))

#get dataframe of differences before and after geo-screening
geo_removed <- setdiff(dop_filt, geo_screened_df)


#visualize in mapview
geo_removed %>% sf::st_as_sf(coords=c("lon_wgs84", "lat_wgs84"), crs=4326, remove=FALSE) %>% 
  mapview()


#########################################################################
##
## 10. Remove locations for dispersers not during active dispersal events
##
##########################################################################

disp_screened <- geo_screened_df %>%
  filter(dispersal_status=="resident" | dispersing ==TRUE)


#########################################################################
##
## 11. Create summary table
##
##########################################################################

#tally locations removed at each step
location_attempts <- nrow(locs_raw)
successful_locations <- nrow(locs_raw_no_na)
missing_dop_diff <- nrow(locs_raw_no_na) - nrow(locs_raw_filt)
capture_effects_diff <- nrow(locs_raw_filt) - nrow(cap_eff)
dop_filter_diff <- nrow(cap_eff) - nrow(dop_filt)
geo_filter_diff <- nrow(dop_filt) - nrow(geo_screened_df)
disp_filter_diff <- nrow(geo_screened_df) - nrow(disp_screened)

remaining <- nrow(disp_screened)

#make table
summary_table_diff <- tibble(`Screening step` = c("Location attempts",
                                 "Successful locations", 
                                 "DOP filter", 
                                 "Study area filter", 
                                 "Capture effects (24 hr)", 
                                 "Post-dispersal filter",
                                 "Total locations for analysis"),
                        Locations = c(location_attempts,
                                              successful_locations,
                                              -dop_filter_diff,
                                              -geo_filter_diff,
                                              -capture_effects_diff,
                                              -disp_filter_diff,
                                              remaining)
                        )

summary_table_diff

#make flex table
diff_table <- flextable(summary_table_diff) %>% 
  colformat_double(digits = 0) %>% 
  colformat_double(i = ~ `Screening Step` == "Success rate", digits = 2) %>% 
  width(width= 2) %>% 
  add_footer_lines(paste0("GPS fix success rate: ", round(fix_success_rate, 2)))
 
#export flextable
save_as_image(diff_table, "gps_screening_diff_table_training_dat_2-13-24.png")
save_as_docx(diff_table, path = "gps_screening_diff_table_training_dat_2-13-24.docx")

#check for duplicate locations
get_dupes(disp_screened, deployment_id, date_time_utc)

###############################################################################

#12. Export screened locations

###############################################################################
#write_csv(disp_screened, "data/gps_data/master_locations/screenedgps_locs_dop_disp_screened_02-13-2024.csv")

