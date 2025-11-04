#### Create full test data frame including new ER data and rejected indivs from step screening ####

# Author: Read Barbee

# Date:2024-02-20
#Last updated: 2024-02-20

################################ Libraries #################################
library(tidyverse)
library(sf)
library(amt)
library(flextable)

################################ Functions #################################

#Create track for each animal
multi_track <- function(d){
  make_track(d, lon_utm, lat_utm, date_time_utc, check_duplicates = TRUE, all_cols = TRUE, crs = 5070)
}


###############################################################################

#1. Import data

###############################################################################

#earth ranger locations post May 2023
er_dat <- read_csv("data/gps_data/test_data/er_locs/test_locs_er_may23_feb24_dop_screened_02-08-24.csv") %>%
  rename(animal_id = subject,
         date_time_utc = timestamp)

#collar deployments
md <- read_csv("data/gps_data/metadata/deployments/collar_deployments_master_2-15-2024.csv")

#dispersals
disp <- read_csv("data/gps_data/metadata/dispersals/nsd/dispersals_master_reduced_02-15-2024.csv")

#locations rejected from step building
rejects <- read_csv("data/gps_data/test_data/step_rejects_pre_May2023/test_dat_step_reject_indivs_pre2024_2-15-24.csv")


#extract tracking dates from collar deployments
tracking_dates <- md %>% 
  mutate(start_date = mdy(start_date),
         end_date = mdy(end_date)) %>% 
  group_by(name) %>% 
  summarize(start_date = min(start_date, na.rm=T),
            end_date = max(end_date))


########################################################################
##
## 2. Earth Ranger data initial formatting
##
##########################################################################

#add coordinates in crs 5070
er_dat <- er_dat %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = 5070) %>%
  mutate(lon_utm = st_coordinates(.)[,1],
         lat_utm = st_coordinates(.)[,2],
         .after = lon) %>% 
  as_tibble() %>% 
  select(-geometry)

########################################################################
##
## 3. Trim EarthRanger locations by tracking dates
##
##########################################################################

er_nested <- er_dat %>% 
  nest_by(animal_id)

#make individual tracks
er_nested$tracks <-map(er_nested$data, multi_track)

#trim tracks by tracking dates

trimmed_tracks <- list()

for(i in 1:nrow(er_nested)){
  name_i <- er_nested$animal_id[i]
  
  start_date <- tracking_dates %>% 
    filter(name == name_i) %>% pull(start_date)
  
  end_date <- tracking_dates %>% 
    filter(name == name_i) %>% pull(end_date)
  
  if(is.na(end_date)==T){
    trimmed_tracks[[i]] <- er_nested$tracks[[i]] %>% tracked_from_to(from = start_date)
    
  } else if(is.na(end_date)==F){
    trimmed_tracks[[i]] <- er_nested$tracks[[i]] %>% tracked_from_to(from = start_date, to = end_date)
  }
  
  print(paste0(i, "/", nrow(er_nested)))
}

er_nested$tracks <- trimmed_tracks

########################################################################
##
## 4. Add metadata to ER tracks
##
##########################################################################
disp_names <- disp %>% distinct(animal_id) %>% pull(animal_id)
disp_dates <- disp %>% select(animal_id, disp_date_nsd)

sexes <- md %>% 
  distinct(name, .keep_all = T) %>% 
  select(name, sex) %>% 
  rename(animal_id = name)

er_dem_dat <- er_nested %>% 
  left_join(sexes, by="animal_id") %>% 
  left_join(disp_dates, by = "animal_id") %>% 
  mutate(dispersal_status = case_when(animal_id %in% disp_names ~ "disperser",
                                      .default = "resident" )) %>% 
  select(animal_id, sex, dispersal_status, disp_date_nsd, data) %>% 
  mutate(disp_date_nsd = mdy(disp_date_nsd))

########################################################################
##
## 5. Remove locations from before dispersal events 
##
##########################################################################

new_dat <- list()
for(i in 1:nrow(er_dem_dat)){
  disp_date <- er_dem_dat$disp_date_nsd[i]
  
  if(is.na(disp_date) == F){
    
    new_dat[[i]] <- er_dem_dat$data[[i]] %>% filter(date_time_utc >= disp_date)
  } else {
    new_dat[[i]] <- er_dem_dat$data[[i]]
  }
  
  print(paste0(i, "/", nrow(er_dem_dat)))
}


er_dem_dat$data <- new_dat


#clean up
er_dem_dat_final <- er_dem_dat %>% 
  unnest(cols=data) %>% 
  ungroup() %>% 
  select(-disp_date_nsd) %>% 
  rename(utm_e = lon_utm,
         utm_n = lat_utm) %>% 
  relocate(lon, .before = lat)

#check for duplicates 
janitor::get_dupes(er_dem_dat_final, animal_id, date_time_utc, lon, lat)

########################################################################
##
## 6. Add locations from step screening rejects
##
##########################################################################

rejects_sel <- rejects %>% 
  select(animal_id, 
         sex, 
         dispersal_status,
         date_time_utc,
         lon_wgs84,
         lat_wgs84,
         lon_utm,
         lat_utm,
         dop) %>% 
  rename(lat = lat_wgs84,
         lon = lon_wgs84,
         utm_e = lon_utm,
         utm_n = lat_utm)

#check for duplicates
janitor::get_dupes(rejects_sel, animal_id, date_time_utc, lon, lat)


#combine er locations and step rejects
test_dat <- bind_rows(er_dem_dat_final, rejects_sel)

#check for duplicates
janitor::get_dupes(test_dat, animal_id, date_time_utc, lon, lat)


#write_csv(test_dat, "test_dat_full_er_and_amt_rejects_2-20-24.csv")

########################################################################
##
## 7. Make summary table for test data
##
##########################################################################

loc_counts <- test_dat %>% 
  group_by(sex, dispersal_status) %>% 
  count()

indiv_counts <- test_dat %>% 
  distinct(animal_id, .keep_all = T) %>% 
  group_by(sex, dispersal_status) %>% 
  count()

summ_tab <- indiv_counts %>% 
  left_join(loc_counts, by=c("sex", "dispersal_status")) %>% 
  rename(individuals = n.x,
         locations = n.y)


ftab <- flextable(summ_tab)

#save_as_docx("table 1" = ftab, path = "test_data_summary_tab_2-20-24.docx")
