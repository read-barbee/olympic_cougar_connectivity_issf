#### Combine all location source files ####

# Author: Read Barbee

# Date:2023-05-16 

# Last updated:2023-09-26 

# Purpose: Combine all location source files:
#1. web source
#2. hist_source
#3. collar_source

#total deployments: 171


###############################################################################

#### Library ####
library(tidyverse)
library(lubridate)
library(janitor)
library(sf)


#########################################################################
##
## 1. Import and format source files
##
##########################################################################
web_source <- read_csv("data/Location_Data/Source_Files/sources/web_source_5-16-2023.csv", col_types = list(fix_type = col_character(),
                                                                                                    collar_id = col_character())) %>% 
  mutate(source = "web")

hist_source <- read_csv("data/Location_Data/Source_Files/sources/hist_source_5-16-2023.csv",col_types = list(fix_type = col_character(),
                                                                                                     collar_id = col_character())) %>% 
  mutate(source = "historic")

collar_source <- read_csv("data/Location_Data/Source_Files/sources/collar_source_5-16-2023.csv",col_types = list(fix_type = col_character(),
                                                                                                         collar_id = col_character())) %>% 
  mutate(source = "collar_download")

deployments_master <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/Deployments/collar_deployments_master_7-11-2023.csv") %>%
  mutate(deployment_id = paste0(name,"_",collar_id), .before=name,
         start_date = mdy(start_date),
         end_date = mdy(end_date)) 
#%>% distinct(deployment_id) %>% pull()

collar_meta <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/Collars/collar_model_list_partial.csv") %>% 
  janitor::clean_names() %>% 
  select(collar_id, collar_brand, collar_type, satellite_system) %>% 
  mutate(collar_id = as.character(collar_id))


#########################################################################
##
## 2. Combine source files into single dataframe
##
##########################################################################
locs_all <- bind_rows(web_source,
                      hist_source,
                      collar_source)


#########################################################################
##
## 3. Check for duplicate points and completeness
##
##########################################################################

#check for duplicate points
get_dupes(locs_all, deployment_id, date_time_utc, latitude, longitude)
get_dupes(locs_all, deployment_id, date_time_utc)


#make sure the location file contains all the deployments in the master deployment list
locs_all_deps <- locs_all %>% distinct(deployment_id) %>% pull()

setdiff(deployments_master %>% pull(deployment_id), locs_all_deps)
setdiff(locs_all_deps, deployments_master %>% pull(deployment_id))

#create local time column
locs_all <- locs_all %>% 
  mutate(date_time_local = with_tz(date_time_utc, tzone="US/Pacific"))

#check to make sure all necessary fields are complete
summary(locs_all)

#21 locations have coordinates but NA for dop because of excel find and replace. All of those NA should be 0 to keep the column numeric for filtering
locs_all %>% filter(!is.na(latitude) & is.na(dop))

#replace all NAs in dop column with 0
locs_all <- locs_all %>% mutate(dop=replace_na(dop, 0))


#########################################################################
##
## 4. Trim locations for each individual by deployment dates
##
##########################################################################

locs_all_nested <- locs_all %>% nest_by(deployment_id)

#Create function to extract deployments based on collar ID and deployment dates and append columns of animal_id and deployment_id 

extract_deployments <- function(data_p, animal_id_p, collar_id_p, start_date_p, end_date_p) {
  if (is.na(end_date_p) == TRUE){ #don't filter by end date if it's not included
    trimmed_track <- data_p %>% 
      filter(collar_id == collar_id_p) %>% 
      filter(date_time_local>= ymd(start_date_p)) %>% 
      mutate(animal_id = animal_id_p,
             deployment_id = paste0(animal_id_p,"_",collar_id_p)) %>% 
      select(deployment_id, animal_id, everything())
  } else{
    trimmed_track <- data_p %>% 
      filter(collar_id == collar_id_p) %>% 
      filter(date_time_local>= ymd(start_date_p) & date_time_local<=ymd(end_date_p)) %>% 
      mutate(animal_id = animal_id_p,
             deployment_id = paste0(animal_id_p,"_",collar_id_p)) %>% 
      select(deployment_id, animal_id, everything())
  }
  
  return(trimmed_track)
}
#create empty list to fill with data frames of trimmed individual tracks
trimmed_deployments <- list()

#loop to generate trimmed tracks for each deployment and append them to list (working)
for (i in 1:(nrow(locs_all_nested)+1)) {
  trimmed_deployments[[i]]<- extract_deployments(
    data_p = locs_all,
    animal_id_p = deployments_master$name[[i]],
    collar_id_p = deployments_master$collar_id[[i]],
    start_date_p = deployments_master$start_date[[i]],
    end_date_p = deployments_master$end_date[[i]])
}

#combine all trimmed deployments into single dataframe
trimmed_all <- bind_rows(trimmed_deployments)

trimmed_all <- trimmed_all %>% 
  filter(!(animal_id=="Moxie" & date_time_utc >ymd("2022-04-17"))) %>% #remove Moxie's post-relocation points
  filter(deployment_id!="Gypsy_87521") %>%  #remove Gypsy's post-relocation points
  filter(!(animal_id=="Junior" & date_time_utc <ymd("2014-02-26")))

trimmed_sf <- trimmed_all %>% sf::st_as_sf(coords =c("longitude", "latitude"), crs = 4326, na.fail=FALSE) %>% nest_by(animal_id)

#Inspect tracks to make sure they're trimmed correctly

#inspect by name
inspect <- trimmed_sf %>% filter(animal_id == "Junior")
mapview::mapview(inspect$data) 

#inspect by number
mapview::mapview(trimmed_sf$data[[1]]) #needs geoscreening: 16, 29, 47, 63, 93, 108, 110


#with_tz doesn't print to csv unless coerced to a character.
trimmed_all <- trimmed_all %>%
  mutate(date_time_utc = as.character(date_time_utc),
         date_time_local = as.character(date_time_local))

utm_coords <- trimmed_all %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326, na.fail=FALSE) %>% 
  st_transform(crs=5070) %>% 
  st_coordinates()




#########################################################################
##
## 5. Create final output
##
##########################################################################

final <- trimmed_all %>% 
  rename(lat_wgs84 = latitude,
         lon_wgs84 = longitude) %>% 
  mutate(lon_utm = utm_coords[,1],
         lat_utm = utm_coords[,2],
         dop = case_when(is.nan(lat_wgs84) ~ NA,
                         .default = dop)) %>% 
  mutate(lon_utm = case_when(is.nan(lon_utm) ~ NA,
                             .default = lon_utm),
         lat_utm = case_when(is.nan(lat_utm) ~ NA,
                             .default = lat_utm)) %>%
  relocate(c(lat_utm, lon_utm), .after=lon_wgs84) %>%
  left_join(collar_meta, by=join_by(collar_id))
  

#write_csv(final, "data/Location_Data/Source_Files/locations_master/gps_locs_master_10-02-2023.csv")


#create a file detailing the source file for each deployment
deployment_sources <- trimmed_all %>% distinct(deployment_id, source)

#write_csv(deployment_sources, "data/Location Data/Metadata/deployment_sources_5-16-2023.csv")


###############################################################################  