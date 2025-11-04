##### Collar Download Data Formatting #####


#Author: Read Barbee

#Creation Date: 2022-10-03
#Updated: 2023-05-16

#Purpose: Filter for collared dates and combine and format files from retrieved collar downloads. 

library(tidyverse)
library(lubridate)
library(janitor)
#library(collar)


################################ Import all downloaded collar files #################################
format1 <-
  list.files(
    path="data/Location Data/Raw Data/Collar Downloads/Lotek/formatted/format1", 
    pattern = "*.csv",
    full.names = TRUE) %>% 
  map_df(~read_csv(.)) %>% 
  clean_names() 

format2 <-
  list.files(
    path="data/Location Data/Raw Data/Collar Downloads/Lotek/formatted/format2", 
    pattern = "*.csv",
    full.names = TRUE) %>% 
  map_df(~read_csv(.)) %>% 
  clean_names() 

format3 <-
  list.files(
    path="data/Location Data/Raw Data/Collar Downloads/Lotek/formatted/format3", 
    pattern = "*.csv",
    full.names = TRUE) %>% 
  map_df(~read_csv(.)) %>% 
  clean_names() 

#Import cougar deployment list for date filtering
cougar_deployments <- read_csv("data/Location Data/Metadata/From Teams/Formatted for R/collar_deployments_master_5-11-2023.csv") %>% 
  clean_names() %>% 
  mutate(deployment_id = paste0(name,"_",collar_id), .before=name)

downloads <- read_csv("data/Location Data/Metadata/collar_download_deployment_list_5-11-2023.csv") %>% pull(deployment_id)

downloaded_deployments <- cougar_deployments %>% 
  filter(deployment_id %in% downloads)


################################ Convert all files to common format #################################

f1_standard <- format1 %>% 
  mutate(collar_id=as.character(collar_id),
         date_time_utc = mdy_hm(gmt_time, tz="UTC"),
         date_time_local = with_tz(date_time_utc, tzone="US/Pacific"),
         fix_type = case_when(satellites >= 4 ~ "3D",
                              satellites == 0 ~ NA_character_,
                              satellites > 0 & satellites < 4 ~ "2D")) %>%
  select(deployment_id:collar_id, date_time_utc, date_time_local, latitude:altitude, dop, fix_type)

f2_standard <- format2 %>% 
  mutate(collar_id=as.character(collar_id),
         date_time_utc = mdy_hm(date_time_gmt, tz="UTC"),
         date_time_local = with_tz(date_time_utc, tzone="US/Pacific"),
         fix_type = case_when(fix_status %in% c("3-D least-squares", "4 or more SV KF") ~ "3D",
                              fix_status %in% c("2-D least-squares", "3_SV KF") ~ "2D",
                              is.na(fix_status)==TRUE ~ NA_character_)) %>%
  select(deployment_id:collar_id, date_time_utc, date_time_local, latitude:altitude, dop, fix_type)

f3_standard <- format3 %>% 
  mutate(collar_id=as.character(collar_id),
         date_time_utc = mdy_hms(paste0(utc_date," ",utc_time), tz="UTC"),
         date_time_local = with_tz(date_time_utc, tzone="US/Pacific"),
         fix_type = case_when(fix_type %in% c("val. GPS-3D", "GPS-3D") ~ "3D",
                              fix_type == "GPS-2D" ~ "2D",
                              fix_type == "No Fix" ~ NA_character_)) %>%
  rename(altitude = height_m) %>% 
  select(deployment_id:collar_id, date_time_utc, date_time_local, latitude, longitude, altitude, dop, fix_type)




#Filter tracks by deployment times
#Create function to extract deployments based on collar ID and deployment dates and append columns of animal_id and deployment_id (Working)

extract_deployments <- function(data_p, animal_id_p, collar_id_p, start_date_p, end_date_p) {
  if (is.na(end_date_p) == TRUE){ #don't filter by end date if it's not included
    trimmed_track <- data_p %>% 
      filter(collar_id == collar_id_p) %>% 
      filter(date_time_local>= as.POSIXct(start_date_p)) %>% 
      mutate(animal_id = animal_id_p,
             deployment_id = paste(animal_id_p, "_", collar_id_p)) %>% 
      select(deployment_id, animal_id, everything())
  } else{
    trimmed_track <- data_p %>% 
      filter(collar_id == collar_id_p) %>% 
      filter(date_time_local>= mdy(start_date_p) & date_time_local<=mdy(end_date_p))
  }
  
  return(trimmed_track)
}

#comine formatted files into single dataframe
collar_downloads_all <- bind_rows(f1_standard,
                                  f2_standard,
                                  f3_standard)

#create empty list to fill with data frames of trimmed individual tracks
trimmed_deployments <- list()



#loop to generate trimmed tracks for each deployment and append them to list (working)
for (i in 1:nrow(downloaded_deployments)) {
  trimmed_deployments[[i]]<- extract_deployments(
    collar_downloads_all,
    downloaded_deployments$name[[i]],
    downloaded_deployments$collar_id[[i]],
    downloaded_deployments$start_date[[i]],
    downloaded_deployments$end_date[[i]])
}


#recombine trimmed tracks into single dataframe

retrieved_collars_final <- bind_rows(trimmed_deployments) 


#write_csv(retrieved_collars_final, "/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Source Files/collar_source_5-16-2023.csv")

#final file contains all deployments in downloaded list
setdiff(retrieved_collars_final$deployment_id, downloads)
setdiff(downloads, retrieved_collars_final$deployment_id)


