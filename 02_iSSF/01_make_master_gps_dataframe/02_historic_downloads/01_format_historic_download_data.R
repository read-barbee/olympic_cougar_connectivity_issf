##### Collar Download Data Formatting #####


#Author: Read Barbee

#Date: 2023-05-11


#Purpose: Filter for collared dates and combine and format files from retrieved collar downloads. Get everything into Movebank format.

library(tidyverse)
library(lubridate)
library(janitor)



#Import historical downloads grouped in the 4 formats
hist_form1 <- list.files(
    path="/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Historic Downloads/format1", 
    pattern = "*.csv",
    full.names = TRUE) %>% 
  map_df(~read_csv(.)) %>% 
  clean_names()

hist_form2 <- list.files(
  path="/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Historic Downloads/format2", 
  pattern = "*.csv",
  full.names = TRUE) %>% 
  map_df(~read_csv(.)) %>% 
  clean_names()

hist_form3 <- list.files(
  path="/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Historic Downloads/format3", 
  pattern = "*.csv",
  full.names = TRUE) %>% 
  map_df(~read_csv(.)) %>% 
  clean_names()

hist_form4 <- list.files(
  path="/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Historic Downloads/format4", 
  pattern = "*.csv",
  full.names = TRUE) %>% 
  map_df(~read_csv(., col_types = list(NAV = col_character()))) %>% 
  clean_names()

hist_form5 <- list.files(
  path="/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Historic Downloads/format5/", 
  pattern = "*.csv",
  full.names = TRUE) %>% 
  map_df(~read_csv(.)) %>% 
  clean_names()

hist_form6 <- list.files(
  path="/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Historic Downloads/format6/", 
  pattern = "*.csv",
  full.names = TRUE) %>% 
  map_df(~read_csv(.)) %>% 
  clean_names()

hist_form7 <- list.files(
  path="/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Historic Downloads/format7/", 
  pattern = "*.csv",
  full.names = TRUE) %>% 
  map_df(~read_csv(.)) %>% 
  clean_names()


makah_form <- list.files(
  path="/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Historic Downloads/makah_format/", 
  pattern = "*.csv",
  full.names = TRUE) %>% 
  map_df(~read_csv(., col_types = list(Sats_Used = col_character()))) %>% 
  clean_names()

pre_formatted <- list.files(
  path="/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Historic Downloads/pre-formatted/", 
  pattern = "*.csv",
  full.names = TRUE) %>% 
  map_df(~read_csv(., col_types = list(fix_type = col_character()))) %>% 
  clean_names()



#######Standardize formatting of historical downloads#########
#Standard fields:
#deployment_id
#animal_id
#collar_id
#date_time_utc
#latitude
#longitude
#altitude
#dop
#fix_type

#use force_tz to avoid parsing failures due to daylight savings time


f1_formatted <- hist_form1 %>% 
  mutate(collar_id= as.character(collar_id),
         date_time_utc = mdy_hms(utc_date_time, tz= "UTC"),
         date_time_local = with_tz(date_time_utc, tzone="US/Pacific"),
         deployment_id = paste0(animal_id_2,"_",collar_id),
         fix_type = case_when(fix_type %in% c("val. GPS-3D", "GPS-3D") ~ "3D",
                              fix_type %in% c("GPS-2D") ~ "2D",
                              fix_type %in% c("No Fix") ~ NA_character_)) %>%
  select(deployment_id, animal_id_2, collar_id, date_time_utc, date_time_local, latitude:fix_type) %>% 
  rename(animal_id = animal_id_2,
         altitude = height_m)



f2_formatted <- hist_form2 %>%
  mutate( animal_id = str_replace(animal_id, " ", "_"),
          collar_id=as.character(collar_id),
          date_time_utc = mdy_hms(paste0(date_gmt, " ", time_gmt), tz="UTC"),
          date_time_local = with_tz(date_time_utc, tzone="US/Pacific"), 
          deployment_id = paste0(animal_id,"_",collar_id),
          fix_type= case_when(fix_status %in% c("3D Fix-V", "3D Fix", "4 or more SV KF", "3-D least-squares") ~ "3D",
                              fix_status %in% c("2D Fix", "2-D least-squares", "3_SV KF") ~ "2D",
                              fix_status %in% c("No Sats","0", "DR") ~ NA_character_)) %>%
  select(deployment_id, animal_id, collar_id, date_time_utc, date_time_local, latitude:altitude, dop, fix_type)


f3_formatted <- hist_form3 %>% 
  mutate(animal_id = str_replace(animal_id, " ", "_"),
         collar_id=as.character(collar_id),
         date_time_utc= ymd_hms(date_time_gmt, tz="UTC"),
         date_time_local = ymd_hms(date_time_local, tz="US/Pacific"),
         deployment_id = paste0(animal_id,"_",collar_id),
         fix_type=case_when(fix_status %in% c("3D-V Fix", "3D Fix") ~ "3D",
                            fix_status %in% c("2D Fix") ~ "2D",
                            fix_status %in% c("No Sats") ~ NA_character_)) %>%
  select(deployment_id, animal_id, collar_id, date_time_utc, date_time_local:altitude, dop, fix_type )

  

f4_formatted <- hist_form4 %>% 
  mutate(collar_id= as.character(collar_id),
         date_time_utc = mdy_hms(paste0(gmt_date," ",gmt_time, tz= "UTC")),
         date_time_local = with_tz(date_time_utc, tzone= "US/Pacific"),
         deployment_id = paste0(animal_id,"_",collar_id),
         fix_type = case_when(nav == "3D" ~ "3D",
                              nav == "2D" ~ "2D",
                              nav =="No" ~ NA_character_)) %>%
  rename(altitude = height) %>% 
  select(deployment_id, animal_id, collar_id, date_time_utc, date_time_local, latitude, longitude, altitude, dop, fix_type)




f5_formatted <- hist_form5 %>% 
  mutate(collar_id= as.character(collar_id),
         date_time_utc = dmy_hms(paste0(gmt_date_dmy," ",gmt_time, tz= "UTC")),
         date_time_local = force_tz(dmy_hms(paste0(lmt_date_dmy, " ", lmt_time)), tzone="US/Pacific"),
         deployment_id = paste0(animal_id,"_",collar_id),
         fix_type = case_when(nav == "3D" ~ "3D",
                              nav == "2D" ~ "2D",
                              nav =="No" ~ NA_character_)) %>%
  rename(altitude = height) %>% 
  select(deployment_id, animal_id, collar_id, date_time_utc, date_time_local, latitude, longitude, altitude, dop, fix_type)

# use force_tz for dates considered invalid on the day of daylight savings change in pacific time zone
# using utc for lmt columng because the existing lmt times are the same as utc in the file %>% 
f6_formatted <- hist_form6 %>% 
  mutate(collar_id = as.character(collar_id),
         date_time_utc = mdy_hms(paste0(utc_date," ",utc_time), tz="UTC"),
         date_time_local = with_tz(mdy_hms(paste0(utc_date," ",utc_time)), tzone="US/Pacific"),
         deployment_id = paste0(animal_id,"_",collar_id),
         fix_type = case_when(
           fix_type %in% c("val. GPS-3D", "GPS-3D") ~ "3D",
           fix_type == "GPS-2D" ~ "2D",
           fix_type == "No Fix" ~ NA_character_)) %>% 
  rename(altitude = height_m) %>% 
  select(deployment_id, animal_id, collar_id, date_time_utc, date_time_local, latitude, longitude, altitude, dop, fix_type)


# using utc for lmt columng because the existing lmt times are the same as utc in the file
f7_formatted <- hist_form7 %>% 
  mutate(collar_id = as.character(collar_id),
         date_time_utc = mdy_hms(paste0(utc_date," ",utc_time), tz="UTC"),
         date_time_local = with_tz(mdy_hms(paste0(utc_date," ",utc_time)), tzone="US/Pacific"),
         deployment_id = paste0(animal_id,"_",collar_id),
         fix_type = case_when(
           fix_type %in% c("val. GPS-3D", "GPS-3D") ~ "3D",
           fix_type == "GPS-2D" ~ "2D",
           fix_type == "No Fix" ~ NA_character_)) %>% 
  rename(altitude = height_m) %>% 
  select(deployment_id, animal_id, collar_id, date_time_utc, date_time_local, latitude, longitude, altitude, dop, fix_type)


makah_formatted <- makah_form %>% 
  mutate(collar_id= as.character(collar_id),
         date_time_utc = mdy_hms(paste0(gmt_date," ",gmt_time, tz= "UTC")),
         date_time_local = force_tz(mdy_hms(paste0(lmt_date, " ", lmt_time)), tzone="US/Pacific"),
         deployment_id = paste0(animal_id,"_",collar_id),
         fix_type = case_when(nav == "3" ~ "3D",
                              nav == "2" ~ "2D")) %>%
  rename(altitude = elevation,
         latitude = latitude_wgs84,
         longitude = longitude_wgs84) %>% 
  select(deployment_id, animal_id, collar_id, date_time_utc, date_time_local, latitude, longitude, altitude, dop,fix_type)
  

pre_formatted2 <- pre_formatted %>% 
  mutate(collar_id = as.character(collar_id),
         lmt_date_time = with_tz(lmt_date_time, tzone="US/Pacific"),
         deployment_id =  paste0(animal_id,"_",collar_id)) %>% 
  rename(date_time_utc = utc_date_time,
         date_time_local = lmt_date_time) %>% 
  select(deployment_id, animal_id, collar_id, date_time_utc, date_time_local, latitude, longitude, altitude, dop,fix_type)
  


########## combine all formatted dfs into single df #############
hist_combined <- bind_rows(f1_formatted, 
                           f2_formatted, 
                           f3_formatted, 
                           f4_formatted, 
                           f5_formatted, 
                           f6_formatted,
                           f7_formatted,
                           makah_formatted,
                           pre_formatted2)


get_dupes(hist_combined, animal_id, collar_id, date_time_utc) #no_dupes, we good

#check deployments (n=33)
hist_combined %>% distinct(deployment_id)

write_csv(hist_combined, "/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Source Files/hist_source_5-16-2023.csv")




