##### Combine Web Data #####


#Author: Read Barbee

#Creation Date: 2023-05-11
#Last Updated: 2023-05-11



#Purpose: Format and combine all lotek and vectronic web data

#Steps:
#Step 1: Import and format web location data
#Step 2: Import and format deployment metadata
#Step 3: Filter location data using deployment metadata
#Step 4: Trim Tracks by deployment dates
#Step 5: Remove downloaded deployments from web data to be replaced by download files later


#Note: uses acquisiiton time instead of scts for vectronic because they are earlier and represent the actual time the collar recorded the data.

#Libraries
library(tidyverse)
library(lubridate)
library(janitor)


################################## Step 1: Import and format web location data ###########################

#Import all lotek web data and format columns
lotek_complete_raw <- read_csv("data/Location Data/Raw Data/Lotek/Lotek_complete_download_raw_2023-05-16.CSV") %>% 
  clean_names() %>% 
  mutate(device_id = as.character(device_id),
         date_time_local = force_tz(date_time_local, tzone="US/Pacific"),
         fix_type= case_when(fix_status %in% c("3D Fix-V", "3D Fix", "4 or more SV KF", "3-D least-squares") ~ "3D",
                             fix_status %in% c("2D Fix", "2-D least-squares", "3_SV KF") ~ "2D",
                             fix_status %in% c("No Sats","0", "DR") ~ NA)) %>% 
  rename(collar_id = device_id,
         date_time_utc= date_time_gmt) %>% 
  select(-c(device_name, fix_status, temp_c:back_v)) %>% 
  relocate(fix_type, .after = dop)


#Import all vectronic data and format columns
vec_complete_raw <- read_csv("data/Location Data/Raw Data/Vectronics/Vectronic_complete_download_raw_2023-05-16.csv") %>%
  clean_names() %>% 
  #filter(mortality_status == "Nothing Detected" | is.na(mortality_status==TRUE)) %>% #remove mortality points
  select(-c(scts_utc, origin, mortality_status:animal)) %>% 
  rename(date_time_utc = acq_time_utc,
         latitude = latitude_deg,
         longitude = longitude_deg,
         altitude = altitude_m) %>% 
  mutate(date_time_local = with_tz(date_time_utc, tzone="US/Pacific"), .after= date_time_utc,
         collar_id = as.character(collar_id),
         fix_type=case_when(fix_type %in% c("3D", "3D Validated") ~ "3D",
                            fix_type %in% c("2D") ~ "2D",
                            fix_type %in% c("No Fix") ~ NA))


#combine lotek and vectronic data into single dataframe
web_raw <- bind_rows(lotek_complete_raw, vec_complete_raw)


################################## Step 2: Import and format deployment metadata ###########################

#Import cougar deployment list from local file and format columns
cougar_deployments <- read_csv("data/Location Data/Metadata/From Teams/Formatted for R/collar_deployments_master_5-11-2023.csv")%>%
  clean_names() %>% 
  mutate(collar_id = as.character(collar_id),
         start_date = mdy(start_date),
         end_date = mdy(end_date),
         deployment_id = paste0(name,"_",collar_id), .before=name) 

#Import list of deployments with downloaded data
downloaded_deployments <- read_csv("data/Location Data/Metadata/collar_download_deployment_list_5-11-2023.csv") %>% pull(deployment_id)

################################## Step 3: Filter location data using deployment metadata ###########################

web_collars <- web_raw %>% distinct(collar_id) %>% pull()
dep_collars <- cougar_deployments %>% distinct(collar_id) %>% pull()

setdiff(web_collars, dep_collars) #77 collar ids are in web data but not in deployment list
setdiff(dep_collars, web_collars) #31 collar ids are in deployment list but not in web data

#Filter location data only for collars in deployment list: removes 77 collar IDs from GPS data. Includes collars from previous projects, bobcat collars, and undeployed collars
web_filtered <- web_raw %>% 
  filter(collar_id %in% dep_collars) 

setdiff(web_raw$collar_id, web_filtered$collar_id)


#Filter deployment list to collars in filtered web data for looping: removes 31 collar IDs, 34 deployments (because of redeployed collars)
dep_filtered <- cougar_deployments %>% 
  filter(collar_id %in% web_filtered$collar_id)

setdiff(cougar_deployments$collar_id, dep_filtered$collar_id)


################################# Step 4: Trim Tracks by deployment dates #################################

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
for (i in 1:nrow(dep_filtered)) {
  trimmed_deployments[[i]]<- extract_deployments(
    data_p = web_filtered,
    animal_id_p = dep_filtered$name[[i]],
    collar_id_p = dep_filtered$collar_id[[i]],
    start_date_p = dep_filtered$start_date[[i]],
    end_date_p = dep_filtered$end_date[[i]])
}

#combine all trimmed deployments into single dataframe
web_all <- bind_rows(trimmed_deployments)


################################# Step 5: Remove downloaded deployments from web data #################################

#remove downloaded deployments from web file to avoid duplicating points (should remove 12 deployments)
web_final <- web_all %>% 
  filter(!(deployment_id %in% downloaded_deployments))

web_final %>% distinct(deployment_id) %>% count() #web final has 121 distinct deployments instead of previous 125 because of additional downloaded deployments received from Kim on 5-15-2023.

#check to make sure only the 17 downloaded collars were removed
setdiff(web_all$deployment_id, web_final$deployment_id)

#only 15 out of 17 removed becaue Lilu 87524 and Gypsy 3468 not on web


#export to csv
#write_csv(web_final, "data/Location Data/Source Files/web_source_5-16-2023.csv")




#make sure that all collar IDs not on the web are accounted for in historic downloads:

historic_downloads <- c("28360", "82112", "82111", "80464", "80462", "81660", "80463", "32156", "33432", "32153", "32248", "38018", "32244", "38015", "32246", "32245", "38019", "28358", "28359", "34639", "28361", "28362","99999991", "99999992", "99999993", "99999994", "99999995", "99999996", "99999997", "99999998" )

not_on_web <- setdiff(dep_collars, web_collars)

setdiff(not_on_web, historic_downloads)

#the only diff is 34638, which is in the collar downloads



