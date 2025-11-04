##### Vectronics Data Formatting #####


#Author: Read Barbee

#Creation Date: 2022-10-04

#Last updated: 2022-10-04

#fixed extract_deployments() function from Vectronics_Web_Data_Formattingv1

#Purpose: format complete location download from vectronic for combination with vectronic web data and collar downloads. Get everything into Movebank format.

##Steps:
#1. Import all collar locations from Vectronics and deployment list from local file
#1b (optional): generate summary counts of collars and deployments to compare between web and local list
#2. Filter Vectronics data for collars in deployment list
#3. Filter out retrieved collars? (could do this via matching after all data are combined) Note: figure out if retrieved collars contain all data on site and more or if there is data on the site that's not on the collars

#4. Trim each track by deployment dates and assign to individual based on deployment list

#5. Combine all deployments in to single dataframe


library(tidyverse)
library(lubridate)
library(janitor)
library(collar)


################################## Step 1: Data Import ###########################

#1. import csv from complete vectronics site download
vec_complete_raw <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Vectronics/Vectronic_complete_download_raw_2022-09-27.csv") %>% clean_names()

#2. import cougar deployment list from local file and 
cougar_deployments <- read_csv("data/Location Data/OCP_Cougar_Deployments_9-30-22.csv")%>% clean_names()

#3. format dates
cougar_deployments <- cougar_deployments %>% 
  mutate(deployment_date = mdy(deployment_date),
         end_date = mdy(end_date))


################################## Step 2 & 3: Filter Vectronic Data #######################

#1. Create vector of IDs for retrieved collars
retrieved_collars <- cougar_deployments %>% 
  filter(file_source == "Collar") %>% 
  select(collar_id) %>% 
  unique() %>% 
  unlist() %>% 
  as.vector()


#2. Filter vectronic data only for collars in deployment list and remove data for retrieved collars
vec_web_filtered <- vec_complete_raw %>% 
  filter(collar_id %in% cougar_deployments$collar_id) %>% 
  filter(!(collar_id %in% retrieved_collars))


#3. subset vectronic collars from deployment list
cougar_deployments_vec <- cougar_deployments %>% 
  filter(collar_brand == "Vectronic")


#4. subset collar deployments not on the website
not_on_web<- setdiff(cougar_deployments_vec$collar_id, vec_web_filtered$collar_id)


#5. remove collar deployments not on the website from deployment list for trimming
cougar_deployments_vec_web <- cougar_deployments_vec %>% 
  filter(!(collar_id %in% not_on_web))

#make sure number of deployments in data match number in local list
setdiff(cougar_deployments_vec_web$collar_id, vec_web_filtered$collar_id)

##################### Step 4: Trim Tracks by Deployment Dates #######################

#Create function to extract deployments based on collar ID and deployment dates and append columns of animal_id and deployment_id (Working)

extract_deployments <- function(data_p, animal_id_p, collar_id_p, start_date_p, end_date_p) {
  if (is.na(end_date_p) == TRUE){ #don't filter by end date if it's not included
    trimmed_track <- data_p %>% 
      filter(collar_id == collar_id_p) %>% 
      filter(acq_time_utc>= as.POSIXct(start_date_p)) %>% 
      mutate(animal_id = animal_id_p,
             deployment_id = paste(animal_id_p,"_",collar_id_p)) %>% 
      select(deployment_id, animal_id, everything())
  } else{
    trimmed_track <- data_p %>% 
      filter(collar_id == collar_id_p) %>% 
      filter(acq_time_utc>= as.POSIXct(start_date_p) & acq_time_utc<=as.POSIXct(end_date_p)) %>% 
      mutate(animal_id = animal_id_p,
             deployment_id = paste(animal_id_p,"_",collar_id_p)) %>% 
      select(deployment_id, animal_id, everything())
  }
  
  return(trimmed_track)
}




#create empty list to fill with data frames of trimmed individual tracks
trimmed_deployments <- list()

#loop to generate trimmed tracks for each deployment and append them to list (working)
for (i in 1:nrow(cougar_deployments_vec_web)) {
  trimmed_deployments[[i]]<- extract_deployments(
    data_p = vec_web_filtered,
    animal_id_p = cougar_deployments_vec_web$name[[i]],
    collar_id_p = cougar_deployments_vec_web$collar_id[[i]],
    start_date_p = cougar_deployments_vec_web$deployment_date[[i]],
    end_date_p = cougar_deployments_vec_web$end_date[[i]])
}



####################Step 5: Combine Trimmed Tracks to Single DF #######################

vec_web_final <- bind_rows(trimmed_deployments) %>% 
  mutate(deployment_id = str_replace_all(deployment_id, fixed(" "), ""))


#write to csv in original format (optional)

#write_csv(vec_web_final, "/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/vec_web_trimmed_combined_2022-10-03.csv")


##################### Step 6: Convert dataframe to standard format ###############3

#Standard fields:
#deployment_id
#animal_id
#collar_id
#date_time_gmt
#latitude
#longitude
#altitude_m
#dop
#fix_type
#temp_c
#main_v
#back_v

vec_web_formatted <- vec_web_final %>% 
  select(-c(scts_utc, origin, mortality_status)) %>% 
  rename(date_time_gmt = acq_time_utc,
         latitude = latitude_deg,
         longitude = longitude_deg,
         back_v = backup_v)


#make sure deployments in data match local list
setdiff(cougar_deployments_vec_web$collar_id, vec_web_formatted$collar_id)


#write_csv(vec_web_formatted, "/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/vec_web_final_2022-10-04.csv")





#Future Improvements:
#1. Add functionality to scrape data from vectronic automatically



