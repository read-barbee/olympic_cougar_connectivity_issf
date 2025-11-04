#### Identify dispersals by NSD ####

# Author: Read Barbee

# Date:2023-10-02

#Purpose: Identify dispersing individuals and the timestamp of the dispersal date using Net Squared Displacement (NSD) plots

#Inputs:
# 1.Cleaned and combined gps locations

#Outputs: 
# 1. Manually enter dispersal classes and timestamps into a spreadsheet


#Note: It is harder to determine the dispersal date for females because they don't disperse as far


################################ Libraries #################################
library(tidyverse)
library(sf)
library(amt)
library(mapview)
library(plotly)

################################ Import screened location data #################################

locs <- read_csv("data/gps_data/master_locations/gps_locs_dop_disp_screened_02-13-2024.csv")


################################ Set up animal indexing #################################

#Get vector of unique animal ids
animal_ids <- unique(locs$animal_id)

#Increment this value to cycle through names
current_id <- animal_ids[1]


################################ Inspect locations one individual at a time #################################

#Filter locations to individual of interest
filt <- locs %>% 
  filter(animal_id==current_id) 

#View locs before making track
mapview(filt %>%  sf::st_as_sf(coords = c("lon_utm", "lat_utm"), crs = 5070))

#Make track
track <- make_track(filt, lon_utm, lat_utm, date_time_utc, check_duplicates = TRUE, all_cols = TRUE, crs = 5070) 

################################ Plot net squared displacement #################################

#calculate net squared displacement from first location over time
track <- add_nsd(track)

#plot nsd over time to identify dispersal date 
nsd_plot <- ggplot(data = track, 
                   aes(x = t_, y = nsd_)) +
  geom_point()


#use plotly to find the dispersal date interactively
ggplotly(nsd_plot)


