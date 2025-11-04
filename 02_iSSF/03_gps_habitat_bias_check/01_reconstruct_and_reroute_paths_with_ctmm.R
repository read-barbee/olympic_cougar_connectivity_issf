#### Reconstruct and Reroute Paths with CTMM ####

# Author: Read Barbee

# Date:2023-07-12 

# Purpose: Impute missing GPS locations using a continuous time movement model

#Inputs:
#1. Screened gps locations

#Outputs: 
#Paths with imputed locations


################################ Libraries #################################
library(tidyverse)
library(terra)
library(amt)
library(janitor)
library(furrr)
library(sf)
library(mapview)
library(ctmm)

################################ Parameters #################################

#number of cores to use for ctmm path imputation
cores = 8

#########################################################################
##
## 1. Import data
##
##########################################################################
# screened gps locations
locs_screened <- read_csv("data/gps_data/master_locations/screened/gps_locs_dop_disp_screened_02-13-2024.csv") %>%
  mutate(date_time_local = with_tz(date_time_local, tzone = "US/Pacific"))

#water body polygons for path rerouting
water_polys <- st_read("data/spatial/op_water_polys_with_salt.gpkg")


#########################################################################
##
## 2. Convert locations to amt tracks
##
##########################################################################

#Nest locations by animal_id, sex, and dispersal status
locs_nested <- locs_screened %>% 
  nest_by(animal_id, sex, dispersal_status) 


#create tracks for each animal
locs_nested$tracks <-map(locs_nested$data, 
                         function(x){
                           make_track(x, 
                                      lon_utm, 
                                      lat_utm, 
                                      date_time_utc, 
                                      check_duplicates = TRUE, 
                                      all_cols = TRUE, 
                                      crs = 5070)})

#########################################################################
##
## 3. Examine Sampling Rates
##
##########################################################################

#calculate sampling rate for each individual
locs_nested$sr <-  map(locs_nested$tracks, summarize_sampling_rate)

#calculate median sampling rate
locs_nested$median_sr <-  map(locs_nested$tracks, 
                              function(x){
                                summ <- summarize_sampling_rate(x)
                                med <- round(summ$median)
                                return(med)})


#calculate sampling rate statistics by individual
sampling_rates <- locs_nested %>% 
  select(-c(data, tracks)) %>% 
  unnest(cols = c(sr)) 


#########################################################################
##
## 4. Reconstruct tracks with ctmm (~40 min 8 clusters)
##
##########################################################################

#create telemetry object to fit ctmm models
ctmm_dat <- locs_nested %>% 
  unnest(cols=c(tracks)) %>%
  mutate(t_ = round_date(t_, unit="hour")) %>% 
  dplyr::select(animal_id, collar_id, t_, lat_wgs84,lon_wgs84) %>% 
  rename(individual.local.identifier = animal_id,
         tag.local.identifier = collar_id,
         timestamp = t_,
         location.long = lon_wgs84,
         location.lat = lat_wgs84) %>% 
  as.telemetry(datum = "+proj=longlat +datum=WGS84 +no_defs +type=crs",
               projection = "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs" )


#initialize furrr
plan(multisession, workers = cores)

#Fit ctmm model and impute points for each individual (1 hour with 9 cores)

system.time(ctmm_sims <- future_map(1:length(ctmm_dat), function(i){
  
  #initialize list for simulations
  ctmm_sims <- list()
  
  #calculate initial values for ctmm model
  guess <- ctmm.guess(ctmm_dat[[i]], interactive = FALSE)
  
  #perform ctmm model selection
  ctmm_fits<- ctmm.select(ctmm_dat[[i]], guess, verbose=TRUE, cores=8)
  
  #select top fitted ctmm model
  mod <- ctmm_fits[[1]]
  
  #simulate missing points based on the model and the observed data (simulation runs through observed pts)
  ctmm_sim <- ctmm::simulate(object = mod, data = ctmm_dat[[i]], complete=TRUE)
  
  #classify points as observed or imputed and assign unique group id to consecutive locs of same type
  sims_marked <- ctmm_sim %>% as_tibble() %>% 
    mutate(imp_status = case_when(ctmm_sim$timestamp %in% ctmm_dat[[i]]$timestamp ~ "observed", .default = "imputed")) %>% 
    mutate(group = data.table::rleid(imp_status))
  
  #count the number of consecutive values in each group
  group_counts <- sims_marked %>% group_by(group) %>% count()
  
  max_streak <- 24 / locs_nested$median_sr[[i]]
  
  #join the group counts to the original data frame and filter out groups of consecutive imputed points spanning 24 hr or more 
  sims_marked <- sims_marked %>% 
    left_join(group_counts, by=join_by(group)) %>% 
    filter(imp_status =="observed" | n < max_streak) %>% 
    mutate(animal_id = ctmm_dat[[i]]@info$identity, .before = t) %>% 
    select(-c(group, n))
  
  ctmm_sims[[i]] <- sims_marked
}, .progress = T))


#bind simulated frames together
ctmm_sims_unrouted <- bind_rows(ctmm_sims)

#export
#write_csv(ctmm_sims_unrouted, "data/Location_Data/imputed_paths/ctmm_sim_imputed_paths_unrouted_12-04-2023.csv")


## map imputed tracks for quality control
og_telem_sf <- ctmm_dat[[1]] %>% as_tibble() %>% st_as_sf(coords = c("x", "y"), crs = 5070) %>% 
  mutate(type= "og")

sim_test_sf <- ctmm_sims[[1]] %>% as_tibble() %>% st_as_sf(coords = c("x", "y"), crs = 5070) %>% 
  mutate(type ="sim")

mapview::mapview(sim_test_sf, zcol="imp_status")
mapview::mapview(bind_rows(sim_test_sf, og_telem_sf), zcol="type")

#Exclude Comet(23) and Cato(20). (too few points)

#########################################################################
##
## 5. Reroute imputed paths around major water bodies (~ 20 min 8 clusters)
##
##########################################################################

#filter water body polygons to only include permanent water bodies
water_polys_filtered <- water_polys %>% filter(WB_PERIOD_ =="PER" )


#For loop to reroute all the imputed paths from ctmm (~ 21 min with 9 cores)
system.time(sims_rerouted <- future_map(1:length(ctmm_sims), function(i){
  
  #name the iteration? Not sure what this does
  cat(i," ")
  
  #convert ctmm sim object to sf_object for pathroutr
  sim_sf <- ctmm_sims[[i]] %>% 
    as_tibble() %>% 
    st_as_sf(coords = c("x", "y"), crs = 5070, remove = FALSE)
  
  #Filter water body polygons to retain those within convex hull of all of individual's points
  water_polys_cropped <- sf::st_buffer(sim_sf, dist = 10000) %>%
    sf::st_union() %>%
    sf::st_convex_hull() %>%
    sf::st_intersection(water_polys_filtered) %>%
    st_collection_extract('POLYGON') %>%
    st_union() %>% 
    st_sf()
  
  #create buffer around barrier objects as visgraph for rerouting function (connects all verticies of barrier polygon with Delaunay triangle mesh and removes any edges that cross the barrier). Essentially it creates a roadmap of traversible terrain
  visgraph <- pathroutr::prt_visgraph(water_polys_cropped, buffer = 15)
  
  sims_rerouted <- list()
  path <- pathroutr::prt_trim(sim_sf, water_polys_cropped)
  sims_rerouted[[i]] <- pathroutr::prt_reroute(path, water_polys_cropped, visgraph) %>% 
    pathroutr::prt_update_points(path) 
}, .progress = T))


#inspect rerouted paths
mapview::mapview(sims_rerouted[[1]], zcol = "imp_status")

#bind for loop iterations into single dataframe
sims_rerouted_full <- bind_rows(sims_rerouted)

#get rerouted coordinates
rerouted_coords <- sims_rerouted_full %>% st_coordinates()

#update coordinates and format for export
sims_rerouted_full_df <- sims_rerouted_full %>% 
  as_tibble() %>%
  mutate(x_rerouted = rerouted_coords[,1],
         y_rerouted = rerouted_coords[,2], .after=t) %>% 
  rename(x_old = x,
         y_old = y,
         longitude_old = longitude,
         latitude_old = latitude) %>% 
  relocate(timestamp, .before=t) %>% 
  select(-c(geometry, longitude_old, latitude_old, t)) 

#export rerouted paths
#write_csv(sims_rerouted_full_df, "data/gps_data/imputed_paths/ctmm_sim_imputed_paths_rerouted_12-04-2023.csv" )
