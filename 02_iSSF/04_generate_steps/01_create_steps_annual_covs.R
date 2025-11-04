#### OCP iSSF Module_02: Prepare Screened GPS Data in amt ####

# Author: Read Barbee

# Date:2023-07-12 
#Last updated: 2025-10-27

# Purpose: Prepare screened GPS data in amt to fit iSSFs.

# Inputs:
#   •	Screened GPS Locations (corrected for error and habitat bias)
#   •	Deployment list
#   •	Disperser list
#   •	Covariate Stack
#
# Outputs:
#   •	amt_step dataframe with fields for all necessary covariates



################################ Libraries #################################
library(tidyverse)
library(terra)
library(amt)
library(janitor)
#library(doParallel)
library(sf)
library(mapview)
library(flextable)


################################ User-Defined Parameters #################################

project_crs <- 5070 #NAD83 Albers Equal Area Projection

resample_int_hours <- 2 #interval to resample tracks

tolerance_mins <- 45 #tolerance for track resampling in minutes

rand_steps <- 10 #number of random steps to generate per actual step

study_area_buffer_dist <- 1000 #meters. 500m allows for walking out on beaches at low tide

#path to annual covariate raster stacks
cov_folder_path <- "/Volumes/Reads_Seagate/ms_thesis/covariates/annual_covs_01-07-24/30m"

#plan(multisession, workers = 8)

#########################################################################
##
## 1. Import and format Barrier Polygons
##
##########################################################################

## Study area boundary ##
op_poly <- st_read("data/spatial/vector/ocp_study_area_poly_wa_only_10-30-23.gpkg")
op_poly_buffer <- st_buffer(op_poly, study_area_buffer_dist)

## Water Polygons ##
water_polys <- st_read("data/spatial/vector/op_water_polys_with_salt.gpkg") %>% st_transform(crs = project_crs)

## Freshwater only ##
water_polys_cropped <-  water_polys %>%
  filter(WB_PERIOD_ =="PER" | OBJECTID != 1) %>% #only include permanent water bodies with area > 100 m2
  filter(SHAPEAREA >= 100) %>%
  sf::st_crop(op_poly) #crop to study area

## Freshwater only (dissolved) ##
water_polys_mask <- water_polys_cropped %>%
  sf::st_union() %>% #dissolve polygons into single vector mask layer
  st_sf()

## Saltwater only ##
ocean <- water_polys %>% filter(OBJECTID == 1)

#########################################################################
##
## 2. Import and format location data
##
##########################################################################
# Mountain lion location data 
locs_screened <- read_csv("data/gps_data/master_locations/screened/gps_locs_dop_disp_screened_02-13-2024.csv")

#check for duplicate locations
get_dupes(locs_screened, deployment_id, date_time_utc)

#########################################################################
##
## 3. Import and format covariate data
##
##########################################################################

cov_files <- list.files(cov_folder_path, pattern = "\\.tif$", full.names = TRUE)

# Create an empty list to store the raster objects
cov_stacks <- list()

# Loop through each .tif file, import it, and assign the file name as the object name
for (tif_file in cov_files) {
  raster_name <- tools::file_path_sans_ext(basename(tif_file))
  raster <- rast(tif_file)
  cov_stacks[[raster_name]] <- raster
}

#import static stack
static_stack <- rast("/Volumes/Reads_Seagate/ms_thesis/covariates/1km_buffer/static_stack_1km_buffer_11-29-23.tif")

#########################################################################
##
## 4. Convert locations to amt tracks
##
##########################################################################

#Nest locations by animal_id 
locs_nested <- locs_screened %>% 
  nest_by(animal_id, sex, dispersal_status) 

#create tracks
locs_nested$tracks <-map(locs_nested$data, function(d){
  make_track(d, 
             lon_utm, 
             lat_utm, 
             date_time_utc, 
             check_duplicates = TRUE, 
             all_cols = TRUE, 
             crs = project_crs) 
  })

#check for duplicate points in tracks
bind_rows(locs_nested$tracks) %>% get_dupes(deployment_id, t_)

#########################################################################
##
## 5. Screen out long round-trips
##
##########################################################################

#source: USFS mountain lion fact sheet: https://www.fs.usda.gov/Internet/FSE_DOCUMENTS/stelprdb5251229.pdf

#max sprint speed = 50 mi/hr = 80 km/hr
#average walking speed = 10 mi/hr = 16 km/hr

#maximum speed of 80 km/hr sustained for 60 seconds
sdr <- calculate_sdr(speed = 80, time = seconds(60), speed_unit = "km/h")

get_displacement(sdr, hours(2)) #compare to actual observed max step length

#copy locs_nested into new object for screening
locs_nested_trip_screened <- locs_nested

#screen for fast roundtrips
locs_nested_trip_screened$tracks <- map(locs_nested_trip_screened$tracks, function(x){
  flag_roundtrips(x, 
                  delta = sdr, 
                  epsilon = 1, 
                  time_unit = "secs")
})


#this filter doesn't remove any locations
locs_nested_trip_screened %>% select(-data) %>% unnest(cols= "tracks") %>% filter(fast_roundtrip_ == TRUE)


#########################################################################
##
## 6. Examine Sampling Rates
##
##########################################################################

#inspect sampling rate for each individual

locs_nested$sr <-  map(locs_nested$tracks, summarize_sampling_rate)

locs_nested$median_sr <-  map(locs_nested$tracks, function(x){
  summ <- summarize_sampling_rate(x)
  med <- round(summ$median)
  return(med)
})

#sampling rate statistics by individual
sampling_rates <- locs_nested %>% 
  select(-c(data, tracks, median_sr)) %>% 
  unnest(cols = c(sr)) 

summary(sampling_rates)

#write_csv(sampling_rates, "individual_sampling_rates_2-15-24.csv")

sr_tab <- sampling_rates %>% 
  select(animal_id, sex, dispersal_status, median) %>% 
  mutate(median=round(median)) %>% 
  mutate(median2 = case_when(median > 4 ~ 999,
                             median <= 4 ~ median)) %>% 
  group_by(sex, dispersal_status, median2) %>% 
  count() %>% 
  arrange(median2) %>% 
  pivot_wider(names_from = median2, values_from = n) %>% 
  rename(`> 4` = `999`)

# Column sums
col_sums <- sr_tab %>%
  ungroup() %>% 
  select(-c(sex, dispersal_status)) %>% 
  summarise(across(everything(), function(x){sum(x, na.rm=T)})) %>% 
  mutate(sex = "Total", dispersal_status = NA, .before = `1`)

# Row sums 
sr_tab2 <- sr_tab %>% 
  bind_rows(col_sums) |> 
  select(-c(sex, dispersal_status)) %>% 
  rowwise() %>% 
  mutate(Total = sum(c_across(everything()), na.rm=T)) %>% 
  ungroup()

sampling_tab <- flextable(sr_tab2) %>%
  add_header_row(values = c(" ", "Median sampling rate (hours)", " "), colwidths = c(2, 5, 1))
  align(align="center", part = "all")

# save_as_docx("table1" = sampling_tab, path = "gps_median_sampling_rates_by_demography_training_dat_2-15-24.docx")


#########################################################################
##
## 7. Identify optimal resampling interval
##
##########################################################################

#function to calculate bursts for each individual---increment by hours to find optimum

test <- locs_nested

test$steps <- map(locs_nested$tracks, 
                    function(x) {
                      x %>% 
                        amt::track_resample(rate = hours(2), tolerance = minutes(45)) %>% #resample to 2 hours
                        amt::filter_min_n_burst() %>% #divide into bursts of min 3 pts
                        amt::steps_by_burst()
                    })

#make list of individuals to remove (because they can't be resampled at the specified interval)
indiv_to_remove <- test %>% 
  filter(nrow(steps)<3) %>% 
  pull(animal_id)

#check sampling rates of individuals to remove -- all ranging from 3-15 hours
sampling_rates %>% filter(animal_id %in% indiv_to_remove)

#1 hour removes 68 individuals (72 with imp points)
#2 hour removes 33 individuals (34 with imp points)
#3 hours removes 34 individuals (35 with imp points)
#4 hours removes 35 individuals (31 with imp points)
#5 hours removes 66 individuals (62 with imp points)
#6 hours removes 9 individuals (8 with imp points)

#6 hours retains the most individuals
# 2 hours may offer the best balance between data resolution and data loss


#########################################################################
##
## 8. Resample tracks and convert to steps
##
#######################################################################

#get resampled count for table
tracks_resampled <- locs_nested
tracks_resampled$tracks <- map(tracks_resampled$tracks, 
                               amt::track_resample, 
                               rate = hours(resample_int_hours), 
                               tolerance = minutes(tolerance_mins))


#define function to convert locatons to steps. 
steps_calc <- function(x) {
  x %>% 
    amt::track_resample(rate = hours(resample_int_hours), tolerance = minutes(tolerance_mins)) %>% #resample 
    amt::filter_min_n_burst() %>% #divide into bursts of min 3 pts
    amt::steps_by_burst() #convert to steps
}

#make a copy of locs nested for step generation
amt_locs <- locs_nested 

#map the step function over each individual and append as a nested dataframe
amt_locs$steps <- map(amt_locs$tracks, steps_calc)

#select the relevant columns and unnest the steps column
amt_steps <- amt_locs %>% 
  select(animal_id:dispersal_status,
         steps) 

#check distribution of step lengths
amt_steps %>% unnest(cols=steps) %>% pull(sl_) %>% quantile(c(0, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 1))

#check again after removing steps over 20km
amt_steps %>% 
  unnest(cols=steps) %>% 
  filter(sl_  <= 20000) %>% 
  pull(sl_) %>% 
  quantile(c(0, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 1))

#########################################################################
##
## 9. Summarize step length and turn angle by individual and demography
##
#######################################################################


#extract step length and turn angle values for each individual
names <- vector()
sl_vals <- list()
ta_vals <- list()
for(i in 1:nrow(amt_steps)){
  names[i] <- amt_steps$animal_id[i]
  sl_vals[[i]] <- amt_steps$steps[[i]] %>% pull(sl_) %>%  
    quantile(c(0, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 1)) %>% 
    enframe() %>% pivot_wider(names_from = name, values_from = value) %>% 
    rename(sl_0 =  `0%`,
           sl_2.5 = `2.5%`,
           sl_5 = `5%`,
           sl_25 = `25%`,
           sl_50 = `50%`,
           sl_75 = `75%`,
           sl_95 = `95%`,
           sl_975 = `97.5%`,
           sl_100 = `100%`)
  
  ta_vals[[i]] <- amt_steps$steps[[i]] %>% pull(ta_) %>%  
    quantile(c(0, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 1), na.rm=T) %>% 
    enframe() %>% pivot_wider(names_from = name, values_from = value) %>% 
    rename(ta_0 =  `0%`,
           ta_2.5 = `2.5%`,
           ta_5 = `5%`,
           ta_25 = `25%`,
           ta_50 = `50%`,
           ta_75 = `75%`,
           ta_95 = `95%`,
           ta_975 = `97.5%`,
           ta_100 = `100%`)
}

#bind extracted values together
sl_vals <- bind_rows(sl_vals)
ta_vals <- bind_rows(ta_vals)

#create a summary frame of individual movement parameters
indiv_move_sum <- bind_cols(sl_vals, ta_vals) %>% mutate(animal_id = names, .before =everything())

#join movement summary to demographic data
indiv_move_summ_full <- sampling_rates %>%
  select(animal_id, sex, dispersal_status) %>% 
  left_join(indiv_move_sum, by = "animal_id") 

#remove individuals for which sl and ta couldn't be calculated and summarize movement by demography
indiv_move_summ_filtered <- indiv_move_summ_full %>% 
  filter(!(is.na(sl_0))) %>% 
  group_by(sex, dispersal_status) %>% 
  summarize(across(where(is.numeric), median)) %>% 
  mutate(across(contains("sl_"), round)) %>% 
  mutate(across(contains("ta_"), function(x){round(x, 2)})) %>% 
  select(sex, dispersal_status, 
         contains(c("sl_2.5", "sl_50", "sl_975")),
         contains(c("ta_2.5", "ta_50", "ta_975"))) %>% 
  ungroup()

#make flextable
move_sum_tab <- flextable(indiv_move_summ_filtered) %>% 
  add_header_row(values = c(" ", "Step length (m)", "Turn angle (radians)"), 
                 colwidths = c(2, 3, 3)) %>% 
  vline(j=5, part = "body")

#write_csv(indiv_move_summ_filtered, "indiv_sl_ta_quantiles_2hr_2-15-24.csv")

# save_as_docx("table1" = move_sum_tab, path = "sl_ta_demographic_quants_training_dat_2-15-24.docx")


#########################################################################
##
## 10. Remove steps < 100m in length and > 20 km in length
##
#######################################################################
sl_screen <- amt_steps

sl_screen$steps <- map(sl_screen$steps, function(x){
  x %>% 
    filter(sl_ >= 100) %>% 
    filter(sl_ <= 20000)
})

########################################################################
##
## 11. Remove individuals with too little data for inference
##
##########################################################################

#Get vector of individuals with < 30steps and/or < 20 days of steps
removal_list <- sl_screen  %>%
  unnest(cols=steps) %>% 
  ungroup() %>% 
  group_by(animal_id) %>% 
  mutate(date_range=interval(start=min(t1_), end=max(t2_)),
         step_days = as.duration(date_range)/ddays(1)) %>%
  summarize(n_steps=n(), step_days = round(first(step_days), 0)) %>%
  filter(n_steps <30 | step_days < 20) %>%
  pull(animal_id)

#remove individuals with no steps and individuals in the removal list
sl_screen2 <- sl_screen %>%
  filter(nrow(steps) != 0) %>% 
  filter(!(animal_id %in% removal_list))

  

#########################################################################
##
## 12. Generate random steps
##
#######################################################################
#Generate random steps. Won't generate random steps in water or outside the study area. Can only generate random steps for individuals with some bursts with > 3 steps
make_random_steps <- function(x){
  
  burst_counts <- x %>% 
    group_by(burst_) %>% 
    count() %>% 
    filter(n>=3)
  
  if(nrow(burst_counts) <= 1){ return(NA)} #must have at least 2 bursts with >=3 steps 
  
  else{
    
    #generate 10 times the amount of desired available steps
    part2 <- x %>% 
      amt::random_steps(n_control = rand_steps*5) %>% #generate random steps per used step
      #amt::extract_covariates(cov_stack, where = "end") %>% #extract covariates at end of step
      amt::time_of_day(include.crepuscule = FALSE) %>% #calculate time of day for each step--not working
      mutate(unique_step = paste(burst_,step_id_,sep="_"), .before=step_id_) %>%  #add column for unique step id }
      mutate(log_sl_ = log(sl_),    
             cos_ta_ = cos(ta_))
    
    #determine whether each random step intersects a water body and/or is within the study area. this is the bottleneck
    part3 <- part2 %>% 
      st_as_sf(coords = c("x2_", "y2_"), crs = project_crs, remove=FALSE) %>% 
      mutate(intersects_freshwater = lengths(st_intersects(., water_polys_mask)) > 0 ) %>% #this is the bottleneck (water polys mask. works with just ocean or just water polys, but not both together)
      mutate(intersects_ocean = lengths(st_intersects(., ocean)) > 0 ) %>%
      mutate(intersects_water = case_when(intersects_ocean==FALSE & intersects_freshwater==FALSE ~FALSE,
                                          .default=TRUE)) %>%
      mutate(intersects_study_area = lengths(st_intersects(., op_poly_buffer)) > 0) %>%
      dplyr::select(unique_step, t1_, t2_, x2_, y2_, intersects_water, intersects_study_area) %>% 
      as.data.frame() %>% 
      select(-geometry)
    
    #join intersection status to original step frame
    part4 <- part2 %>% 
      left_join(part3, by=join_by(unique_step, t1_, t2_, x2_, y2_))
    
    #split used and random steps to filter the random steps by intersection
    tf_split <- split(part4, part4$case_)
    
    #Keep only random steps within the study area that don't intersect water and randomly sample the desired number from those remaining
    tf_split$`FALSE` <- tf_split$`FALSE` %>% 
      filter(intersects_study_area == TRUE) %>% 
      filter(intersects_water == FALSE) %>% 
      slice_sample(n = rand_steps, by = unique_step)
    
    #put the used and random steps back together in a single frame
    final <- bind_rows(tf_split) %>% arrange(t2_, desc(case_))
    
    return(final)
  }
}

amt_steps_r <- sl_screen2

amt_steps_r$steps <- map(amt_steps_r$steps, make_random_steps)

#get list of individuals for which random steps could not be calculated and remove them
removal_indices <- which(is.na(amt_steps_r$steps))
indivs_to_remove2 <- amt_steps_r[c(removal_indices),] %>% pull(animal_id)

amt_steps_r2 <- amt_steps_r %>% filter(!(animal_id %in% indivs_to_remove2))

#########################################################################
##
## 13. Extract covariate values
##
#######################################################################

amt_steps_all_covs <- amt_steps_r2

#static covariates
amt_steps_all_covs$steps <- map(amt_steps_all_covs$steps, extract_covariates, covariates = static_stack)

#annual covariates (make sure to adjust column indices if more are added)
extract_annual_covs <- function(steps, cov_stack_list){
  #create year column
  steps2 <- steps %>%
    mutate(year = as.factor(year(t1_)), .before = t1_)
  
  #split data by year
  steps_split <- split(steps2, steps2$year)
  
  #define function to select the raster stack for the relevant year to extract from and extract the covariates
  extract_fun <- function(steps_split, year){
    stack <- cov_stack_list[[paste0("cov_stack_", as.character(year))]] # "_water_masked"
    covs <- steps_split %>% extract_covariates(stack)
    
    return(covs)
  }
  
  #apply the extraction function to each year
  steps_split_covs <- list()
  for(i in 1:length(steps_split)){
    names <- names(steps_split)
    steps_split_covs[[i]] <- extract_fun(steps_split[[i]], names[i])
  }
  
  names(steps_split_covs) <- names
  
  remove_year_names <- function(steps_split_covs){
    old_names <- names(steps_split_covs)
    new_names <- vector()
    for(i in 1:length(old_names)){
      if(str_detect(old_names[i], coll("20"))){
        new_names[i] <- substr(old_names[i], 1, nchar(old_names[i]) - 5)
      } else{
        new_names[i] <- old_names[i]
      }
    }
    
    names(steps_split_covs) <- new_names
    
    return(steps_split_covs)
  }
  
  #apply renaming function
  covs_renamed <- map(steps_split_covs, remove_year_names)
  
  #bind all years together 
  steps_covs_final <- bind_rows(covs_renamed)
  
  return(steps_covs_final)
  
}

amt_steps_all_covs$steps <- map(amt_steps_all_covs$steps, extract_annual_covs, cov_stack_list = cov_stacks)

#########################################################################
##
## 14. Interpolate values in water-masked extent of GEE products that are actually on land
##
#########################################################################

#points near bodies of water don't have hii values, npp, gpp, or perc veg values because of coarse water-masking in GEE products. Solve by interpolation

#list of covariates with coarse water masking to interpolate values for
covs_to_interp_300m <- c("perc_nonveg_annual", "perc_nontree_veg_annual", "perc_tree_cov_annual", "hii_annual", "roads_hii_annual", "infra_hii_annual", "landuse_hii_annual", "popdens_hii_annual", "mtpi")

#Interpolate precip, npp and gpp with 600 m radius b/c of coarser scale
covs_to_interp_600m <- c("precip_annual", "npp_annual", "gpp_annual")

#function to interpolate values from specified buffer radius around points
interp_cov_vals <- function(cov, df, buffer_size){

#prep data and add na status column
col_name <- cov 
column <- df[[col_name]] 

#bypass function if there are no NA values
if(sum(is.na(column)) == 0){
  return(df)
} else if(sum(is.na(column)) > 0){
df2 <- df %>% 
  mutate(na_status = is.na(!!sym(col_name))) 

#split data by na_status
na_split <- split(df2, df2$na_status)

#select the na values for interpolation
na_only <- na_split$`TRUE` %>% 
  st_as_sf(coords=c("x2_", "y2_"), crs = project_crs, remove = FALSE) %>% 
  st_buffer(buffer_size)

#split na_values by year
na_year <- split(na_only, na_only$year)

#interpolate missing values as average of cells within buffer radius for each year
for(i in 1:length(na_year)){
 if(nrow(na_year[[i]])>0){
  year <- names(na_year)[i]
  stack <- cov_stacks[[paste0("cov_stack_", as.character(year))]]
  
  if(sum(str_detect(names(stack), coll(col_name)))==0){
    rast <- static_stack[[str_detect(names(static_stack), coll(col_name))]]
  } else{
  rast <- stack[[str_detect(names(stack), coll(col_name))]]
  }
  
  #extract values
  na_year[[i]][[col_name]] <- terra::extract(rast, na_year[[i]], fun=function(x){mean(x, na.rm=T)})[,2]
  #replace NaN values with NA
  na_year[[i]][[col_name]] <- ifelse(is.nan(na_year[[i]][[col_name]]), NA, na_year[[i]][[col_name]])
  }
}

#substitute interpolated values and recombine datafarme
na_only_sub <- bind_rows(na_year) %>% as.data.frame() %>% select(-geometry)
na_split$`TRUE` <- na_only_sub
out <- bind_rows(na_split) %>% arrange(t2_, desc(case_)) #animal_id

return(out)
}
}

#function to map interpolation function across all covariates
cov_interp_map <- function(steps_interp){
#Interpolate missing values at 300m radius ~13 min
for(i in 1:length(covs_to_interp_300m)){
  steps_interp <- interp_cov_vals(covs_to_interp_300m[i], df = steps_interp, buffer_size = 300)
print(paste0(i, "/", length(covs_to_interp_300m)))
}

#Interpolate missing values at 600m radius
for(i in 1:length(covs_to_interp_600m)){
  steps_interp<- interp_cov_vals(covs_to_interp_600m[i], df = steps_interp, buffer_size = 600)
  print(paste0(i, "/", length(covs_to_interp_600m)))
}

return(steps_interp)
}

steps_interp <- amt_steps_all_covs

#map interpolation across all indidivual step frames ~ 15 min
steps_interp$steps <- map(amt_steps_all_covs$steps, cov_interp_map, .progress = T)

#Check remaining NA points
steps_interp %>% unnest(cols=c(steps)) %>% select(-na_status) %>% DataExplorer::plot_missing()

#########################################################################
##
## 15. Round mean values for landuse and landcover to nearest integer
##
##########################################################################

steps_unscaled <- steps_interp %>% 
  mutate(land_cover_usfs_annual = round(land_cover_usfs_annual),
         land_use_usfs_annual = round(land_use_usfs_annual))


#########################################################################
##
## 16. Create final step data frame for export
##
##########################################################################
#add unique id and rearrange fields to final format
steps_final <- steps_unscaled %>% 
  mutate(unique_id = 1:nrow(steps_unscaled), .before= animal_id) %>% 
  select(-na_status)
  
#inspect points
steps_final %>% filter(is.na(popdens_hii_annual)) %>% 
  filter(case_==FALSE) %>% 
  st_as_sf(coords=c("x2_", "y2_"), crs=project_crs, remove=FALSE) %>% 
  mapview()

#last check for dupes
steps_final %>% 
  filter(case_==TRUE) %>% 
  get_dupes(animal_id, t1_)

steps_final %>% group_by(animal_id, step_id_) %>% count() %>% pull(n) %>% unique()

steps_final %>% group_by(animal_id, step_id_) %>% count() %>% filter(n!=11)

#export
#write_csv(steps_final, "data/Location_Data/Steps/2h_steps_unscaled_no_imp_annual_cov_02-13-2024.csv")


#########################################################################
##
## 17. Summary Table
##
#########################################################################

total_locations <- locs_screened %>% 
  group_by(sex, dispersal_status) %>% 
  count() %>% 
  unite("sex_disp", sex, dispersal_status) %>% 
  pivot_wider(names_from = sex_disp, values_from = n) %>% 
  mutate(Total = sum(c_across(everything()))) %>% 
  select(Female_resident, Female_disperser, Male_resident, Male_disperser, Total)

total_individuals <- locs_screened %>% 
  distinct(animal_id, .keep_all = T) %>% 
  group_by(sex, dispersal_status) %>% 
  count() %>% 
  unite("sex_disp", sex, dispersal_status) %>% 
  pivot_wider(names_from = sex_disp, values_from = n) %>% 
  mutate(Total = sum(c_across(everything()))) %>% 
  select(Female_resident, Female_disperser, Male_resident, Male_disperser, Total)

track_resample <- tracks_resampled %>% 
  select(-data) %>% 
  unnest(cols = tracks) %>% 
  group_by(sex, dispersal_status) %>% 
  count() %>% 
  unite("sex_disp", sex, dispersal_status) %>% 
  pivot_wider(names_from = sex_disp, values_from = n) %>% 
  mutate(Total = sum(c_across(everything()))) %>% 
  select(Female_resident, Female_disperser, Male_resident, Male_disperser, Total)

step_conversion <- amt_steps %>% 
  unnest(cols=steps) %>%  
  group_by(sex, dispersal_status) %>% 
  count() %>% 
  unite("sex_disp", sex, dispersal_status) %>% 
  pivot_wider(names_from = sex_disp, values_from = n) %>% 
  mutate(Total = sum(c_across(everything()))) %>% 
  select(Female_resident, Female_disperser, Male_resident, Male_disperser, Total)
  

step_filter_100m <- amt_steps %>% 
  unnest(cols = steps) %>% 
  filter(sl_ >= 100) %>% 
  group_by(sex, dispersal_status) %>% 
  count() %>% 
  unite("sex_disp", sex, dispersal_status) %>% 
  pivot_wider(names_from = sex_disp, values_from = n) %>% 
  mutate(Total = sum(c_across(everything()))) %>% 
  select(Female_resident, Female_disperser, Male_resident, Male_disperser, Total)

step_filter_20km <- amt_steps %>% 
  unnest(cols = steps) %>% 
  filter(sl_ >= 100) %>% 
  filter(sl_ <= 20000) %>% 
  group_by(sex, dispersal_status) %>% 
  count() %>% 
  unite("sex_disp", sex, dispersal_status) %>% 
  pivot_wider(names_from = sex_disp, values_from = n) %>% 
  mutate(Total = sum(c_across(everything()))) %>% 
  select(Female_resident, Female_disperser, Male_resident, Male_disperser, Total)

insuff_data <- sl_screen2 %>% 
  unnest(cols = steps) %>% 
  group_by(sex, dispersal_status) %>% 
  count() %>% 
  unite("sex_disp", sex, dispersal_status) %>% 
  pivot_wider(names_from = sex_disp, values_from = n) %>% 
  mutate(Total = sum(c_across(everything()))) %>% 
  select(Female_resident, Female_disperser, Male_resident, Male_disperser, Total)
  
insuff_data_rs <- amt_steps_r2 %>% 
  unnest(cols=steps) %>% 
  filter(case_==TRUE) %>% 
  group_by(sex, dispersal_status) %>% 
  count() %>% 
  unite("sex_disp", sex, dispersal_status) %>% 
  pivot_wider(names_from = sex_disp, values_from = n) %>% 
  mutate(Total = sum(c_across(everything()))) %>% 
  select(Female_resident, Female_disperser, Male_resident, Male_disperser, Total)

remaining_individuals <- amt_steps_r2 %>% 
  unnest(cols = steps) %>% 
  distinct(animal_id, .keep_all = T) %>% 
  group_by(sex, dispersal_status) %>% 
  count() %>% 
  unite("sex_disp", sex, dispersal_status) %>% 
  pivot_wider(names_from = sex_disp, values_from = n) %>% 
  mutate(Total = sum(c_across(everything()))) %>% 
  select(Female_resident, Female_disperser, Male_resident, Male_disperser, Total)


tab_dat <- bind_rows(total_individuals,
                     total_locations,
                     track_resample,
                     step_conversion,
                     step_filter_100m,
                     step_filter_20km,
                     insuff_data,
                     insuff_data_rs,
                     remaining_individuals)

row_names <- c("Total individuals", 
               "Total locations",
                "N locations after filtering to 2-hour intervals (+/- 45 min)",
                "N steps after converting tracks",
                "N steps after omitting step lengths < 100m apart", 
                "N steps after omitting step lengths > 20km apart", 
                "N steps after excluding individuals with < 30 steps or < 30 days of steps", 
                "N steps after excluding individuals with insufficient data to fit movement distributions", 
                "Total individuals for analysis"
               )

summary_table <- tab_dat %>% 
  mutate(`Processing Step` = row_names, .before = Female_resident)


table <- flextable(summary_table) %>% 
  separate_header %>% 
  width(j = "Processing Step", width = 2.5) %>% 
  align(align = "center", part = "all") %>% 
  align(j = "Processing Step", align = "left")  

#save_as_image(table, "step_processing_table_training_dat_2-13-24.png")
#save_as_docx("table1" = table, path = "step_processing_table_training_dat_2-13-24.docx")



#make same table with differences instead of counts
diff_tab <- summary_table %>% 
  filter(!(`Processing Step` %in% c("Total individuals", "Total individuals for analysis"))) %>% 
  select(-`Processing Step`) %>% 
  mutate_all(~ . - lag(.)) %>% 
  filter(!is.na(Total))

diff_tab_full <- bind_rows(total_individuals,
                           total_locations,
                           diff_tab,
                           insuff_data_rs,
                           remaining_individuals
                           )


row_names_diff <- c("Total individuals", 
                    "Total locations",
                    "Filter locations to 2-hour intervals (+/- 45 min)",
                    "Convert tracks to steps",
                    "Omit step lengths < 100m apart", 
                    "Omit step lengths > 20km apart", 
                    "Exclude individuals with < 30 steps or < 30 days of steps", 
                    "Exclude individuals with insufficient data to fit movement distributions",
                    "Total steps for analysis",
                    "Total individuals for analysis"
)

summary_table_diff <- diff_tab_full %>% 
  mutate(`Processing Step` = row_names_diff, .before = Female_resident)
  
diff_table <- flextable(summary_table_diff) %>% 
  separate_header %>% 
  width(j = "Processing Step", width = 2.5) %>% 
  align(align = "center", part = "all") %>% 
  align(j = "Processing Step", align = "left")  

#save_as_image(diff_table, "step_processing_diff_table_training_dat_2-13-24.png")
#save_as_docx("table1" = diff_table, path = "step_processing_diff_table_training_dat_2-13-24.docx")



