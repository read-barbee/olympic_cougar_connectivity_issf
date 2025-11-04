#### Compare covariate distributions for original and imputed steps ####

# Author: Read Barbee

# Date:2025-10-27 

# Purpose: Extract covariate values and plot distributions of known and imputed locations for each covariate value. Differences in distributions may indicate a habitat bias in GPS fixes. 

#Inputs:
#1. Imputed GPS tracks

#Outputs:
#1. Used/imputed histogram plots for each covariate.


################################ Libraries #################################
library(tidyverse)
library(sf)
library(terra)

#########################################################################
##
## 1. Import Location Data
##
##########################################################################

sims_rerouted <- read_csv("data/gps_data/imputed_paths/ctmm_sim_imputed_paths_rerouted_12-04-2023.csv")

imp_rerouted_sf <- sims_rerouted %>% st_as_sf(coords=c("x_rerouted", "y_rerouted"), crs = 5070)


#########################################################################
##
## 2. Import Covariates
##
##########################################################################
cov_folder_path <- "/Volumes/Reads_Seagate/ms_thesis/covariates/annual_covs_01-07-24/30m"

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
## 2. Extract covariates by year
##
##########################################################################

extract_annual_covs <- function(points, cov_stack_list){
  #create year column
  points2 <- points %>%
    mutate(year = as.factor(year(timestamp)))
  
  #split data by year
  points_split <- split(points2, points2$year)
  
  #define function to select the raster stack for the relevant year to extract from and extract the covariates
  extract_fun <- function(points_split, year){
    stack <- cov_stack_list[[paste0("cov_stack_", as.character(year))]] # "_water_masked"
    covs <- terra::extract(stack, points_split, bind = TRUE) %>% as_tibble()
    
    return(covs)
  }
  
  #apply the extraction function to each year
  points_split_covs <- list()
  for(i in 1:length(points_split)){
    names <- names(points_split)
    points_split_covs[[i]] <- extract_fun(points_split[[i]], names[i])
  }
  
  names(points_split_covs) <- names
  
  remove_year_names <- function(points_split_covs){
    old_names <- names(points_split_covs)
    new_names <- vector()
    for(i in 1:length(old_names)){
      if(str_detect(old_names[i], coll("20"))){
        new_names[i] <- substr(old_names[i], 1, nchar(old_names[i]) - 5)
      } else{
        new_names[i] <- old_names[i]
      }
    }
    
    names(points_split_covs) <- new_names
    
    return(points_split_covs)
  }
  
  #apply renaming function
  covs_renamed <- map(points_split_covs, remove_year_names)
  
  #bind all years together 
  point_covs_final <- bind_rows(covs_renamed)
  
  return(point_covs_final)
  
}

#extract covariates at all locations
static_vals <- terra::extract(static_stack, imp_rerouted_sf, bind = TRUE) %>% as_tibble()

annual_vals <- extract_annual_covs(imp_rerouted_sf, cov_stacks)

all_covs <- annual_vals %>% 
  left_join(static_vals, by = c("fid", "animal_id", "timestamp", "x_old", "y_old", "vx", "vy", "imp_status")) %>% 
  relocate(elevation:distance_water, .before = ndvi_annual)


#########################################################################
##
## 3. Prepare dataframe for plotting
##
##########################################################################

vals_pivot <- all_covs %>% 
  select(!contains("usfs"), -aspect, -dist_very_small_roads_annual) %>% 
  pivot_longer(c(elevation:dens_all_roads_annual), names_to = "cov", values_to= "value") %>% 
  mutate(imp_status = as.factor(imp_status)) %>% 
  mutate(imp_status =fct_relevel(imp_status, "observed", "imputed")) %>% 
  mutate(cov = fct_recode(cov,
                          "Aspect (Eastness)" = "aspect_eastness",
                          "Aspect (Northness)" = "aspect_northness",
                          "Density All Roads" = "dens_all_roads_annual",
                          "Density Paved Roads" = "dens_paved_roads_annual",
                          "Density Unpaved Roads" = "dens_unpaved_roads_annual",
                          "Distance to All Roads" = "dist_all_roads_annual",
                          "Distance to Paved Roads" = "dist_paved_roads_annual",
                          "Distance to Unpaved Roads" = "dist_unpaved_roads_annual",
                          "Distance to Major Roads" = "dist_major_roads_annual",
                          "Distance to Minor Roads" = "dist_minor_roads_annual",
                          #"Distance to Very Small Roads" = "dist_very_small_roads_annual",
                          "Distance to Non-Motorized Roads" = "dist_non_motorized_roads_annual",
                          "Distance to Water" = "distance_water",
                          "Elevation" = "elevation",
                          "Slope" = "slope",
                          "mTPI" = "mtpi",
                          "TPI" = "tpi",
                          "TRI" = "tri",
                          "EVI" = "evi_annual",
                          "NDVI" = "ndvi_annual",
                          "GPP" = "gpp_annual",
                          "NPP" = "npp_annual",
                          "Precipitation" = "precip_annual",
                          "% Tree Cover" = "perc_tree_cov_annual",
                          "% Non-tree Veg" = "perc_nontree_veg_annual",
                          "% Non-vegetated" = "perc_nonveg_annual",
                          "Human Impact Index" = "hii_annual",
                          "Human Landuse Index (HII)" = "landuse_hii_annual",
                          "Human Infrastructure Index (HII)" = "infra_hii_annual",
                          "Human Population Density (HII)" = "popdens_hii_annual",
                          "Road Impact (HII)" = "roads_hii_annual"
                          )) %>% 
  mutate(cov =fct_relevel(cov,
                          "Elevation",
                          "Slope",
                          "Aspect (Northness)",
                          "Aspect (Eastness)",
                          "TRI",
                          "TPI",
                          "mTPI",
                          "NDVI",
                          "EVI",
                          "GPP",
                          "NPP",
                          "Precipitation",
                          "% Tree Cover",
                          "% Non-tree Veg",
                          "% Non-vegetated" ,
                          "Distance to Water",
                          "Human Impact Index",
                          "Human Landuse Index (HII)",
                          "Human Infrastructure Index (HII)",
                          "Human Population Density (HII)",
                          "Road Impact (HII)" ,
                          "Density All Roads" ,
                          "Density Paved Roads" ,
                          "Density Unpaved Roads" ,
                          "Distance to All Roads" ,
                          "Distance to Paved Roads",
                          "Distance to Unpaved Roads",
                          "Distance to Major Roads" ,
                          "Distance to Minor Roads",
                         #"Distance to Very Small Roads" ,
                          "Distance to Non-Motorized Roads"
                        
  )) %>%
  arrange(cov)

#density
ggplot(vals_pivot, aes(x = value, y = ..density.., fill = as.factor(imp_status)))+ geom_histogram(position="identity", alpha=0.7) + 
  facet_wrap(vars(cov), scales = "free") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) + 
  #theme(axis.title.x=element_text(size=16)) +
  labs(fill = "Location Type", y = "Density", x = "Observed vs. Imputed Distributions (All Years)")

#ggsave("observed_imputed_distributions_annual_covs_3-13-24.png", width = 25, height = 15, units = "in") 


#frequency
 ggplot(vals_pivot, aes(x = value, fill = as.factor(imp_status)))+ 
   geom_histogram(position="identity", alpha=0.7) + 
  facet_wrap(vars(cov), scales = "free") +
  xlab("Observed vs Imputed Distributions") + theme(axis.title.x=element_text(size=16)) 
 
 #map example of an imputed path
 imp_rerouted_sf %>% 
   filter(animal_id=="Al") %>% 
 mapview::mapview(zcol = "imp_status")
 
 
