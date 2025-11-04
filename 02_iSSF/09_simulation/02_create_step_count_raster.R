#### Create and step count raster ####

# Author: Read Barbee

# Date:2025-10-31 

# Purpose: Create and validate step count raster
#Inputs:
#1. Simulated points

#Outputs: 
#1. Raster of simulated step counts per cell


################################ Libraries #################################
library(tidyverse)
library(sf)
library(terra)
library(data.table)
library(future.apply)


point_path <- "/Volumes/Reads_Seagate/ms_thesis/ssf_sims"

point_files <- list.files(point_path, pattern = ".csv", full.names = T)

#scaled prediction raster
cov_stack_pred <- rast("data/spatial/raster/pred_stack_2023_top_ssf_train_scaled_uncapped_9-18-25.tif")

#study area polygon for masking
op_poly <- st_read("data/spatial/vector/ocp_study_area_poly_wa_only_10-30-23.gpkg")


########################################################################
##
## 1. Tally number of steps taken per raster cell
##
##########################################################################

## template rasters for count calculation #
tmp_rast30 <- cov_stack_pred[[1]] #vals
values(tmp_rast30) <- 0#  initialize with zeros


# Precompute total number of cells
n_cells <- ncell(tmp_rast30)

# Initialize an empty numeric vector to store counts
total_counts <- numeric(n_cells)

# Loop through each point file and accumulate counts
for (i in seq_along(point_files)) {
  pts <- fread(point_files[i], select = c("x", "y"))
  
  # Convert point coordinates to cell indices
  cells <- cellFromXY(tmp_rast30, pts)
  
  # Count occurrences per cell
  counts <- tabulate(cells, nbins = n_cells)
  
  # Accumulate counts
  total_counts <- total_counts + counts
  
  cat("Processed", i, "of", length(point_files), "files\n")
}

# Assign final counts to raster
values(tmp_rast30) <- total_counts

#set ocean values as NA
count_rast_ocean_masked <- mask(tmp_rast30, op_poly)

# Write out the combined count raster
writeRaster(count_rast_ocean_masked,
            "results/comparison_layers/ssf/simulation/ssf_sim_paths_raw_counts_200000p_4380s_both_sides_I5_30m_uncapped_10-31-25.tif",
            overwrite = TRUE)


########################################################################
##
## 2. Create binned raster
##
##########################################################################

#extract raster values
pred_vals <- terra::values(count_rast_ocean_masked)

#identify quantile breaks
breaks <- quantile(pred_vals, probs = 0:10/10, na.rm = T) #obtain 10 quantile 

#classify raster by breaks
binned <- terra::classify(count_rast_ocean_masked, rcl=as.vector(breaks), include.lowest=TRUE)

#set factor levels
levels(binned) <- 1:10

#writeRaster(binned,  "results/comparison_layers/ssf/simulation/ssf_sim_paths_binned_counts_200000p_4380s_both_sides_I5_30m_uncapped_10-31-25.tif", overwrite=T)


