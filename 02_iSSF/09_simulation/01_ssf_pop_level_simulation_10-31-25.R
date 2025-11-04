#### iSSF simulation ####

# Author: Read Barbee

# Date:2025-10-31

# Purpose:

#1. Extract fixed effects coefficients from top SSF model
#2. Create raster surface by applying model to landscape
#3. Simulate tracks using population-level fixed effects



################################ Libraries #################################
library(tidyverse)
library(terra)
library(sf)
library(INLA)
library(flextable)

#simulation function
source("r_scripts/02_iSSF/07_simulation/sim_path_function_thompson_1-7-23.R")


################################ Import data #################################

#model data for demography and individual IDs
mod_dat <- read_csv("data/fitted_models/top_ssf_inla_corrected_stratum_data_2-16-24.rds")

#extract individual IDS and demography
animal_info <- mod_dat %>% distinct(animal_id, .keep_all = TRUE) %>% 
  select(id1, animal_id:dispersal_status) %>% 
  rename(ID = id1)

#scaling parameters from training data to use for prediction surface (only once)
scaling_parameters <- read_csv("data/fitted_models/scaling_parameters_top_ssf_mod_2-16-24.csv")

#scaled prediction raster
cov_stack_pred <- rast("data/spatial/raster/pred_stack_2023_top_ssf_train_scaled_uncapped_9-18-25.tif")

#step data for movement parameters 
steps <- read_csv("data/gps_data/05_steps/2h_steps_unscaled_no_imp_annual_cov_02-13-2024.csv")

###############################################################################

#1. Scale prediction raster using scaling parameters from training data (only once)

###############################################################################

# cov_stack_pred_raw <- rast("/Volumes/Reads_Seagate/ms_thesis/prediction_rasters/pred_stack_2023_raw_uncapped.tif")

# pred_stack_mod <- cov_stack_pred_raw[[names(cov_stack_pred_raw) %in% scaling_parameters$covariate]]
# 
# scaled_rasts <- list()
# 
# for(i in 1:nlyr(pred_stack_mod)){
#   cov <- names(pred_stack_mod)[i]
#   
#   center_p <- scaling_parameters %>% 
#     filter(covariate == !!cov) %>% 
#     pull(center)
#   
#   scale_p <- scaling_parameters %>% 
#     filter(covariate == !!cov) %>% 
#     pull(scale)
#   
#   scaled_rasts[[i]] <- (pred_stack_mod[[i]] - center_p)/scale_p
#   
#   print(paste0(i, "/", nlyr(pred_stack_mod)))
# }
# 
# scaled_stack <- rast(scaled_rasts)
# names(scaled_stack) <- names(pred_stack_mod)
# 
# scaled_stack$`I(perc_tree_cov_annual^2)` <- scaled_stack$perc_tree_cov_annual^2
# 

#writeRaster(scaled_stack, "/Volumes/Reads_Seagate/ms_thesis/prediction_rasters/pred_stack_2023_top_ssf_train_scaled_uncapped_9-18-25.tif")

rm(cov_stack_pred_raw)
rm(pred_stack_mod)
rm(scaled_rasts)


#########################################################################
##
## 2.Import or fit top global model
##
##########################################################################

#current top model 
top_mod <- readRDS("data/fitted_models/top_ssf_inla_corrected_stratum_2-16-24.rds")

coeff_names <-  top_mod$names.fixed

coeff_names <- coeff_names[!(coeff_names %in% c("sl_", "log_sl_", "cos_ta_"))]

#########################################################################
##
## 3.Extract population-level coefficients and sl and ta distributions
##
##########################################################################
pop_coeffs <- top_mod$summary.fixed %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  filter(!(rowname %in% c("sl_", "log_sl_", "cos_ta_"))) %>% 
  select(rowname, mean)

used <- steps %>% filter(case_== TRUE)

sl_dist <- amt::fit_distr(used$sl_, "gamma")

ta_dist <- amt::fit_distr(used$ta_, "vonmises")

sl_shape <- sl_dist$params$shape
sl_scale <- sl_dist$params$scale

sl_rate <- 1/sl_scale #gamma rate parameter is 1/scale

ta_kappa <- ta_dist$params$kappa
ta_mu <- ta_dist$params$mu

#########################################################################
##
## 4.Make population level map
##
##########################################################################
names_linear_no_movement <- pop_coeffs %>% 
  #filter(rowname!= "I(perc_tree_cov_annual^2)") %>% 
  pull(rowname)

cov_stack_sel <- cov_stack_pred[[names_linear_no_movement]]

#Generate predictions by multiplying each covariate layer by its coefficient estimate from the model
pred <- list()
for(j in 1:nrow(pop_coeffs)){
  
  cf <- pop_coeffs$mean[j]
  
  names <- pop_coeffs$rowname
  
    pred[[j]] <- cov_stack_sel[[names[j]]] * cf #breaks when mean coeff is 0
    print(paste0(j, "/", nrow(pop_coeffs)))
}

predictions <- Reduce("+", pred) #sum layers together
pred_vals <- terra::values(predictions) #convert raster to values for calculations
preds_exp <- plogis(pred_vals) #backtransform predictions from logit to prob scale

#Clamp extreme values (Sells et al. 2022)
x.min <- quantile(preds_exp, 0.025, na.rm = T) 
x.max <- quantile(preds_exp, 0.975, na.rm = T)
preds_exp[preds_exp > x.max] <- x.max
preds_exp[preds_exp < x.min] <- x.min
x.range <- x.max - x.min
predictions. <- (preds_exp - x.min)/x.range

vals <- cov_stack_sel[[1]] %>% terra::setValues(predictions.) #project those values onto the map

#export
#writeRaster(vals, "/Users/tb201494/Desktop/Panthera/pumas/projects/olympic_cougar_connectivity/results/comparison_layers/ssf/probability_of_use/pop_level_map_uncapped_train_scaled_normalized_11-04-25.tif")

#vals <- rast("predictive surfaces/ssf/predictions/ssf_projected_prob_use_standardized_max_cap_30m_3-08-24.tif")

#########################################################################
##
## 5.Set up simulation
##
##########################################################################

move_pars <- c("gamma_shape" = sl_shape, "gamma_rate" = sl_rate, "vm_mu" = ta_mu, "vm_kappa" = ta_kappa)

R_rsf <-  vals #predictions 

r_mask <- rast("data/spatial/op_water_mask2_01-07-24.tif")

start_zone2 <- vect("data/spatial/ocp_study_area_poly_wa_only_10-30-23.gpkg")

#########################################################################
##
## 6.Simulate
##
##########################################################################


# ~ 50 min for 2 paths of 5000 steps on one core
# ~ 4 min for 100 paths of 100 steps over 4 cores
#~ 3 hours for 100 paths of 5000 steps over 5 cores
#  5.32 hours for 280 paths of 5000 steps over 7 cores
#1.6 hours for 9000 paths of 5000 points over 6 cores?

for(i in 1:20){
  
indivs = 8
paths_per_indiv = 1250
steps_per_path = 4380

system.time(paths <- sim_paths(
         n_lists = indivs, 
         n_paths_per_list = paths_per_indiv, 
         n_steps_per_path = steps_per_path, 
         n_rand = 10,
         move_pars = move_pars,
         R_rsf = R_rsf, 
         R_step = NULL,
         R_log_step = NULL,
         R_cos_angle = NULL,
         R_daylight = NULL,
         R_mask = r_mask, 
         step_cos_angle_par = 0,
         log_step_cos_angle_par = 0,
         step_day_par = 0,
         log_step_day_par = 0,
         cos_angle_day_par = 0,
         x0_bounds = start_zone2, 
         daylight_var = numeric(steps_per_path),
         season_var = rep(1, steps_per_path),
         n_check_mask = 1,
         step_factor = 1,
         n_cores = indivs, 
         n_print = 1000,
         ret_list = (indivs > 1)) )

#visualise a sample of simulated paths
#map_test <- paths %>% bind_rows() %>% unite("unique_path", i_rep, path_id) %>% st_as_sf(coords = c("x", "y"), crs = 5070) 

#mapview::mapview(map_test, zcol = "unique_path") + mapview::mapview(raster::raster(vals))

s_paths <- bind_rows(paths) %>% unite("unique_path", i_rep, path_id) %>% arrange(unique_path)

path_ids <- s_paths %>% distinct(unique_path) %>% 
  mutate(path_id = 1:nrow(.))

paths_out <- s_paths %>% 
  left_join(path_ids, join_by(unique_path)) %>% 
  select(path_id, i_strata:direction) %>% 
  rename(step_id = i_strata)

point_name <- paste0("pop_sim_paths_10000p_4380s_both_sides_I5_round_", i, "_10-06-25.csv")

write_csv(paths_out, point_name)

print(i)
}
