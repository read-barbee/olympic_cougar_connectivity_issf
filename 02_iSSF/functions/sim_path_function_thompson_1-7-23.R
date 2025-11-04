# Simulate random paths based on SSF input
#
# n_lists: integer > 0; number of independent lists of paths to make (important for parallelization as each list will be done in parallel, unless specified otherwise using the "parallel" argument)
# n_paths_per_list: integer > 0; number of paths per list (these will all be included in the same data.frame with a path_id column)
# n_steps_per_path: integer > 0; how long is each path?
# n_rand: integer > 0; when picking random points to choose from (including initial locations), how many do we pick? Higher values will have more computational cost but will be more accurate
# move_pars: named numeric vector; contains parameters necessary for simulating random step lengths and turning angles, names include "gamma_shape", "gamma_rate", "vm_mu", "vm_kappa"
# R_rsf: SpatRaster object; contains environmental data (i.e., the untransformed output from an RSF or something similar)
# R_step: SpatRaster object; contains environmental data that are meant to interact with the step length (i.e., rate parameter of gamma distribution). Note that this should include the "intercept" (i.e., the "step-only" term) and thus the intercept can't be supplied elsewhere
# R_log_step: SpatRaster object; contains environmental data that are meant to interact with the log of the step length (i.e., shape parameter of gamma distribution). Note that this should include the "intercept" (i.e., the "log-step-only" term) and thus the intercept can't be supplied elsewhere
# R_cos_angle: SpatRaster object; contains environmental data that are meant to interact with the cosine of the turning angles (i.e., concentration parameter of the von Mises distribution). Note that this should include the "intercept" (i.e., the "angle-only" term) and thus the intercept can't be supplied elsewhere
# R_daylight: SpatRaster object; contains linear combination of environmental covariates that are interacting with daylight_var
# R_mask: SpatRaster object; should be 0 or 1 (or NA) everywhere; determines whether certain areas are "in bounds" for steps
# step_cos_angle_par: numeric; interaction parameter between step length ("rate" parameter of gamma distribution) and cosine of turning angle ("kappa" parameter of von Mises distribution)
# log_step_cos_angle_par: numeric; interaction parameter between log(step length) ("shape" parameter of gamma distribution) and cosine of turning angle ("kappa" parameter of von Mises distribution)
# step_day_par: numeric; interaction parameter between step length ("rate" parameter of gamma distribution) and daylight_var
# log_step_day_par: numeric; interaction parameter between log(step length) ("shape" parameter of gamma distribution) and daylight_var
# cos_angle_day_par: numeric; interaction parameter between cos(turn angle) ("kappa" parameter of von Mises distribution) and daylight_var
# x0_bounds: SpatVector polygon object; space from which to sample initial locations
# daylight_var: numeric vector with length equal to n_steps_per_path; contains information on some quantity (e.g., time of day) at each time step
# season_var: integer vector; all values should be between 1 and the # of layers in the covariate rasters ("R_..."). Tells us which index to grab from at each point in time.
# n_check_mask: integer > 0; how many points along a step do we check for being out of bounds?
# step_factor: numeric > 0; how much do we multiply randomly generated steps by to get them to their true values? Can be relevant if move_pars are on a different scale (e.g., metres versus kilometres)
# n_cores: integer > 0; how many cores do we divide the work into? Should be at most n_lists. If it's less than or equal to 1 we don't use parallelization
# n_print: integer > 0; how often do we print an update on the progress of the function? If it's greater than n_steps_per_path then we never print anything
# ret_list: logical; only really relevant when n_lists == 1, and in that case, tells the function whether to return a list (with 1 element) or a data.frame
#
# Returns a data.frame (or list of data.frames) for paths
sim_paths = function(n_lists = 1, 
                     n_paths_per_list = 1, 
                     n_steps_per_path = 100, 
                     n_rand = 10,
                     move_pars = NULL,
                     R_rsf = NULL, 
                     R_step = NULL,
                     R_log_step = NULL,
                     R_cos_angle = NULL,
                     R_daylight = NULL,
                     R_mask = NULL,
                     step_cos_angle_par = 0,
                     log_step_cos_angle_par = 0,
                     step_day_par = 0,
                     log_step_day_par = 0,
                     cos_angle_day_par = 0,
                     x0_bounds = NULL, # TO DO: Could make an option for this to be a list?
                     daylight_var = numeric(n_steps_per_path),
                     season_var = rep(1, n_steps_per_path),
                     n_check_mask = 1,
                     step_factor = 1,
                     n_cores = n_lists, 
                     n_print = n_steps_per_path + 1,
                     ret_list = (n_lists > 1)) {
  
  require(doParallel)
  require(terra)
  require(circular)
  require(tidyverse)
  
  if (n_cores > 1) registerDoParallel(cores = min(n_cores, n_lists))
  
  if (is.null(R_rsf) & is.null(R_mask)) {
    # Make a blank raster with the default resolution
    R_rsf = rast(vals = 0, crs = "+init:epsg=5070")
    # vals = 1 here because we want the animal to be able to go everywhere
    R_mask = rast(vals = 1, crs = "+init:epsg=5070")
  } else if (is.null(R_rsf)) {
    # Make a new raster using the template from R_mask
    R_rsf = rast(vals = 0, extent = ext(R_mask), resolution = res(R_mask), crs = crs(R_mask))
  } else if (is.null(R_mask)) {
    # Make a mask that's 1 everywhere (i.e., everywhere's good)
    R_mask = rast(vals = 1, extent = ext(R_rsf), resolution = res(R_rsf), crs = crs(R_rsf))
  }
  
  if (is.null(x0_bounds)) x0_bounds = vect(ext(R_mask), crs = crs(R_mask))
  if (is.null(move_pars)) move_pars = c(gamma_shape = 1, gamma_rate = 1, vm_mu = 0, vm_kappa = 0)
  
  i_dist_prop = seq(1/n_check_mask, 1, length.out = n_check_mask) # get proportions to check mask
  
  if (n_cores > 1) {
    `%dothis%` = `%dopar%`
  } else {
    `%dothis%` = `%do%`
  }
  
  list_of_all_paths = foreach(yyy = 1:n_lists) %dothis% {
    
    # Random start locations within boundary
    xy = spatSample(x0_bounds, size = n_rand * n_paths_per_list)
    xy = crds(xy)
    
    # Extract values from the raw RSF raster to pick the random initial locations in a way that represents habitat quality
    t_ssf <- terra::extract(R_rsf[[season_var[1]]], xy)[ ,1]
    t_zero <- terra::extract(R_mask, xy)[ ,1]
    
    ud <- exp(t_ssf) * t_zero
    # convert NA values to 0
    ud[!is.finite(ud)] = 0
    keep <- sample(1:nrow(xy), size = n_paths_per_list, replace = FALSE, prob = ud) 
    xy_st <- xy[keep, ]
    
    tmp0 <- tibble(path_id = 1:n_paths_per_list, i_strata = 1, x_prev = as.numeric(xy_st[, 1]), y_prev = as.numeric(xy_st[, 2]), 
                   x = NA, y = NA, t_var = daylight_var[1], direction_prev = round(runif(n_paths_per_list, -pi, pi), 3))
    
    path_list = vector(mode = "list", length = n_steps_per_path + 1)
    
    for (i in 1:n_steps_per_path) {
      
      tmp1 <- tmp0 %>% mutate(i_strata = i, t_var = daylight_var[i])
      
      path_nums = rep(tmp1$path_id, each = n_rand)
      r_steplengths = rgamma(length(path_nums), move_pars["gamma_shape"], move_pars["gamma_rate"])
      r_turn_angles = as.numeric(rvonmises(length(path_nums), circular(move_pars["vm_mu"]), move_pars["vm_kappa"]))
      move_metrics_df = data.frame(step = r_steplengths * step_factor, angle = r_turn_angles, path_id = path_nums)
      
      tmp2a = inner_join(tmp1, move_metrics_df, by = "path_id", multiple = "all") %>% 
        mutate(direction = (direction_prev + angle) %% (2*pi)) %>% 
        mutate(x = x_prev + step * cos(direction), y = y_prev + step * sin(direction))
      
      # Filter steps that cross inhospitible habitat - check midway points for fast steps as they are longer
      zero_mx <- matrix(NA, nrow = nrow(tmp2a), ncol = length(i_dist_prop))
      for (tt in 1:length(i_dist_prop)){
        tt_x <- i_dist_prop[tt] * tmp2a$step * cos(tmp2a$direction) + tmp2a$x_prev
        tt_y <- i_dist_prop[tt] * tmp2a$step * sin(tmp2a$direction) + tmp2a$y_prev
        zero_mx[ ,tt] <- terra::extract(R_mask, cbind(tt_x, tt_y))[,1]
      }
      zero_mx[!is.finite(zero_mx)] = 0
      tmp2b <- tmp2a %>% mutate(step_good = ifelse(rowSums(zero_mx) == length(i_dist_prop), 1, 0))
      # Identify whether we need to "re-roll" any paths because none of the steps ended up inside the mask
      tmp2b = tmp2b %>% group_by(path_id) %>% mutate(reroll = sum(step_good)) %>% ungroup
      
      while(any(tmp2b$reroll == 0)) {
        
        tmp2b_good = tmp2b %>% dplyr::filter(reroll > 0)
        tmp2b_bad = tmp2b %>% dplyr::filter(reroll == 0)
        
        # for debugging purposes
        # message("Entering while loop for re-rolling at time ", i, ". Have to re-roll ", nrow(tmp2b_bad) / n_rand, " paths.")
        
        paths_bad = tmp2b_bad$path_id
        r_steplengths_bad = rgamma(length(paths_bad), move_pars["gamma_shape"], move_pars["gamma_rate"])
        r_turn_angles_bad = as.numeric(rvonmises(length(paths_bad), circular(move_pars["vm_mu"]), move_pars["vm_kappa"]))
        move_metrics_df_bad = data.frame(step = r_steplengths_bad * step_factor, angle = r_turn_angles_bad, path_id = paths_bad)
        
        paths_all_bad = unique(paths_bad) # get all paths that need to be re-rolled
        
        tmp2a_bad = inner_join(tmp1[tmp1$path_id %in% paths_all_bad, ], move_metrics_df_bad, by = "path_id", multiple = "all") %>% 
          mutate(direction = (direction_prev + angle) %% (2*pi)) %>% 
          mutate(x = x_prev + step * cos(direction), y = y_prev + step * sin(direction))
        
        zero_mx <- matrix(NA, nrow = nrow(tmp2a_bad), ncol = length(i_dist_prop))
        for (tt in 1:length(i_dist_prop)){
          tt_x <- i_dist_prop[tt] * tmp2a_bad$step * cos(tmp2a_bad$direction) + tmp2a_bad$x_prev
          tt_y <- i_dist_prop[tt] * tmp2a_bad$step * sin(tmp2a_bad$direction) + tmp2a_bad$y_prev
          zero_mx[ ,tt] <- terra::extract(R_mask, cbind(tt_x, tt_y))[,1]
        }
        zero_mx[!is.finite(zero_mx)] = 0
        tmp2b_bad_new <- tmp2a_bad %>% mutate(step_good = ifelse(rowSums(zero_mx) == length(i_dist_prop), 1, 0)) %>%
          group_by(path_id) %>% mutate(reroll = sum(step_good)) %>% ungroup
        
        tmp2b = rbind(tmp2b_good, tmp2b_bad_new) %>% arrange(path_id)
        
      }
      
      # fix steps of length 0 because they'll give problems later on in the code (we take the log)
      tmp2b = tmp2b %>% mutate(step = ifelse(step == 0, 1e-5, step))
      
      # Get covariate values at all proposed steps
      xy_vals_extract = cbind(tmp2b$x, tmp2b$y)
      
      log_p = terra::extract(R_rsf[[season_var[i]]], xy_vals_extract)[ ,1]
      if (!is.null(R_daylight)) log_p = log_p + terra::extract(R_daylight[[season_var[i]]], xy_vals_extract)[, 1] * daylight_var[i]
      if (!is.null(R_step)) log_p = log_p + terra::extract(R_step[[season_var[i]]], xy_vals_extract)[, 1] * tmp2b$step
      if (!is.null(R_log_step)) log_p = log_p + terra::extract(R_log_step[[season_var[i]]], xy_vals_extract)[, 1] * log(tmp2b$step)
      if (!is.null(R_cos_angle)) log_p = log_p + terra::extract(R_cos_angle[[season_var[i]]], xy_vals_extract)[, 1] * cos(tmp2b$angle)
      
      # make these coefficients seasonal if necessary
      step_day_par_i = ifelse(length(step_day_par) > 1, step_day_par[[season_var[i]]], step_day_par)
      log_step_day_par_i = ifelse(length(log_step_day_par) > 1, log_step_day_par[[season_var[i]]], log_step_day_par)
      cos_angle_day_par_i = ifelse(length(cos_angle_day_par) > 1, cos_angle_day_par[[season_var[i]]], cos_angle_day_par)
      log_step_cos_angle_par_i = ifelse(length(log_step_cos_angle_par) > 1, log_step_cos_angle_par[[season_var[i]]], log_step_cos_angle_par)
      step_cos_angle_par_i = ifelse(length(step_cos_angle_par) > 1, step_cos_angle_par[[season_var[i]]], step_cos_angle_par)
      
      tmp2b <- tmp2b %>% mutate(logp = log_p + daylight_var[i] * 
                                  (step_day_par_i * step + log_step_day_par_i * log(step) + cos_angle_day_par_i * cos(angle)) + 
                                  log(step) * cos(angle) * log_step_cos_angle_par_i + 
                                  step * cos(angle) * step_cos_angle_par_i)
      
      # Add in step effects
      tmp2c = tmp2b %>% mutate(exp_lp = exp(logp) * step_good)
      tmp2c$exp_lp[!is.finite(tmp2c$exp_lp)] = 0
      tmp2c = tmp2c %>% group_by(path_id) %>% mutate(sum_exp_lp = sum(exp_lp)) %>% ungroup
      # In the event that the probabliities all round off and become 0 (seems quite rare) we just pick from them equally (they are all equally)
      tmp2c = tmp2c %>% mutate(exp_lp = ifelse(sum_exp_lp == 0, step_good, exp_lp))
      # Pick a proposed step based on exp_lp
      tmp3b <- tmp2c %>% group_by(path_id) %>% slice_sample(n = 1, weight_by = exp_lp) %>% ungroup
      
      tmp_save <- tmp3b %>% dplyr::select(path_id, i_strata, x, y, t_var, step, direction)
      path_list[[i]] = tmp_save
      
      # for use in the next iteration
      tmp4 <- tmp3b %>% 
        mutate(x_prev = x, y_prev = y, x = NA, y = NA, direction_prev = direction) %>% 
        dplyr::select(all_of(names(tmp0)))
      tmp0 <- tmp4   # use tmp0 in next iteration
      
      if (i %% n_print == 0) message("Completed step ", i, " of ", n_steps_per_path)
      
    }  # End of i (individual steps)
    
    path_all <- bind_rows(path_list)
    path_all %>% mutate(i_rep = yyy, .before = path_id)
    
  }
  
  if (n_lists == 1 & !ret_list) return(list_of_all_paths[[1]]) # just give a data.frame
  list_of_all_paths
  
}