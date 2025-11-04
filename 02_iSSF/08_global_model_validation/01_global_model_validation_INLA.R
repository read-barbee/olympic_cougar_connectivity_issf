#### Global Model Selection and Validation ####

# Author: Read Barbee

# Date:2023-09-18
#Last Updated: 2024-02-14

# Purpose:


################################ libraries #################################
library(tidyverse)
library(beepr)
library(INLA)
library(terra)
library(sf)

#source k-folds function
source("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/r_scripts/02_iSSF/08_cross_validation/k_fold_cv_function_INLA_v4_quad_02-06-24.R")


#########################################################################
##
##  Specify Model Parameters
##
##########################################################################

#define linear terms
params <- c("npp_annual",
            "popdens_hii_annual" ,
            "dens_all_roads_annual",
            "perc_nonveg_annual",
            "perc_tree_cov_annual",
            "elevation",
            "sl_",
            "log_sl_",
            "cos_ta_"
            )


#define quadratic terms
quad_params <- c("perc_tree_cov_annual")
                 


#prior params from https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.html?sequence=40#inla-1
mean.beta <- 0
prec.beta <- 1e-4 

#create quadratic terms if any
quad_terms <- vector()

for(i in 1:length(quad_params)){
  quad_terms[i] <- paste0("I(",quad_params[i],"^2)")
}

#quad_terms <- NULL

#create random effects terms for the INLA formula for the quadratic effects
# rand_terms_quad <- vector()
# 
# for(i in 1:length(quad_terms)){
#   rand_terms_quad[i] <- paste0("(0 + ", quad_terms[i], " | animal_id)")
# }

########################################################################
##
## 1. Import and format step data
##
########################################################################

#no imputation
steps <- read_csv("/Users/tb201494/Desktop/Panthera/pumas/projects/olympic_cougar_connectivity/data/2h_steps_unscaled_no_imp_annual_cov_02-13-2024.csv")

#round and define factors
steps <- steps %>% 
  mutate(across(c(perc_tree_cov_annual, perc_nonveg_annual, perc_nontree_veg_annual), round)) %>% 
  mutate(across(c(popdens_hii_annual, 
                  infra_hii_annual, 
                  landuse_hii_annual,
                  roads_hii_annual, 
                  hii_annual), function(x){return(round(x, digits=2))})) %>% 
  mutate(across(c(#animal_id, 
    sex, 
    dispersal_status, 
    year, 
    #case_, 
    tod_end_, 
    intersects_water, 
    intersects_study_area, 
    land_cover_usfs_annual, 
    land_use_usfs_annual, 
    land_use_change_usfs_annual, 
    dispersing, 
    disp_qual), as.factor)) 

#remove steps when dispersers were not dispersing and drop NA factor values
steps_filt <- steps %>% 
  filter(dispersal_status=="resident" | dispersing ==TRUE) %>% 
  mutate(land_cover_usfs_annual = case_when(land_cover_usfs_annual == 15 ~ NA,
                                            .default = land_cover_usfs_annual),
         land_use_usfs_annual = case_when(land_use_usfs_annual == 7 ~ NA,
                                          .default = land_use_usfs_annual),
         land_use_change_usfs_annual = case_when(land_use_change_usfs_annual == 5 ~ NA,
                                                 .default = land_use_change_usfs_annual)) %>% 
  select(-c(disp_qual:dispersing))



steps_scaled <- steps_filt %>% 
  mutate(across(c(elevation:precip_annual, popdens_hii_annual:dens_all_roads_annual, sl_, log_sl_, cos_ta_), scale)) 

scaling_attr <- steps_scaled %>% 
  select(matches(params)) %>%
  names() %>% 
  map_dfr(., ~ tibble(
    covariate = .x,
    center = attr(steps_scaled[[.x]], "scaled:center"),
    scale  = attr(steps_scaled[[.x]], "scaled:scale"))
    )
  
#write_csv(scaling_attr, "/Users/tb201494/Desktop/Panthera/pumas/projects/olympic_cougar_connectivity/data/fitted_models/scaling_parameters_top_ssf_mod_2-16-24.csv")


steps_scaled <- steps_scaled %>% 
  mutate(across(c(elevation:precip_annual, popdens_hii_annual:dens_all_roads_annual, sl_, log_sl_, cos_ta_), as.numeric))

#########################################################################
##
## 2. Fit global model in INLA
##
##########################################################################


n_indiv = steps_scaled %>% distinct(animal_id) %>% count() %>% pull()

all_terms <- paste(c(params, quad_terms))

rand_terms_inla <- vector()
for(i in 1:length(all_terms)){
  rand_terms_inla[i] <- paste0("f(", paste0("id", i), ", ", all_terms[i], ", ", "values = ", paste0("1:",n_indiv), ", ", 'model = "iid", ', "hyper = list(theta = list(initial = log(1), fixed = FALSE, ", 'prior = "pc.prec", ', "param = c(3, .05))))")
}


#convert response from categorical to numeric
dat <- steps_scaled %>% mutate(case_ = as.numeric(case_))

#create separate animal_id columns for each random effect
for(i in 1:length(rand_terms_inla)){
  name <- as.name(paste0("id", i))
  dat[[name]] <- as.numeric(factor(dat$animal_id))
}

#write_csv(dat, "top_ssf_dat_2-16-24.rds")

#construct the model formula with quadratic terms
form <- as.formula(paste("case_ ~ -1 + ",  
                         #fixed effects
                         paste(c(params, quad_terms), collapse = " + "), "+",
                         #random intercept (strata)
                         "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T))) + ",
                         #random slopes
                         paste(rand_terms_inla, collapse = " + ")))

  
#fit the model    
fit <-  inla(form, family ="Poisson", data=dat,
             control.fixed = list(
               mean = mean.beta,
               prec = list(default = prec.beta)),
             control.compute = list(waic=TRUE, dic = TRUE, cpo = TRUE, config = TRUE)) #
  
#saveRDS(fit, "top_ssf_2-16-24.rds")

#########################################################################
##
## 3. K-fold cross-validation
##
##########################################################################
cov_stack <- rast("/Users/tb201494/Desktop/prediction_rasters_2023/ssf_cov_stack_pred_scaled_hard_capped_01-24-2024.tif")

#only select cov stack layers relevant for prediction (terms in model)
#cov_stack_sel <- cov_stack[[params]]

#28 min for 4 folds and 8 covariates
#41 min with 8 linear and 1 quad, 5 folds
system.time(cv <- k_fold_inla(dat, form, cov_stack, n_folds = 4, plots = TRUE, seed = 777)); beep("fanfare")

#global 4-fold cv =  1,  0.9636364,  0.9878788,  0.9151515
#mean rho = 0.9666667
#rho = 0.9787879 with movement parameters 
cv

mean(unlist(cv[1:4]))

#########################################################################
##
## 4. External validation with new locations
##
##########################################################################

#extract coefficient names and values from the fitted model
mod_sum <- summary(fit6)
coeffs <- mod_sum$fixed %>% as.data.frame() #%>% rownames_to_column("param")
names <- rownames(coeffs)
names <- names[str_detect(names, coll("sl_"))==FALSE]
names <- names[str_detect(names, coll("ta_"))==FALSE]

names_linear <- names[str_detect(names, coll("2"))==FALSE]


cov_stack_sel <- cov_stack[[names_linear]]

#Generate predictions by multiplying each covariate layer by its coefficient estimate from the model
pred <- list()
for(j in 1:length(names)){
  
  cf <- coeffs[names[j], "mean"]
  
  if( cf == 0){
    cf <- cf + (.Machine$double.eps * 1e11) #add random jitter if mean coeff is = 0
  }
  
  if(str_detect(names[j], coll("2")) == T){
    tmp_str1 <- str_remove(names[j], coll("I("))
    tmp_str2 <- str_remove(tmp_str1, coll("^2)"))
    pred[[j]] <- cov_stack_sel[[tmp_str2]]^2 * cf
  } else{
    pred[[j]] <- cov_stack_sel[[names[j]]] * cf #breaks when mean coeff is 0
  }
  print(paste0(j, "/", length(names)))
}
predictions <- Reduce("+", pred) #sum layers together
pred_vals <- terra::values(predictions) #convert raster to values for calculations
preds_exp <- plogis(pred_vals) #backtransform predictions from logit to prob scale
x.min <- quantile(preds_exp, 0.025, na.rm = T) # omit the extremes
x.max <- quantile(preds_exp, 0.975, na.rm = T)
preds_exp[preds_exp > x.max] <- x.max
preds_exp[preds_exp < x.min] <- x.min
x.range <- x.max - x.min
predictions. <- (preds_exp - x.min)/x.range
vals <- cov_stack_sel[[1]] %>% terra::setValues(predictions.)

breaks <- quantile(predictions., probs = 0:10/10, na.rm = T) #obtain 10 quantile values for preds
breaks_j <- breaks + (seq_along(breaks) * .Machine$double.eps) #add jitter to breaks
binned <- terra::classify(vals, rcl=as.vector(breaks_j)) #classify the map based on breaks

#Subset used steps from the test data
test_steps <- read_csv("data/Location_Data/Source_Files/locations_master/test_locs_er_may23_feb24_dop_screened_02-08-24.csv")

#convert end points of used steps to sf points and extract values from the binned raster predictions
test_sf <- st_as_sf(test_steps, coords = c("lon", "lat"), crs=4326) %>% 
  st_transform(crs = 5070)

test_vals <- as.numeric(terra::extract(binned, test_sf)[,2]) %>% na.omit()

#plot proportion of test points in each bin. 
ggplot() +
  geom_bar(aes(x=test_vals, y = after_stat(prop))) 

n_bins <- length(unique(test_vals))

#calculate pearson correlation between the bin number and the proportion of used locations in each bin
cor_cv_res <- cor.test(x = 1:n_bins, y = as.numeric(table(test_vals))/(length(test_vals)/n_bins), method = "spearman", exact = FALSE)$estimate

cor_cv_res 

