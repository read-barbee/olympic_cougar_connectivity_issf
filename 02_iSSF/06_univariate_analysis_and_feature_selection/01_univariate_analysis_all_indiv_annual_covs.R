#### OCP iSSF Module_04: Univariate Analysis ####

# Author: Read Barbee

# Date:2023-09-27 

#updated: 2025-10-27

# Purpose: Univariate feature selection for step selection models


################################ Libraries #################################
library(tidyverse)
library(amt)
library(janitor)
library(DataExplorer)
library(GGally)
library(sjPlot)
library(doParallel)
library(beepr)
library(INLA)
library(terra)


#register doParallel backend. Note:the cores argument implies forking with doMC backend, but you can check with the getDoParName() function
#doParallel::registerDoParallel(cores = 9) 

source("r_scripts_clean/02_iSSF/09_cross_validation/k_fold_cv_function_INLA_v4_quad_02-06-24.R")

set.seed(777)

#########################################################################
##
## 1. Import  data
##
##########################################################################

#steps
steps <- read_csv("data/gps_data/05_steps/2h_steps_unscaled_no_imp_annual_cov_02-13-2024.csv")

#cov stack for CV predictions
pred_stack <- rast("/Volumes/Reads_Seagate/ms_thesis/covariates/prediction_rasters/uncapped/pred_stack_2023_raw_uncapped.tif")

###############################################################################

#2. Format steps

###############################################################################
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


#########################################################################
##
## 3. Check covariate distributions
##
##########################################################################

#data overview
introduce(steps_filt)
plot_intro(steps_filt)

#check for missing values
plot_missing(steps_filt)


#continuous histograms all
plot_histogram(steps_filt %>% select(sl_, ta_, elevation:dens_all_roads_annual))

#continuous histograms used
plot_histogram(steps_filt %>% select(case_, sl_, ta_, elevation:dens_all_roads_annual) %>% filter(case_==TRUE))

#continuous histograms unused
plot_histogram(steps_filt %>% select(case_, sl_, ta_, elevation:dens_all_roads_annual) %>% filter(case_==FALSE))

#categorical bar plots
plot_bar(steps_filt %>% select(land_cover_usfs_annual, land_use_usfs_annual, land_use_change_usfs_annual))

plot_bar(steps_filt %>% select(land_cover_usfs_annual, land_use_usfs_annual, land_use_change_usfs_annual, dispersal_status), by="dispersal_status")

#full report with all DataExplorer metrics
#create_report(steps, y="case_")

#used and available step distributions overlaid on one another
steps_pivot <- steps_filt %>% 
  select(c(animal_id:dispersal_status, case_, elevation:precip_annual, popdens_hii_annual:dens_all_roads_annual, -aspect)) %>% 
  pivot_longer(elevation:dens_all_roads_annual, names_to = "cov", values_to = "value")

ggplot(steps_pivot, aes(x = value, y = ..density.., fill = as.factor(case_)))+ geom_histogram(position="identity", alpha=0.7) + 
  facet_wrap(vars(cov), scales = "free") +
  xlab("Used vs Available Distributions") + theme(axis.title.x=element_text(size=16))


#########################################################################
##
## 4. Scale continuous covariates for model comparison and correlation analysis
##
##########################################################################

steps_scaled <- steps_filt %>% 
  mutate(across(c(elevation:precip_annual, popdens_hii_annual:dens_all_roads_annual), scale)) %>% 
  mutate(across(c(elevation:precip_annual, popdens_hii_annual:dens_all_roads_annual), as.numeric))

#plot_histogram(steps_scaled %>% select(sl_, ta_, gpp:power_hii))

#plot_boxplot(steps_scaled %>% select(case_, sl_, ta_, gpp:power_hii), by = "case_")

#########################################################################
##
## 5. Identify highly correlated/redundant covariates
##
##########################################################################

#Correlation plot
plot_correlation(steps_scaled %>% 
                   select(sl_, ta_, elevation:dens_all_roads_annual) %>% 
                   select_if(is.numeric) %>% 
                   na.omit())

#select covariates of interest and dummy code factors
cont <- steps_scaled %>% 
  select(elevation:dens_all_roads_annual, -aspect) %>% 
  dummify(select = c("land_cover_usfs_annual", "land_use_usfs_annual", "land_use_change_usfs_annual")) %>%
  mutate(across(elevation:dens_all_roads_annual, as.numeric)) %>% 
  select(!contains("NA"))
  

#create a correlation matrix
cor_matrix <- cor(cont, use="pairwise.complete.obs", method = "pearson")

# Set correlation threshold
threshold <- 0.0

high_cor_pairs <- which(abs(cor_matrix) >= threshold, arr.ind = TRUE)
high_cor_pairs <- high_cor_pairs[high_cor_pairs[, 1] != high_cor_pairs[, 2], ]

var1 <- vector()
var2 <- vector()
correlation <- vector()

for (i in 1:nrow(high_cor_pairs)) {
  var1[i] <- rownames(cor_matrix)[high_cor_pairs[i, 1]]
  var2[i] <- colnames(cor_matrix)[high_cor_pairs[i, 2]]
  correlation[i] <- cor_matrix[high_cor_pairs[i, 1], high_cor_pairs[i, 2]]
  #print(paste("Variables:", var1, "and", var2, "- Correlation:", correlation))
}

#ordered list of correlated covariates
cor_dat <- tibble(variables = paste0(var1, "_", var2), corr = correlation) %>% 
  #filter(corr < 0.977) %>% 
  arrange(desc(corr))

rows_to_remove <- seq(from = 2, to = nrow(cor_dat), by = 2)
cor_dat_no_dupes <- cor_dat[-rows_to_remove, ]

#write_csv(cor_dat_no_dupes, "feature_selection/ssf_pairwise_cov_corr_no_imp_all_indiv_annual_cov_02-02-24.csv")

#the GGAlly package corrplot is more informative for numeric variables
ggcorr(steps_scaled %>% 
         select(elevation:dist_all_roads_annual), label = TRUE)


#########################################################################
##
## 6. Prep data for Univariate Muff models (INLA)
##
##########################################################################
#subset dataframe to only include animal_id, case_, and all predictors
covs <- steps_scaled %>% 
  select(case_, animal_id, step_id_, x2_, y2_, sl_, log_sl_, cos_ta_, elevation:dens_all_roads_annual, -aspect) %>%
  mutate(across(c(case_, land_cover_usfs_annual:land_use_change_usfs_annual), as.numeric)) %>% 
  mutate(animal_id = as.numeric(as.factor(animal_id))) %>% 
  mutate(animal_id2 = as.numeric(as.factor(animal_id))) %>% 
  mutate(animal_id3 = as.numeric(as.factor(animal_id))) %>% 
  mutate(animal_id4 = as.numeric(as.factor(animal_id))) %>% 
  mutate(animal_id5 = as.numeric(as.factor(animal_id)))


#########################################################################
##
## 7. Fit null model
##
##########################################################################

null_form <- case_ ~ -1 + sl_ + log_sl_ + cos_ta_ + 
  f(step_id_, 
    model = "iid",
    hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
  f(animal_id, sl_, values = 1:84, model = "iid", 
  hyper = list(theta = list(initial = log(1), fixed = FALSE, 
                            prior = "pc.prec", param = c(3, 0.05)))) +
  f(animal_id2, log_sl_, values = 1:84, model = "iid", 
    hyper = list(theta = list(initial = log(1), fixed = FALSE, 
                              prior = "pc.prec", param = c(3, 0.05)))) +
  f(animal_id3, cos_ta_, values = 1:84, model = "iid", 
    hyper = list(theta = list(initial = log(1), fixed = FALSE, 
                              prior = "pc.prec", param = c(3, 0.05))))


run_null <- function(dat){
  #convert response from categorical to numeric
  dat <- dat %>% mutate(case_ = as.numeric(case_))
  
  form <- null_form
  
  fit <- inla(form, family = "poisson", data = dat, control.compute=list(dic=TRUE, 
                                                                         cpo=TRUE, 
                                                                         waic=TRUE, 
                                                                         return.marginals.predictor=TRUE))
  
  return(fit)
}

null_fit <- run_null(covs)

# null_cv <- k_fold_inla(dat = covs, form = null_form, cov_stack = pred_stack, n_folds = 4, seed=777)

null_fit <- list(null_fit)

names(null_fit) <- "null"


#########################################################################
##
## 8. Fit univaraite models for continuous covariates
##
##########################################################################

cov_names <- covs %>% 
  select(-c(case_:y2_, sl_, log_sl_, cos_ta_, contains("animal_id"))) %>% 
  select(!contains('usfs')) %>% #remove categorical variables
  names()

#get number of individuals from data
n_indiv = covs %>% distinct(animal_id) %>% count() %>% pull()

#initialize lists to store model formulas and waic scores
mod_names <- list()
waic_list <- list()
dic_list <- list()
cpo_list <- list()
cv_list <- list()
beta_est <- list()
low_95 <- list()
high_95 <- list()
low_90 <- list()
high_90 <- list()

#fit univariate model for each covariate
system.time(for (i  in 1:length(cov_names)) { 
  
  cov <- cov_names[i]
  
  #construct the model formula with linear terms only
  form <- as.formula(paste("case_ ~ -1 + sl_ + log_sl_ + cos_ta_ + ",  
                           #fixed effects
                           cov, "+",
                           #random intercept (strata)
                           "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6), fixed=T))) + ",
                           #random slopes
                           paste0("f(animal_id, sl_,", 
                                  " values = ", 
                                  paste0("1:",n_indiv), 
                                  ", model = \"iid\", hyper = list(theta = list(initial = log(1),",
                                  "fixed = FALSE, prior = \"pc.prec\", param = c(3, .05)))) + "),
                           paste0("f(animal_id2, log_sl_,", 
                                  " values = ", 
                                  paste0("1:",n_indiv), 
                                  ", model = \"iid\", hyper = list(theta = list(initial = log(1),",
                                  "fixed = FALSE, prior = \"pc.prec\", param = c(3, .05)))) + "),
                           paste0("f(animal_id3, cos_ta_,", 
                                  " values = ", 
                                  paste0("1:",n_indiv), 
                                  ", model = \"iid\", hyper = list(theta = list(initial = log(1),",
                                  "fixed = FALSE, prior = \"pc.prec\", param = c(3, .05)))) + "),
                           paste0("f(animal_id4, ", cov, 
                                  ", values = ", 
                                  paste0("1:",n_indiv), 
                                  ", model = \"iid\", hyper = list(theta = list(initial = log(1),",
                                  "fixed = FALSE, prior = \"pc.prec\", param = c(3, .05))))")))
  
  
mod <- inla(form, family = "poisson", data = covs, control.compute=list(dic=TRUE, 
                                                                         cpo=TRUE, 
                                                                         waic=TRUE, 
                                                                         return.marginals.predictor=FALSE,
                                                                         config = TRUE))
 
 mod_names[[i]] <- paste(mod$names.fixed, collapse="+")
 
 beta_est[[i]] <- mod$summary.fixed$mean
 
 low_95[[i]] <- mod$summary.fixed$`0.025quant`
 high_95[[i]] <- mod$summary.fixed$`0.975quant`
 
 posterior_samples <- inla.posterior.sample(mod, n = 1000)
 posterior_samples2 <- as.numeric(inla.posterior.sample.eval(paste0(cov), posterior_samples))
 
 low_90[[i]] <- quantile(posterior_samples2, c(0.05, 0.95))[1]
 high_90[[i]] <- quantile(posterior_samples2, c(0.05, 0.95))[2]
 
 waic_list[[i]] <- mod$waic$waic
 
 dic_list[[i]] <- mod$dic$dic
 
 cpo_list[[i]] <- mean(mod$cpo$cpo)
 
 k_fold <- k_fold_inla(dat = covs, form = form, cov_stack = pred_stack, n_folds = 4, plots = FALSE, seed = 777)
 
 cv_list[[i]] <- mean(k_fold)
  
  print(paste0(i, "/", length(cov_names)))
 
}) ; beep("fanfare")


mod_names <- unlist(mod_names)
beta_est <- unlist(beta_est)
waic_list <- unlist(waic_list)
dic_list<- unlist(dic_list)
cpo_list <- unlist(cpo_list)
#cv_list <- unlist(cv_list)

cv_list2 <- lapply(cv_list, function(x) if (is.null(x)) NA else x)

beta_est2 <- beta_est[seq(4, length(beta_est), 4)]


uni_tab <- tibble(name = mod_names,
                  estimate = beta_est2,
                  waic = waic_list,
                  dic = dic_list,
                  cpo = cpo_list,
                  mean_4_fold_rho = unlist(cv_list2))

write_csv(uni_tab, "ssf_uni_fits_all_indiv_muff_movement_2-10-24.csv")


names(uni_fits) <- cov_names

ggplot(mod_table %>% arrange(mean_rho)) +
  geom_point(aes(x = mean_rho, y = waic))

ggplot(mod_table %>% arrange(cpo)) +
  geom_point(aes(x = cpo, y = waic))

#########################################################################
##
## 9. Fit univaraite models for categorical covariates
##
##########################################################################

covs_cat <- covs %>% 
  select(case_, animal_id, step_id_, land_cover_usfs_annual:land_use_change_usfs_annual) %>%
  mutate(across(land_cover_usfs_annual:land_use_change_usfs_annual, as.factor)) %>% 
  dummify(select = c("land_cover_usfs_annual", "land_use_usfs_annual", "land_use_change_usfs_annual")) %>% 
  select(!contains("NA")) %>% 
  rename(land_cover_usfs_trees = land_cover_usfs_annual_1,
         land_cover_usfs_tall_tree_shrub = land_cover_usfs_annual_2,
         land_cover_usfs_tree_shrub = land_cover_usfs_annual_3,
         land_cover_usfs_gfh_tree = land_cover_usfs_annual_4,
         land_cover_usfs_barren_tree = land_cover_usfs_annual_5,
         land_cover_usfs_tall_shrub = land_cover_usfs_annual_6,
         land_cover_usfs_shrubs = land_cover_usfs_annual_7,
         land_cover_usfs_gfh_shrub = land_cover_usfs_annual_8,
         land_cover_usfs_barren_shrub = land_cover_usfs_annual_9,
         land_cover_usfs_gfh = land_cover_usfs_annual_10,
         land_use_usfs_agr = land_use_usfs_annual_1,
         land_use_usfs_dev = land_use_usfs_annual_2,
         land_use_usfs_forest = land_use_usfs_annual_3,
         land_use_usfs_non_forest_wet = land_use_usfs_annual_4,
         land_use_usfs_other = land_use_usfs_annual_5,
         land_use_usfs_rangeland = land_use_usfs_annual_6,
         land_use_change_stable = land_use_change_usfs_annual_1,
         land_use_change_slow_loss = land_use_change_usfs_annual_2,
         land_use_change_fast_loss = land_use_change_usfs_annual_3,
         land_use_change_gain = land_use_change_usfs_annual_4)
  

cov_names_cat <- covs_cat %>% select(-c(case_:step_id_)) %>% names()

# ~12 minutes
uni_fits_cat <- list()
system.time(for (i  in 1:length(cov_names_cat)) {
  
  cov <- cov_names_cat[i]
  
  #construct the model formula with linear terms only
  form <- as.formula(paste("case_ ~ -1 + ",  
                           #fixed effects
                           cov, "+",
                           #random intercept (strata)
                           "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6), fixed=T))) + ",
                           #random slopes
                           paste0("f(animal_id, ", cov, 
                                  ", values = ", 
                                  paste0("1:",n_indiv), 
                                  ", model = \"iid\", hyper = list(theta = list(initial = log(1),",
                                  "fixed = FALSE, prior = \"pc.prec\", param = c(1, .05))))")))
  
  
  uni_fits_cat[[i]] <- inla(form, family = "poisson", data = covs_cat, control.compute=list(dic=TRUE, 
                                                                                    cpo=TRUE, 
                                                                                    waic=TRUE, 
                                                                                    return.marginals.predictor=TRUE))
  
  print(paste0(i, "/", length(cov_names_cat)))
  
}) ; beep("fanfare")


names(uni_fits_cat) <- cov_names_cat


#########################################################################
##
## 10. Fit univariate models with quadratic effects
##
##########################################################################

#takes about 13 minutes to fit 25 models;
cov_names_quad <- covs %>% select(-c(case_:step_id_, x2_, y2_, sl_, log_sl_, cos_ta_, land_cover_usfs_annual:land_use_change_usfs_annual, contains("animal_id"))) %>% names()

mod_names_quad <- list()
waic_list_quad <- list()
dic_list_quad <- list()
cpo_list_quad <- list()
cv_list_quad <- list()
beta_est_quad <- list()
low_95_quad <- list()
high_95_quad <- list()
low_90_quad <- list()
high_90_quad <- list()

#uni_fits_quad <- list() RESUME EDITING HERE
system.time(for (i  in 1:length(cov_names_quad)){
  
  cov <- cov_names_quad[i]
  
  form <- as.formula(paste("case_ ~ -1 + sl_ + log_sl_ + cos_ta_ + ",  
                           #fixed effects
                           cov, "+", "I(",cov,"^2) + ",
                           #random intercept (strata)
                           "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6), fixed=T))) + ",
                           #random slopes
                           paste0("f(animal_id, sl_,", 
                                  " values = ", 
                                  paste0("1:",n_indiv), 
                                  ", model = \"iid\", hyper = list(theta = list(initial = log(1),",
                                  "fixed = FALSE, prior = \"pc.prec\", param = c(3, .05)))) + "),
                           paste0("f(animal_id2, log_sl_,", 
                                  " values = ", 
                                  paste0("1:",n_indiv), 
                                  ", model = \"iid\", hyper = list(theta = list(initial = log(1),",
                                  "fixed = FALSE, prior = \"pc.prec\", param = c(3, .05)))) + "),
                           paste0("f(animal_id3, cos_ta_,", 
                                  " values = ", 
                                  paste0("1:",n_indiv), 
                                  ", model = \"iid\", hyper = list(theta = list(initial = log(1),",
                                  "fixed = FALSE, prior = \"pc.prec\", param = c(3, .05)))) + "),
                           paste0("f(animal_id4, ", cov, 
                                  ", values = ", 
                                  paste0("1:",n_indiv), 
                                  ", model = \"iid\", hyper = list(theta = list(initial = log(1),",
                                  "fixed = FALSE, prior = \"pc.prec\", param = c(1, .05)))) + "),
                           paste0("f(animal_id5, ", "I(", cov,  "^2)", 
                                  ", values = ", 
                                  paste0("1:",n_indiv), 
                                  ", model = \"iid\", hyper = list(theta = list(initial = log(1),",
                                  "fixed = FALSE, prior = \"pc.prec\", param = c(1, .05))))")))
  
  mod <- inla(form, family = "poisson", data = covs, control.compute=list(dic=TRUE, 
                                                                          cpo=TRUE, 
                                                                          waic=TRUE, 
                                                                          return.marginals.predictor=TRUE,
                                                                          config = TRUE))
  
  mod_names_quad[[i]] <- paste0(mod$names.fixed[1], "2")
  
  beta_est_quad[[i]] <- mod$summary.fixed[2,]$mean
  
  low_95_quad[[i]] <- mod$summary.fixed[2,]$`0.025quant`
  high_95_quad[[i]] <- mod$summary.fixed[2,]$`0.975quant`
  
  posterior_samples <- inla.posterior.sample(mod, n = 1000)
  posterior_samples2 <- as.numeric(inla.posterior.sample.eval(paste0("I(", cov, "^2)"), posterior_samples))
  
  low_90_quad[[i]] <- quantile(posterior_samples2, c(0.05, 0.95))[1]
  high_90_quad[[i]] <- quantile(posterior_samples2, c(0.05, 0.95))[2]
  
  waic_list_quad[[i]] <- mod$waic$waic
  
  dic_list_quad[[i]] <- mod$dic$dic
  
  cpo_list_quad[[i]] <- mean(mod$cpo$cpo)
  
  k_fold <- k_fold_inla(dat = covs, form = form, cov_stack = pred_stack, n_folds = 4, plots = FALSE, seed = 777)
  
  cv_list_quad[[i]] <- mean(k_fold)
  
  print(paste0(i, "/", length(cov_names_quad)))
  
}); beep("fanfare")


#rename qudaratic covariates with "2" at the end
quad_names <- vector()
for(i in 1:length(cov_names_quad)){
  old_name <- cov_names_quad[i]
  quad_names[i] <- paste0(old_name, "2")
}

names(uni_fits_quad) <- quad_names

mod_table_quad <- tibble(mod = unlist(mod_names_quad),
                         est = unlist(beta_est_quad),
                         waic = unlist(waic_list_quad),
                         dic = unlist(dic_list_quad),
                         cpo = unlist(cpo_list_quad)
                         ) #%>% #cv = unlist(cv_list_quad)
  cbind(cv_list_flat)

write_csv(mod_table_quad, "mod_table_quad.csv")


temp <- as.data.frame(unlist(cv_list_quad))

write_csv(temp, "uni_quad_cv.csv")

#combine linear and quadratic fits
uni_fits_all <- c(uni_fits, uni_fits_cat, uni_fits_quad, null_fit)


#########################################################################
##
## 11. Create model summary table
##
##########################################################################

mod_names_all <- names(uni_fits_all)



mod_table <- list()
for (i in 1:length(uni_fits_all)){
  name <- mod_names_all[i]
  
  if(str_detect(name, "2")){
    tmp1 <-  uni_fits_all[[i]]$summary.fixed %>% 
      rownames_to_column() %>%  
      as_tibble() %>% filter(str_detect(rowname, "2")==T)
    tau_mean = uni_fits_all[[i]]$summary.hyperpar$mean[2]
    tau_sd = uni_fits_all[[i]]$summary.hyperpar$sd[2]
  } else if (name=="null"){
    tmp1 <-  uni_fits_all[[i]]$summary.fixed %>% 
      rownames_to_column() %>%  
      as_tibble()
    tau_mean = NA
    tau_sd = NA
  } else{
    tmp1 <-  uni_fits_all[[i]]$summary.fixed %>% 
      rownames_to_column() %>%  
      as_tibble() %>% filter(str_detect(rowname, coll("("))==F)
    
    tau_mean = uni_fits_all[[i]]$summary.hyperpar$mean
    tau_sd = uni_fits_all[[i]]$summary.hyperpar$sd
  }
  
  tmp2 <-  tmp1 %>% 
    dplyr::mutate(term = name, .before = mean) %>% 
    dplyr::mutate(tau_mean = tau_mean,
                  tau_sd = tau_sd,
                  waic = uni_fits_all[[i]]$waic$waic,
                  p.eff_waic = uni_fits_all[[i]]$waic$p.eff,
                  dic = uni_fits_all[[i]]$dic$dic,
                  p.eff_dic = uni_fits_all[[i]]$dic$p.eff,
                  mean_cpo = mean(uni_fits_all[[i]]$cpo$cpo),
                  mean_pit = mean(uni_fits_all[[i]]$cpo$pit),
                  mlik = uni_fits_all[[i]]$mlik[1])
  
  if(name=="null"){
    mod_table[[i]] <- tmp2 %>% dplyr::mutate(tau_mean = NA, tau_sd = NA) %>% 
      dplyr::select(term, mean, everything(), -rowname)
  } else{
    mod_table[[i]] <-tmp2 %>% 
      dplyr::select(term, mean, everything(), -rowname)
    }
    
  
  print(i)
}

mod_table <- bind_rows(mod_table)

#remove categorical variables
mod_table_filt <- mod_table %>% filter(!(term %in% c("land_cover_usfs_annual", "land_use_usfs_annual", "land_use_change_usfs_annual", "land_cover_usfs_annual2", "land_use_usfs_annual2", "land_use_change_usfs_annual2")))


#write_csv(mod_table_filt "feature_selection/uni_muff_fits_no_imp_all_indiv_annual_covs_inla_11-30-23.csv")

#remove redundant linear terms from quadratic models
#muff_mod_summ_filt <- muff_mod_summ %>% filter(str_detect(term, "I") ==TRUE | str_detect(name, "2")==FALSE)

# write_csv(muff_mod_summ_filt, "feature_selection/uni_muff_summary_quad_mods_no_imp_residents_9-18-23.csv")

#########################################################################
##
## 12. Examine individual random slope estimates
##
##########################################################################

summary(ndvi_quad)

hist(ndvi_quad$summary.random$animal_id$mean)
hist(ndvi_quad$summary.random$animal_id2$mean)


test_dat <- ndvi_quad$summary.random$animal_id %>% as.data.frame() %>% rename(lower = `0.025quant`, upper = `0.975quant`) %>% 
  mutate(cross_zero = case_when(lower < 0 & upper > 0 ~ TRUE,
                                .default = FALSE))

ggplot(test_dat, aes(x = ID, y = mean, color = cross_zero)) +
  geom_point(size = 1.5, colour = "#112446") +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_hline(aes(yintercept = 0, color = "red"))+
  theme_minimal()


summary(slope_quad)

hist(slope_quad$summary.random$animal_id$mean)
hist(slope_quad$summary.random$animal_id2$mean)


slope_dat <- slope_quad$summary.random$animal_id %>% as.data.frame() %>% rename(lower = `0.025quant`, upper = `0.975quant`) %>% 
  mutate(cross_zero = case_when(lower < 0 & upper > 0 ~ TRUE,
                                .default = FALSE))

ggplot(slope_dat, aes(x = ID, y = mean, color = cross_zero)) +
  geom_point(size = 1.5, colour = "#112446") +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_hline(aes(yintercept = 0), color = "red")+
  theme_minimal()

covs2 <- covs %>% mutate(fitted = plogis(ndvi_quad$summary.fitted.values$mean))

#rsf style fits
hist(covs2$fitted, scale="frequency", breaks="Sturges", col="darkgray")
ggplot(covs2, aes(x = fitted, fill = as.factor(case_)))+ geom_histogram(binwidth=0.05, position="identity", alpha=0.7) + xlab("Predicted Probability of Lion Use") + theme(axis.title.x=element_text(size=16))

