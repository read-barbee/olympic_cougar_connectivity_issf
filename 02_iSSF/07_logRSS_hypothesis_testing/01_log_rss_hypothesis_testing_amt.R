#### iSSF Hypothesis Testing ####

# Author: Read Barbee

# Date:2023-09-15 

# Purpose: Examine log-RSS responses of each individual to each covariate and look for differences between sexes and dispersal groups


################################ Libraries #################################
library(tidyverse)
library(terra)
library(amt)
library(janitor)
library(GGally)
library(flextable)


#########################################################################
##
##  Specify Model Parameters
##
##########################################################################

params <- c("npp_annual",
            "popdens_hii_annual" ,
            "dens_all_roads_annual",
            "perc_nonveg_annual",
            "perc_tree_cov_annual",
            "elevation")

quad_terms <- paste0("I(", params, "^2)")

#########################################################################
##
##  Functions
##
##########################################################################

#run the global model
run_global <- function (dat)
  {
  #library(survival)  survival::clogit
  mod <- amt::fit_issf(formula = as.formula(paste("case_ ~", 
                                                  paste(c(params, quad_terms), collapse = " + "),
                                                  " + sl_ + log_sl_ + cos_ta_ + strata(step_id_)")),
                       data = dat,
                       na.action = "na.omit",
                       model = TRUE)
  return(mod)
}

#generate dataframe 1 for log_RSS calculation
set_df_s1 <- function ()
{
  #data frame varying elevation from min value to max value encountered by individual, holding all other covariates at the mean
  s1 <- data.frame(matrix(0, nrow = 200, ncol = length(params) + 3))
  colnames(s1) <- c(params, "sl_", "log_sl_", "cos_ta_")
  
  s1$sl_ <- 100
  
  s1$log_sl_ <- log(100)
  
  s1$cos_ta_ <- 1
  
  return(s1)
}

#generage dataframe 2 for log_RSS calculation
set_df_s2 <- function ()
{
  #data frame varying elevation from min value to max value encountered by individual, holding all other covariates at the mean
  s2 <- data.frame(matrix(0, nrow = 1, ncol = length(params) + 3))
  colnames(s2) <- c(params, "sl_", "log_sl_", "cos_ta_")
  
  s2$sl_ <- 100
  
  s2$log_sl_ <- log(100)
  
  s2$cos_ta_ <- 1
  
  return(s2)
}

#classify results as mostly positive, mostly negative, or unclear
classify_results <- function (curr_param, steps)
{
  classifications <- data.frame(animal_id = NA, uncertain = NA,
                                positive = NA)
  classifications_ <- classifications
  for (i in 1:nrow(steps)) {
    mod <- steps$global_fit[[i]]
    data <- steps$steps[[i]]
    
    rss_name <- as.name(paste0(curr_param, "_rss"))
    
    l_rss <- steps[i, rss_name][[1]] %>% as.data.frame()
    classifications_$animal_id <- steps$animal_id[[i]]
    
    count_pos <- l_rss %>% summarize(count_pos = sum(log_rss >0))
    
    classifications_$positive <- ifelse(count_pos > 100,T, F)
    
    count_overlap <- l_rss %>% 
      rowwise %>% 
      dplyr::summarize(count_overlap = sum(lwr <=0 && upr >= 0)) %>% 
      dplyr::filter(count_overlap == 1)
    
    count_overlap <- nrow(count_overlap)
    
    classifications_$uncertain <- ifelse(count_overlap > 100, T, F)
    
    classifications <- rbind(classifications, classifications_)
  }
  classifications <- classifications[-1, ]
  classifications <- classifications %>% dplyr::group_by(animal_id) %>%
    dplyr::arrange(animal_id)
  classifications$plot.order <- 1:nrow(classifications)
  return(classifications)
}

#calculate log_rss object for elevation for each individual incorporating quadratic terms, sl and ta
l_rss <- function(indiv, dat, curr_param)
  {
  indiv_dat <- dat %>% 
    filter(animal_id == indiv) %>% 
    unnest(cols=c(steps)) #%>% 
    #na.omit()
  
  s1 <- set_df_s1()
  s2 <- set_df_s2()
  
  cp_name <- as.name(curr_param)
  
  s1[,curr_param] <- seq(from = min(indiv_dat[,cp_name], na.rm=T), to = max(indiv_dat[,cp_name], na.rm=T), length.out = 200)
  
  
  indiv_dat_nested <- dat %>%
    filter(animal_id == indiv)
  
  model <- indiv_dat_nested$global_fit[[1]]
  
  ### Working. variable names have to be the same across all data frames and model
  l_rss_indiv <- amt::log_rss(model, s1, s2, ci = "se", ci_level = 0.95)
  
  return(l_rss_indiv$df)
}

#########################################################################
##
## 1. Import and format step data
##
##########################################################################

#steps
steps <- read_csv("data/gps_data/05_steps/2h_steps_unscaled_no_imp_annual_cov_02-13-2024.csv") %>% 
  mutate(step_id_ = paste0(animal_id, "_", step_id_))

#########################################################################
##
## 2. Reduce land cover and land use categories
##
##########################################################################
steps <- steps %>% 
  nest_by(animal_id, sex, dispersal_status) %>% 
  rename(steps=data)

#plot_bar(steps %>% select(gpp:calving_season))

#########################################################################
##
## 3. Fit global iSSF model to each individual
##
##########################################################################

#50 convergence issues with imputed set for clogit. 19 warnings for fit_issf
global_fits <- steps %>%  
  pull(steps) %>% 
  map(run_global)

steps$global_fit <- global_fits


#########################################################################
##
## 4. Calculate log-RSS and classify individual responses to each covariate
##
##########################################################################

indivs <- steps %>% pull(animal_id)


for(i in 1:length(params)){
  name <- as.name(paste0(params[i], "_rss"))
  steps[[name]] <- map(indivs, l_rss, dat=steps, curr_param=params[i])
  print(paste0(i,"/",length(params)))
}


#convert steps to datatable for ease of processing
steps2 <- data.table::as.data.table(steps)


#########################################################################
##
## 5. Classify RSS Values
##
##########################################################################

#classify individual log_rss trends as per Sells et al 2022
classified <- map(params, classify_results, steps=steps)
names(classified) <- params

#extract sex and dispersal infromation from steps
metadat <- steps %>% select(animal_id, sex, dispersal_status)

#join sex and dispersal status to log_rss classification tables
for(i in 1:length(classified)){
  classified[[i]] <- classified[[i]] %>% left_join(metadat, by=join_by(animal_id)) %>% 
    relocate(c(sex, dispersal_status), .after = animal_id)
}

classified

#########################################################################
##
## 6. Summarize RSS trends
##
##########################################################################

rows <- list()
names <- names(classified)

for (i in 1:length(classified)){
  classes <- classified[[i]] %>% mutate(class =case_when(uncertain==TRUE ~ "Uncertain",
                                                         uncertain==FALSE & positive==TRUE ~ "Positive",
                                                         uncertain==FALSE & positive==FALSE ~ "Negative")) 
  
  class_perc <- classes %>% group_by(sex, dispersal_status, class) %>% summarize(count =n()) %>% mutate(freq = (count / sum(count))) %>% 
    mutate(perc = round((freq * 100), digits=2))
  
  rows[[i]] <- class_perc %>% 
    mutate(perc_count = paste0(perc, " (",count,")")) %>% 
    unite("class", sex:class, sep="_") %>% 
    select(-c(count:perc)) %>% 
    pivot_wider(names_from = class, values_from = perc_count) %>% 
    mutate(name = names[i]) %>% 
    select(name, everything())
}

summary_table <- bind_rows(rows) %>% 
  mutate(across(everything(), ~replace_na(., "0.0"))) %>% 
  select(name,
         Female_resident_Positive, Female_resident_Negative, Female_resident_Uncertain,
         Male_resident_Positive, Male_resident_Negative, Male_resident_Uncertain,
         Female_disperser_Positive, Female_disperser_Negative, Female_disperser_Uncertain,
         Male_disperser_Positive, Male_disperser_Negative, Male_disperser_Uncertain)
  
#write_csv(summary_table, "/Users/tb201494/Desktop/Panthera/pumas/projects/olympic_cougar_connectivity/rss_summary_no_imp_9-17-25.csv")


tab <- flextable(summary_table)

#save_as_docx(tab, path = "/Users/tb201494/Desktop/Panthera/pumas/projects/olympic_cougar_connectivity/rss_summary_no_imp_9-17-25.docx")

#########################################################################
##
##7. Make RSS plots
##
##########################################################################


#pivot data for faceting
rss_plot_dat <- steps2 %>% 
  pivot_longer(npp_annual_rss:elevation_rss, names_to= "cov", values_to = "rss_val")


#initialize blank lists for loop
ls = list()
ls2=list()


#extract covariate ranges and log_rss values for each covariate and individual

for(i in 1:nrow(rss_plot_dat)){
  cov <- str_remove(rss_plot_dat$cov[i], "_rss")
  name <- as.name(paste0(cov, "_x1"))
  ls[[i]] <- rss_plot_dat$rss_val[[i]][[name]]
  ls2[[i]] = rss_plot_dat$rss_val[[i]]$log_rss
  
  print(paste0(i,"/",nrow(rss_plot_dat)))
}


#add lists from for loop to data frame and facet plot by sex (Melodie's RSS values don't make any sense)

rss_plot_dat$cov_vals <- ls
rss_plot_dat$rss_vals <- ls2

rss_plot_dat <- rss_plot_dat %>% 
  mutate(dem = paste0(sex, "_", dispersal_status)) %>% 
  mutate(cov = fct_recode(cov,
                           "Net Primary Productivity (NPP)" = "npp_annual_rss",
                           "Human population density (people/km2)" = "popdens_hii_annual_rss",
                           "Road density (km/km2)" = "dens_all_roads_annual_rss",
                           "% Non-vegetated area" = "perc_nonveg_annual_rss",
                           "% Tree cover" = "perc_tree_cov_annual_rss",
                           "Elevation (m)" = "elevation_rss")) %>% 
  mutate(cov = factor(cov, levels = c("Elevation (m)",
                                            "Net Primary Productivity (NPP)",
                                            "% Tree cover",
                                            "(% Tree cover)^2",
                                            "% Non-vegetated area",
                                            "Road density (km/km2)",
                                            "Human population density (people/km2)"
  ))) %>% 
  mutate(dem = fct_recode(dem,
                          "Male resident" = "Male_resident",
                          "Male disperser" = "Male_disperser",
                          "Female resident" = "Female_resident",
                          "Female disperser" = "Female_disperser")) %>% 
  mutate(dem = factor(dem, levels = rev(c("Male resident",
                                          "Male disperser",
                                          "Female resident",
                                          "Female disperser"))))
  
rm(global_fits)
rm(ls)
rm(ls2)
rm(steps)

#plot of individual log-rss curves for each univariate model by sex
rss_uni_sex <- rss_plot_dat %>% 
  # filter(animal_id!="Melodie") %>% #remove individuals with outlying values skewing plots
  # filter(animal_id!="Sampson") %>%
  # filter(animal_id!="Kingsley") %>%
  # filter(animal_id!="Bunny") %>%
  select(-rss_val) %>% 
  unnest(cols=c(cov_vals, rss_vals)) %>% 
  ggplot(., aes(x = cov_vals, y = rss_vals)) +
  geom_smooth(aes(col=sex, pch=animal_id), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Covariate Values") +
  ylab("log-RSS vs Mean Covariate Value") +
  theme_bw() +
  facet_wrap(~cov, scales = "free")

plotly::ggplotly(rss_uni_sex)

#ggsave(filename= "mf_rss_uni__5-06-2023.png", plot= rss_uni_sex)

#plot of individual log-rss curves for each univariate model by dispersal status
rss_uni_disp <- rss_plot_dat %>% 
  # filter(animal_id!="Melodie") %>% #remove individuals with outlying values skewing plots
  # filter(animal_id!="Sampson") %>%
  # filter(animal_id!="Kingsley") %>%
  # filter(animal_id!="Bunny") %>%
  select(-rss_val) %>% 
  unnest(cols=c(cov_vals, rss_vals)) %>% 
  ggplot(., aes(x = cov_vals, y = rss_vals)) +
  geom_smooth(aes(col=dispersal_status, pch=animal_id), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Covariate Values") +
  ylab("log-RSS vs Mean Covariate Value") +
  theme_bw() +
  facet_wrap(~cov, scales = "free")

plotly::ggplotly(rss_uni_disp)

#ggsave(filename= "disp_rss_uni_5-06-2023.png", plot= rss_uni_disp)


#faceted by dispersal status
rss_disp_facet_uni <- rss_plot_dat %>% 
  # filter(animal_id!="Melodie") %>% #remove individuals with outlying values skewing plots
  # filter(animal_id!="Sampson") %>%
  # filter(animal_id!="Kingsley") %>%
  # filter(animal_id!="Bunny") %>%
  select(-rss_val) %>% 
  unnest(cols=c(cov_vals, rss_vals)) %>% 
  filter(cov!="landuse_rss") %>% 
  ggplot(., aes(x = cov_vals, y = rss_vals)) +
  geom_smooth(aes(pch=animal_id, color=dispersal_status),linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Covariate Values") +
  ylab("log-RSS vs Mean Covariate Value") +
  theme_bw() +
  facet_wrap(vars(cov, dispersal_status), scales = "free")

plotly::ggplotly(rss_disp_facet_uni)

#ggsave(filename= "disp_rss_facet_imp_7-14-2023.png", plot= rss_disp_facet_uni)

#faceted by sex
rss_sex_facet_uni <- rss_plot_dat %>% 
  #filter(animal_id!="Melodie") %>% #remove individuals with outlying values skewing plots
  #filter(animal_id!="Sampson") %>%
  #filter(animal_id!="Kingsley") %>%
  #filter(animal_id!="Bunny") %>%
  select(-rss_val) %>% 
  unnest(cols=c(cov_vals, rss_vals)) %>% 
  filter(cov!="landuse_rss") %>% 
  ggplot(., aes(x = cov_vals, y = rss_vals)) +
  geom_smooth(aes(pch=animal_id, color=sex),linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Covariate Values") +
  ylab("log-RSS vs Mean Covariate Value") +
  theme_bw() +
  facet_wrap(vars(cov, sex), scales = "free")

plotly::ggplotly(rss_sex_facet_uni)

#ggsave(filename= "sex_rss_facet_imp_7-14-2023.png", plot= rss_sex_facet_uni)


#plot of individual log-rss curves for each univariate model by sex and disp status
rss_uni_sex_disp <- rss_plot_dat %>% 
  filter(animal_id!="Hana") %>% #remove individuals with outlying values skewing plots
  filter(animal_id!="Medusa") %>% 
  filter(animal_id!="TaTaw") %>% 
  filter(animal_id!="Jamestown") %>% #strangely positive reaction to road density
  select(-rss_val) %>% 
  unnest(cols=c(cov_vals, rss_vals)) %>% 
  ggplot(., aes(x = cov_vals, y = rss_vals)) +
  geom_smooth(aes(col=dem, pch=animal_id), linewidth = 1) + 
  scale_color_manual(values = rev(c("#0072B2", "#56B4E9", "#D55E00", "#E69F00"))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Covariate Values") +
  ylab("log-RSS vs Mean Covariate Value") +
  theme_minimal() +
  guides(colour = guide_legend(reverse = TRUE)) +
  facet_wrap(~cov, scales = "free")

#plotly::ggplotly(rss_uni_sex_disp)


ggsave(plot = rss_uni_sex_disp, "/Users/tb201494/Desktop/Panthera/pumas/projects/olympic_cougar_connectivity/figures/log_rss_spaghetti_by_demography_v3_9-17-25.png", width = 10, height = 5, dpi = 300)

