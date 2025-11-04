#### Validate step count raster ####

# Author: Read Barbee

# Date:2025-10-31 

# Purpose: Validate step count raster
#Inputs: 
#1. Quantile-binned step count raster

#Outputs: 
#1. Validation plots and spearman rank correlations

################################ Libraries #################################
library(tidyverse)
library(terra)
library(sf)

################################ Import data #################################

#binned step count raster
binned <- rast("results/comparison_layers/ssf/simulation/ssf_sim_paths_binned_counts_200000p_4380s_both_sides_I5_30m_uncapped_10-31-25.tif")

#probability of use raster
pred_rast <- rast("results/predictive surfaces/ssf/predictions/pop_level_map_standardized_train_scaled_9-18-25.tif")

#expert probability of use
binned <- rast("results/comparison_layers/expert/wsdot_habitat_binned_3-14-24.tif")

#extract raster values
pred_vals <- terra::values(pred_rast)

#identify quantile breaks
breaks <- quantile(pred_vals, probs = 0:10/10, na.rm = T) #obtain 10 quantile 

#classify raster by breaks
binned <- terra::classify(pred_rast, rcl=as.vector(breaks), include.lowest=TRUE)

#set factor levels
levels(binned) <- 1:10

#test gps locations
locs <- read_csv("data/gps_data/06_test_data/test_dat_full_er_and_amt_rejects_2-20-24.csv")

#########################################################################
##
##  1. External validation with GPS locations
##
##########################################################################

locs <- locs %>% filter(animal_id!="Cato")

eval_res <- locs %>% st_as_sf(coords = c("utm_e", "utm_n"), crs = 5070)

names(binned) <- "binned"

#extract binned predicted probability of use values for each used location
loc_preds <- terra::extract(binned, eval_res, bind = T) |> 
  as_tibble() |> 
  filter(!is.na(binned)) |> 
  mutate(binned = as.factor(binned))

loc_preds2 <- loc_preds %>% 
  complete(binned, fill = list(perc = 0)) %>%        
  group_by(binned) %>% 
  count() %>% 
  mutate(perc = (n/nrow(loc_preds))) %>% 
  ungroup()

#count points in each bin by demography
loc_preds_dem <- loc_preds %>%
  mutate(dem_class = paste0(sex, "_", dispersal_status), .before = sex) %>%
  count(binned, dem_class) %>%
  group_by(dem_class) %>%
  mutate(
    total = sum(n),
    perc = n / total
  ) %>%
  ungroup() %>%
  mutate(dem_class = fct_recode(dem_class,
                                "Female Disperser" = "Female_disperser",
                                "Female Resident" = "Female_resident",
                                "Male Disperser" = "Male_disperser",
                                "Male Resident" = "Male_resident"))

#calculatae spearman correlations by demography
corr_labels <- loc_preds_dem %>%
  group_by(dem_class) %>%
  summarise(corr = cor(as.numeric(binned), perc, method = "spearman"), .groups = "drop") %>%
  mutate(label = paste0("r = ", round(corr, 2)))


#plot proportion of test points in each bin. 
ggplot(data = loc_preds2) +
  geom_col(aes(x = binned, y = perc, fill = binned)) +
  geom_text(
    aes(x = binned, y = perc, label = scales::percent(perc, accuracy = 0.01)),
    position = position_dodge(width = 1),
    check_overlap = TRUE,
    vjust = -0.5
  ) +
  scale_fill_manual(values = viridis::magma(10),
                    na.value = NA,
                    guide = guide_legend(reverse = FALSE),
                    na.translate = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = "iSSF Class", y = "Percentage of Locations") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 16))

ggsave("expert_uncapped_30m_ext_val_bar_10-31-25.png", width = 7, height = 5, units = "in")


 n_bins <- length(unique(loc_preds))

#calculate pearson correlation between the bin number and the proportion of used locations in each bin
cor_cv_res <- cor.test(x = as.numeric(loc_preds2$binned), y = loc_preds2$perc, method = "spearman", exact = FALSE)$estimate

cor_cv_res #0.96



#########################################################################
##
##  3. Pie charts following Sells et al. 2022
##
##########################################################################
loc_preds_df <- tibble(index = 1:length(loc_preds), vals = loc_preds) %>% 
  group_by(vals) %>% 
  summarize(prop = n()/nrow(.)) %>% 
  mutate(ypos = cumsum(prop)- 0.5*prop) %>% 
  mutate(csum = rev(cumsum(rev(prop))), 
         pos = prop/2 + lead(csum, 1),
         pos = if_else(is.na(pos), prop/2, pos))

#remove unnecessary percentage labels
prop2 <- loc_preds_df$prop
prop2[1:3] <- NA

ggplot(loc_preds_df, aes(x="", y=prop, fill=vals)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = viridis::magma(10), na.value = NA, guide = guide_legend(reverse = FALSE), na.translate=FALSE)+
  geom_text(aes(x=1.2, label = scales::percent(prop2, accuracy = 0.01)), position = position_stack(vjust =0.5), check_overlap=T, vjust = 0) + 
  guides(fill = guide_legend(title = "Group")) +
  scale_y_continuous(breaks = loc_preds_df$pos, labels = as.character(loc_preds_df$vals)) +
  coord_polar("y", start=0, direction = -1) +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 15), 
        legend.position = "none", # Removes the legend
        panel.background = element_rect(fill = "white"))#+
# ggtitle("Population: % Locations per iSSF Class 1 - 10") 

#ggsave("ssf_sim_100m_ext_val_pie_3-12-24.png", width = 7, height = 5, units = "in")


#########################################################################
##
##  3. Ggplot of binned simulated surface
##
##########################################################################

ggplot()+
  tidyterra::geom_spatraster(data=binned, mapping=aes(), maxcell = 700000) +
  scale_fill_manual(values = viridis::magma(10), na.value = NA, guide = guide_legend(reverse = FALSE), na.translate=FALSE) +
  coord_sf(xlim = c(-124.8, -121),  ylim = c(45.5, 48.8), crs = 4326)


#########################################################################
##
##  3. External Individual Validataion
##
##########################################################################

ext_locs <- locs %>% st_as_sf(coords = c("utm_e", "utm_n"), crs = 5070)

indiv_names <- ext_locs %>% distinct(animal_id) %>% pull(animal_id)

#extract percentages of points in each bin for each individual
indiv_percs <- list()
for(i in 1:length(indiv_names)){
  
  indiv <- indiv_names[i]
  
  eval_indiv <- eval_res %>% 
    filter(animal_id == indiv)
  
  indiv_percs[[i]] <- terra::extract(binned, eval_indiv)[,2] %>% 
    as.numeric() %>% 
    na.omit() %>% 
    as.factor() %>% 
    enframe() %>% 
    mutate(value = factor(value, levels = as.character(1:10))) %>% 
    group_by(value) %>% 
    count() %>% 
    ungroup() %>% 
    complete(value, fill = list(n=0), explicit = FALSE) %>%
    mutate(perc = (n/ nrow(eval_indiv))*100) %>% 
    ungroup() %>% 
    select(value, perc) %>% 
    rename(class=value,
           !!indiv := perc)
  
  print(paste0(i, "/", length(indiv_names)))
}

#combine list into single dataframe
indiv_perc_tab <- reduce(indiv_percs, function(x, y) left_join(x, y, by = "class"))

#format for ggplot
indiv_perc_gg <- indiv_perc_tab %>% 
  pivot_longer(-class, names_to = "animal_id", values_to = "perc") %>% 
  left_join(ext_locs %>% 
              as.data.frame() %>% 
              distinct(animal_id, .keep_all=T) %>% 
              select(animal_id, sex, dispersal_status), by="animal_id") %>% 
  relocate(c(sex, dispersal_status), .after=animal_id) %>% 
  arrange(animal_id, class) %>% 
  mutate(dispersal_status = str_to_title(dispersal_status))

#calcualte pop levcel stats for points and error bars
group_stats <- indiv_perc_tab %>% 
  group_by(class) %>% 
  summarize(quant_025 = quantile(c_across(everything()), 0.025),
            mean = mean(c_across(everything())),
            quant_975 = quantile(c_across(everything()), 0.975))


## Plot all individuals ##
ggplot() +
  geom_line(data = indiv_perc_gg, aes(x = class, y=perc, group = animal_id, col = animal_id)) + 
  geom_point(data = group_stats, aes(x = class, y = mean)) +
  geom_errorbar(data = group_stats, aes(x =class, ymin = quant_025, ymax = quant_975, width = 0.2)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "iSSF Class", y = "Percentage of Locations") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 16))

#plotly::ggplotly(plot)

ggsave("expert_uncapped_30m_ext_val_indiv_95_ci_10-31-25.png", width = 7, height = 5, units = "in")

## Facet by sex and dispersal_status ##
ggplot(data = indiv_perc_gg) +
  geom_line(aes(x = class, y=perc, group = animal_id, col = sex)) + 
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "iSSF Class", y = "Percentage of Locations") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 16)) +
  facet_wrap(vars(sex, dispersal_status))

ggsave("expert_uncapped_30m_ext_val_indiv_sex_disp_facet_10-31-25.png", width = 7, height = 5, units = "in")

## Facet by sex and dispersal_status (Bar) ##
ggplot(data = loc_preds_dem) +
  geom_col(aes(x = binned, y = perc, fill = binned)) +
  geom_text(
    aes(x = binned, y = perc, label = scales::percent(perc, accuracy = 0.01)),
    position = position_dodge(width = 1),
    check_overlap = TRUE,
    vjust = -0.5
  ) +
  scale_fill_manual(values = viridis::magma(10),
                    na.value = NA,
                    guide = guide_legend(reverse = FALSE),
                    na.translate = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = "iSSF Class", y = "Percentage of Locations") +
  theme_minimal() +
  facet_wrap(vars(dem_class), scales = "free_y") +
  geom_text(
    data = corr_labels,
    aes(x = -Inf, y = Inf, label = label),   # use panel-relative corners
    hjust = -0.3, vjust = 1.7,               # nudge inward from the corner
    size = 5,
    fontface = "bold"
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 16),
    strip.background = element_rect(fill = "grey85", color = NA)
  )




ggsave("ssf_sim_uncapped_30m_ext_val_sex_disp_facet_bar_10-31-25.png", width = 12, height = 7, units = "in")

#########################################################################
##
##  3. Habitat covariates within quantile bins
##
##########################################################################

covs_raw_23 <- rast("/Users/tb201494/Desktop/annual_covs_01-07-24/cov_stack_2023.tif")
covs_static <- rast("/Users/tb201494/Desktop/1km_buffer/static_stack_1km_buffer_11-29-23.tif")

new_names <- vector()
for(i in 1:length(names(covs_raw_23))){
  new_names[i] <- substr(names(covs_raw_23)[i], 1, nchar(names(covs_raw_23)[i]) -5)
}

names(covs_raw_23) <- new_names

covs_all <- c(covs_static, covs_raw_23)

#writeRaster(covs_all, "covs_all_2023_raw_uncapped_2-22-24.tif")

cov_sel_23 <- covs_all[[names_linear_no_movement]]

#levels(binned) <- 1:10

cov_sel_200 <- terra::resample(cov_sel_23, binned)

binned_comb <- c(binned, cov_sel_200)

plot_df <- as_tibble(binned_comb) %>% 
  rename(class = sum) %>% 
  filter(!is.na(class)) %>% 
  pivot_longer(-class, names_to = "cov", values_to = "value") %>% 
  mutate(cov = as.factor(str_remove(cov, coll("_annual")))) %>% 
  group_by(cov)%>%
  mutate(median = median(value, na.rm=T)) %>%
  ungroup() %>% 
  mutate(cov = fct_recode(
    cov,
    "Elevation (m)" = "elevation",
    "NPP" = "npp",
    "% Non-vegetated Area" = "perc_nonveg",
    "% Tree Cover" = "perc_tree_cov",
    "Population Density (HII)" = "popdens_hii",
    "Road Density (km/km2)" = "dens_all_roads"
  )) %>% 
  mutate(cov=fct_relevel(cov, 
                         "Elevation (m)", 
                         "NPP", 
                         "% Tree Cover", 
                         "% Non-vegetated Area",
                         "Population Density (HII)", 
                         "Road Density (km/km2)"))


#horizontal facet
ggplot(plot_df) +
  geom_boxplot(aes(x = class, y = value, fill=class, col="grey")) +
  scale_fill_manual(values = viridis::magma(10), na.value = NA, guide = guide_legend(reverse = FALSE), na.translate=FALSE) +
  scale_color_manual(values = "#6C6B66") +
  facet_wrap(vars(cov), scales = "free") +
  geom_hline(aes(yintercept = median, group = cov), colour = 'red') +
  labs(x = "iSSF Class", y = "Landscape Parameter Value") +
  theme(legend.position = "none")

#ggsave("ssf_sim_100m_max_cap_habitat_boxplot_3-12-24.png", width = 10, height = 5, units = "in")

#vertical facet
ggplot(plot_df) +
  geom_boxplot(aes(x = class, y = value, fill=class, col="grey")) +
  scale_fill_manual(values = viridis::magma(10), na.value = NA, guide = guide_legend(reverse = FALSE), na.translate=FALSE) +
  scale_color_manual(values = "#6C6B66") +
  facet_wrap(vars(cov), scales = "free", dir="v") +
  geom_hline(aes(yintercept = median, group = cov), colour = 'red') +
  labs(x = "iSSF Class", y = "Landscape Parameter Value") +
  theme(legend.position = "none")

#ggsave("ssf_sim_100m_max_cap_habitat_boxplot_vertical_3-12-24.png", width = 10, height = 8, units = "in")


#########################################################################
##
##  3. Model coefficient plots and tables
##
##########################################################################

coeffs <- top_mod$summary.fixed %>% 
  rownames_to_column("Model term") %>% 
  as_tibble() %>% 
  select(-c(mode, kld)) %>% 
  rename(Mean = mean,
         `Std. Dev` = sd,
         `2.5%` = `0.025quant`,
         `50%` = `0.5quant`,
         `97.5%` = `0.975quant`) %>% 
  mutate(`Model term` = fct_recode(
    `Model term`,
    "Elevation" = "elevation",
    "NPP" = "npp_annual",
    "% Non-vegetated area" = "perc_nonveg_annual",
    "% Tree cover" = "perc_tree_cov_annual",
    "(% Tree cover)^2" = "I(perc_tree_cov_annual^2)",
    "Population density" = "popdens_hii_annual",
    "Road density" = "dens_all_roads_annual",
    "Step length" = "sl_",
    "ln(Step length)" = "log_sl_",
    "cos(Turn angle)" = "cos_ta_"
  )) %>% 
  mutate(`Model term`=fct_relevel(`Model term`, 
                                  "Elevation", 
                                  "NPP", 
                                  "% Tree cover", 
                                  "(% Tree cover)^2",
                                  "% Non-vegetated area",
                                  "Population density", 
                                  "Road density",
                                  "Step length",
                                  "ln(Step length)",
                                  "cos(Turn angle)")) %>% 
  arrange(`Model term`)

coeffs_exp <- coeffs %>% 
  mutate(across(-`Model term`, exp))


#Tables

coeff_tab_raw <- flextable(coeffs) %>% 
  colformat_double(digits = 4) %>% 
  width(width=1.7)

coeff_tab_exp <- flextable(coeffs_exp) %>% 
  colformat_double(digits = 4) %>% 
  width(width=1.7)

# save_as_docx("raw coefficients" = coeff_tab_raw,
#              "exponentiated coefficients" = coeff_tab_exp,
#              path = "top_ssf_fixed_effect_coeff_tables_2-22-24.docx")


#Plots 

#raw coeffs
ggplot(coeffs %>% mutate(`Model term` = fct_rev(`Model term`))) +
  geom_point(aes(x = Mean, y=`Model term`)) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`, y=`Model term`), height=0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Mean \u03B2 estimate for fixed effects")

#exponentiated coeffs
ggplot(coeffs_exp %>% mutate(`Model term` = fct_rev(`Model term`))) +
  geom_point(aes(x = Mean, y=`Model term`)) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`, y=`Model term`), height=0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "Mean \u03B2 estimate for fixed effects")




