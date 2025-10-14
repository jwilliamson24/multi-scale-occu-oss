## =================================================
##
## Title: msc-occu-yearly-plots
## Author: Jasmine Williamson 
## Date Created: 08/04/2025
##
## Description: Make plots that show yearly occu estimates
##
## =================================================
setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/multi-scale-occu-oss")


## Setup 
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(cowplot)

# Load data
  e_covs <- read.csv("data/enes.prepost.multiscale.occu.csv") 
  o_covs <- read.csv("data/oss.prepost.multiscale.occu.csv")
  load("data/msc-enes-data-workspace-V2.RData")
  load("data/multiscale_output_082525_enes_small.RData")
  E = a2
  load("data/multiscale_output_082525_oss_small.RData")
  O = a2
  

# Combine
  summary(E)
  summary(O)
  E2 = runjags::combine.mcmc(E)
  O2 = runjags::combine.mcmc(O)   


#### Plot predicted occupancy for each year and treatment 
  
  b <- O2  # # # # choose species # # # #
  
  # number of posterior samples
  n.samples = nrow(b)
  
  # dataframe with only valid year-treatment combos
  valid_combos <- data.frame(
    year = c(1, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9),
    treatment = c("UU", "UU", "UU", "HU", "UU", "HU", 
                  "UU", "HU", "UU", "HU", "UU", "HU", # 2013-2019, just UU, HU
                  "UU", "HU", "BU", "HB", "BS", # 2023, all trts
                  "UU", "HU", "BU", "HB", "BS")) # 2024, all trts
  
  #CI_range <- c(0.05,0.95) # 90%
  CI_range <- c(0.025, 0.975) # 95%
  #CI_range <- c(0.075, 0.925) # 85%

  
  
# BU    
  # get years where BU treatment occurred
  BU_years <- valid_combos[valid_combos$treatment == "BU", "year"]
  BU_data <- BU_years
  
  # set all other covs at mean
  HB = 0
  HU = 0
  BS = 0
  lat = 0
  long = 0
  elev = 0
  
  # create matrices to stick estimates in
  logit_psi = matrix(NA, n.samples, length(BU_data))
  
  # psi predictions for BU
  # Sample from posterior for the sequence of values for BU - raw linear predictor for each MCMC sample 
  for (i in 1:n.samples){
    for (j in 1:length(BU_data)){
      # create the linear predictors for dominant species
      year_effect_col <- paste0("beta0.psi.year[", BU_data[j], "]")
      logit_psi[i,j] = b[,year_effect_col][[i]] + b[,'beta1.psi.BU'][[i]] * 1 + 
        b[,'beta2.psi.HB'][[i]] * HB + b[,'beta3.psi.HU'][[i]] * HU +
        b[,'beta4.psi.BS'][[i]] * BS + b[,'beta5.psi.lat'][[i]] * lat +  
        b[,'beta6.psi.lon'][[i]] * long +  b[,'beta8.psi.elev'][[i]] * elev 
      }}

  # create array 
  BU_psi = matrix(NA, n.samples, length(BU_data))
  
  # transform psi off logit-scale back to probability scale - for each sample
  BU_psi <- plogis(logit_psi)
  
  # calculate means and credible intervals - mean predicted occupancy probability across all posterior samples
  BU_psi_means = colMeans(BU_psi) # - taking the mean of all samples (transformed occu prob)
  BU_psi_CIs <- apply(BU_psi,2,quantile, c(CI_range), na.rm=TRUE)
  
  # stuff into df
  BU_psi_preds <- data.frame(year = BU_data, 
                             predicted = BU_psi_means, 
                             treatment = "BU",
                             LCI = BU_psi_CIs[1,],
                             UCI = BU_psi_CIs[2,]) 
  
  
# HB
  # get years where HB treatment occurred
  HB_years <- valid_combos[valid_combos$treatment == "HB", "year"]
  HB_data <- HB_years
  
  # set all other covs at mean
  BU = 0
  HU = 0
  BS = 0
  lat = 0
  long = 0
  elev = 0
  
  # create matrices to stick estimates in
  logit_psi = matrix(NA, n.samples, length(HB_data))
  
  # psi predictions for HB
  # Sample from posterior for the sequence of values for HB
  for (i in 1:n.samples){
    for (j in 1:length(HB_data)){
      # create the linear predictors for dominant species
      year_effect_col <- paste0("beta0.psi.year[", HB_data[j], "]")
      logit_psi[i,j] = b[,year_effect_col][[i]] + b[,'beta1.psi.BU'][[i]] * BU + 
        b[,'beta2.psi.HB'][[i]] * 1 + b[,'beta3.psi.HU'][[i]] * HU +
        b[,'beta4.psi.BS'][[i]] * BS + b[,'beta5.psi.lat'][[i]] * lat +  
        b[,'beta6.psi.lon'][[i]] * long +  b[,'beta8.psi.elev'][[i]] * elev 
    }}
  
  HB_psi = matrix(NA, n.samples, length(HB_data)) # create array 
  HB_psi <- plogis(logit_psi) # transform psi off logit-scale back to probability scale
  HB_psi_means = colMeans(HB_psi) # calculate means and credible intervals
  HB_psi_CIs <- apply(HB_psi,2,quantile, c(0.05,0.95), na.rm=TRUE)
  
  HB_psi_preds <- data.frame(year = HB_data, 
                             predicted = HB_psi_means, 
                             treatment = "HB",
                             LCI = HB_psi_CIs[1,],
                             UCI = HB_psi_CIs[2,])
  
  
# HU   
  # get years where HU treatment occurred
  HU_years <- valid_combos[valid_combos$treatment == "HU", "year"]
  HU_data <- HU_years
  
  # set all other covs at mean
  BU = 0
  HB = 0
  BS = 0
  lat = 0
  long = 0
  elev = 0
  
  # create matrices to stick estimates in
  logit_psi = matrix(NA, n.samples, length(HU_data))
  
  # psi predictions for HU
  # Sample from posterior for the sequence of values for HB
  for (i in 1:n.samples){
    for (j in 1:length(HU_data)){
      # create the linear predictors for dominant species
      year_effect_col <- paste0("beta0.psi.year[", HU_data[j], "]")
      logit_psi[i,j] = b[,year_effect_col][[i]] + b[,'beta1.psi.BU'][[i]] * BU + 
        b[,'beta2.psi.HB'][[i]] * HB + b[,'beta3.psi.HU'][[i]] * 1 +
        b[,'beta4.psi.BS'][[i]] * BS + b[,'beta5.psi.lat'][[i]] * lat +  
        b[,'beta6.psi.lon'][[i]] * long +  b[,'beta8.psi.elev'][[i]] * elev 
    }}
  
  HU_psi = matrix(NA, n.samples, length(HU_data))   # create array 
  HU_psi <- plogis(logit_psi)   # transform psi off logit-scale back to probability scale
  HU_psi_means = colMeans(HU_psi)   # calculate means and credible intervals
  HU_psi_CIs <- apply(HU_psi,2,quantile, c(CI_range), na.rm=TRUE)
  
  HU_psi_preds <- data.frame(year = HU_data, 
                             predicted = HU_psi_means, 
                             treatment = "HU",
                             LCI = HU_psi_CIs[1,],
                             UCI = HU_psi_CIs[2,])  
  
  
# BS
  # Get years where BS treatment occurred
  BS_years <- valid_combos[valid_combos$treatment == "BS", "year"]
  BS_data <- BS_years
  
  # set all other covs at mean
  BU = 0
  HB = 0
  HU = 0
  lat = 0
  long = 0
  elev = 0
  
  # create matrices to stick estimates in
  logit_psi = matrix(NA, n.samples, length(BS_data))
  
  # psi predictions for BS
  # Sample from posterior for the sequence of values for BS 
  for (i in 1:n.samples){
    for (j in 1:length(BS_data)){
      # create the linear predictors for dominant species
      year_effect_col <- paste0("beta0.psi.year[", BS_data[j], "]")
      logit_psi[i,j] = b[,year_effect_col][[i]] + b[,'beta1.psi.BU'][[i]] * BU + 
        b[,'beta2.psi.HB'][[i]] * HB + b[,'beta3.psi.HU'][[i]] * HU +
        b[,'beta4.psi.BS'][[i]] * 1 + b[,'beta5.psi.lat'][[i]] * lat +  
        b[,'beta6.psi.lon'][[i]] * long +  b[,'beta8.psi.elev'][[i]] * elev 
    }}
  
  BS_psi = matrix(NA, n.samples, length(BS_data))   # create array 
  BS_psi <- plogis(logit_psi)   # transform psi off logit-scale back to probability scale
  BS_psi_means = colMeans(BS_psi)  # calculate means and credible intervals
  BS_psi_CIs <- apply(BS_psi,2,quantile, c(CI_range), na.rm=TRUE)
  
  BS_psi_preds <- data.frame(year = BS_data, 
                             predicted = BS_psi_means, 
                             treatment = "BS",
                             LCI = BS_psi_CIs[1,],
                             UCI = BS_psi_CIs[2,])
  
  
# UU  
  # Get years where UU (untreated) occurred
  UU_years <- valid_combos[valid_combos$treatment == "UU", "year"]
  UU_data <- UU_years
  
  # set all other covs at mean (all treatments = 0 for untreated)
  BU = 0
  HB = 0
  HU = 0
  BS = 0
  lat = 0
  long = 0
  elev = 0
  
  # create matrices to stick estimates in
  logit_psi = matrix(NA, n.samples, length(UU_data))
  
  # psi predictions for UU
  for (i in 1:n.samples){
    for (j in 1:length(UU_data)){
      # create the linear predictors for dominant species
      year_effect_col <- paste0("beta0.psi.year[", UU_data[j], "]")
      logit_psi[i,j] = b[,year_effect_col][[i]] + b[,'beta1.psi.BU'][[i]] * BU + 
        b[,'beta2.psi.HB'][[i]] * HB + b[,'beta3.psi.HU'][[i]] * HU +
        b[,'beta4.psi.BS'][[i]] * BS + b[,'beta5.psi.lat'][[i]] * lat +  
        b[,'beta6.psi.lon'][[i]] * long +  b[,'beta8.psi.elev'][[i]] * elev 
    }}
  
  UU_psi = matrix(NA, n.samples, length(UU_data))   # create array 
  UU_psi <- plogis(logit_psi)   # transform psi off logit-scale back to probability scale
  UU_psi_means = colMeans(UU_psi)  # calculate means and credible intervals
  UU_psi_CIs <- apply(UU_psi,2,quantile, c(CI_range), na.rm=TRUE)
  
  UU_psi_preds <- data.frame(year = UU_data, 
                             predicted = UU_psi_means, 
                             treatment = "UU",
                             LCI = UU_psi_CIs[1,],
                             UCI = UU_psi_CIs[2,])  
  
  
  
# Combine - only include dataframes for treatments that have data
  all_treatment_preds <- list()
  
  if(nrow(BU_psi_preds) > 0) all_treatment_preds[["BU"]] <- BU_psi_preds
  if(nrow(HB_psi_preds) > 0) all_treatment_preds[["HB"]] <- HB_psi_preds  
  if(nrow(HU_psi_preds) > 0) all_treatment_preds[["HU"]] <- HU_psi_preds
  if(nrow(BS_psi_preds) > 0) all_treatment_preds[["BS"]] <- BS_psi_preds
  if(nrow(UU_psi_preds) > 0) all_treatment_preds[["UU"]] <- UU_psi_preds
  
  year_treatment_preds <- do.call(rbind, all_treatment_preds)  
  print(year_treatment_preds)
  row.names(year_treatment_preds) <- NULL
  
  
  # save as you go to have one for each spp
  # year_treatment_preds_e <- year_treatment_preds
  # year_treatment_preds_o <- year_treatment_preds
  
  
  

#### 2023-2024 Predictions for paper ----------------------------------------
  
  # Average predictions across years 8 and 9 for each treatment
  hb_avg <- year_treatment_preds_e[year_treatment_preds_e$treatment == "HB" & year_treatment_preds_e$year %in% c(8, 9), ]
  hu_avg <- year_treatment_preds_e[year_treatment_preds_e$treatment == "HU" & year_treatment_preds_e$year %in% c(8, 9), ]
  bu_avg <- year_treatment_preds_e[year_treatment_preds_e$treatment == "BU" & year_treatment_preds_e$year %in% c(8, 9), ]
  uu_avg <- year_treatment_preds_e[year_treatment_preds_e$treatment == "UU" & year_treatment_preds_e$year %in% c(8, 9), ]
  
  # Calculate average predictions
  hb_pred_avg <- mean(hb_avg$predicted)
  hu_pred_avg <- mean(hu_avg$predicted)
  bu_pred_avg <- mean(bu_avg$predicted)
  uu_pred_avg <- mean(uu_avg$predicted)
  
  # Calculate average CI widths (for uncertainty estimation)
  hb_ci_width <- mean(hb_avg$UCI - hb_avg$LCI)
  hu_ci_width <- mean(hu_avg$UCI - hu_avg$LCI)
  bu_ci_width <- mean(bu_avg$UCI - bu_avg$LCI)  
  uu_ci_width <- mean(uu_avg$UCI - uu_avg$LCI)    
#### Plot --------------------------------------------------------------------- 

    p2 <- ggplot(year_treatment_preds, aes(x = year, y = predicted, color = treatment,
                                         shape = treatment)) +
    geom_pointrange(aes(ymin = LCI, ymax = UCI), 
                    position = position_dodge(width = 0.6),
                    size = 0.8) +
    labs(x = "Year", 
         y = "Predicted Occupancy Probability",
         title = "Occupancy Estimates by Year and Treatment - OSS") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = "bottom",
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()) +
    scale_x_continuous(breaks = unique(year_treatment_preds$year))
  
  
  #ggsave("figures/o-yearly-preds-shape.png", plot = p2, dpi = 300, bg = "white")
  

  
  
  
#### Plot with custom dodge widths *Claude help---------------

#ENES
  
  # Set the order of treatments
  year_treatment_preds_e$treatment <- factor(year_treatment_preds_e$treatment, 
                                           levels = c("UU", "BU", "HU", "HB", "BS"))
  
  # Version with custom spacing parameters for easy adjustment
  compress_factor <- 0.62  # Smaller = more compressed (try 0.3 to 0.6)
  year_8_position <- 5.5    # Position of year 8
  year_9_position <- 6.8    # Position of year 9
  dodge_1_7 <- 0.4        # Dodge width for years 1-7
  dodge_8_9 <- 0.9        # Dodge width for years 8-9
  
  year_treatment_preds_e$x_custom <- case_when(
    year_treatment_preds_e$year <= 7 ~ year_treatment_preds_e$year * compress_factor,
    year_treatment_preds_e$year == 8 ~ year_8_position,
    year_treatment_preds_e$year == 9 ~ year_9_position
  )
  
  data_1_7_custom_e <- year_treatment_preds_e[year_treatment_preds_e$year <= 7, ]
  data_8_9_custom_e <- year_treatment_preds_e[year_treatment_preds_e$year >= 8, ]
  
  p2_adjustable_e <- ggplot() +
    geom_pointrange(data = data_1_7_custom_e,
                    aes(x = x_custom, y = predicted, color = treatment,
                        ymin = LCI, ymax = UCI),
                    position = position_dodge(width = dodge_1_7),
                    size = 0.5) +
    geom_pointrange(data = data_8_9_custom_e,
                    aes(x = x_custom, y = predicted, color = treatment,
                        ymin = LCI, ymax = UCI),
                    position = position_dodge(width = dodge_8_9),
                    size = 0.5) +
    scale_color_manual(name = "Treatment",
                      values = c("BS" = "#EA7317", 
                                  "HB" = "#FEC601",
                                  "HU" = "#62B6CB", 
                                  "UU" = "#2364AA",
                                  "BU" = "#69995D"),
                       labels = c("BS" = "Burn-Salvage",
                                  "BU" = "Burn-only", 
                                  "HB" = "Harvest-Burn", 
                                  "HU" = "Harvest-only",
                                  "UU" = "Control")) +
    labs(x = NULL, 
         y = "Predicted "*psi~"",
         title = expression("Occupancy Estimates by Year and Treatment - "*italic("E. eschscholtzii"))) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 16, margin = margin(t = 10)),
          axis.title.y = element_text(size = 16, margin = margin(r = 10)),
          legend.position = "right",
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()) +
    scale_x_continuous(breaks = c(1:7 * compress_factor, year_8_position, year_9_position),
                       labels = c(2013, 2014, 2015, 2016, 2017, 2018, 2019, 2023, 2024),
                       limits = c(0.6, 7.5)) +  
    geom_vline(xintercept = (7 * compress_factor + year_8_position) / 2-0.2, 
               linetype = "dotted", alpha = 0.6) +
    geom_vline(xintercept = (2 * compress_factor + 3 * compress_factor) / 2,  # Between years 2-3
               linetype = "dotted", alpha = 0.6)
  
  p2_adjustable_e
   
  # ggsave("figures/e-yearly-preds-v2.png", 
  #        plot = p2_adjustable_e, 
  #        width = 10,
  #        height = 4,
  #        dpi = 300, 
  #        bg = "white")
  

    
# OSS
  
  # Set the order of treatments
  year_treatment_preds_o$treatment <- factor(year_treatment_preds_o$treatment, 
                                             levels = c("UU", "BU", "HU", "HB", "BS"))
  
  # Version with custom spacing parameters for easy adjustment
  compress_factor <- 0.62  # Smaller = more compressed (try 0.3 to 0.6)
  year_8_position <- 5.5    # Position of year 8
  year_9_position <- 6.8    # Position of year 9
  dodge_1_7 <- 0.4        # Dodge width for years 1-7
  dodge_8_9 <- 0.9        # Dodge width for years 8-9
  
  year_treatment_preds_o$x_custom <- case_when(
    year_treatment_preds_o$year <= 7 ~ year_treatment_preds_o$year * compress_factor,
    year_treatment_preds_o$year == 8 ~ year_8_position,
    year_treatment_preds_o$year == 9 ~ year_9_position
  )
  
  data_1_7_custom_o <- year_treatment_preds_o[year_treatment_preds_o$year <= 7, ]
  data_8_9_custom_o <- year_treatment_preds_o[year_treatment_preds_o$year >= 8, ]
  
  p2_adjustable_o <- ggplot() +
    geom_pointrange(data = data_1_7_custom_o,
                    aes(x = x_custom, y = predicted, color = treatment,
                        ymin = LCI, ymax = UCI),
                    position = position_dodge(width = dodge_1_7),
                    size = 0.5) +
    geom_pointrange(data = data_8_9_custom_o,
                    aes(x = x_custom, y = predicted, color = treatment,
                        ymin = LCI, ymax = UCI),
                    position = position_dodge(width = dodge_8_9),
                    size = 0.5) +
    scale_color_manual(name = "Treatment",
                       values = c("BS" = "#EA7317", 
                                  "HB" = "#FEC601",
                                  "HU" = "#62B6CB", 
                                  "UU" = "#2364AA",
                                  "BU" = "#69995D"),
                       labels = c("BS" = "Salvage",
                                  "BU" = "Burn-only", 
                                  "HB" = "Harvest-Burn", 
                                  "HU" = "Harvest-only",
                                  "UU" = "Control")) +
    labs(x = NULL, 
         y = "Predicted "*psi~"",
         title = expression("Occupancy Estimates by Year and Treatment - "*italic("B. wrighti"))) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 16, margin = margin(t = 10)),
          axis.title.y = element_text(size = 16, margin = margin(r = 10)),
          legend.position = "right",
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()) +
    scale_x_continuous(breaks = c(1:7 * compress_factor, year_8_position, year_9_position),
                       labels = c(2013, 2014, 2015, 2016, 2017, 2018, 2019, 2023, 2024),
                       limits = c(0.6, 7.5)) +  
    geom_vline(xintercept = (7 * compress_factor + year_8_position) / 2-0.2, 
               linetype = "dotted", alpha = 0.6) +
    geom_vline(xintercept = (2 * compress_factor + 3 * compress_factor) / 2,  # Between years 2-3
               linetype = "dotted", alpha = 0.6)
  
  p2_adjustable_o
  
  # ggsave("figures/o-yearly-preds-v2.png", 
  #        plot = p2_adjustable_o, 
  #        width = 10,
  #        height = 4,
  #        dpi = 300, 
  #        bg = "white")

  
  
#### Stacked plot *ChatGPT help -------------------------------------------------
  
  # spacing preferences
  compress_factor <- 0.6  # Smaller = more compressed (try 0.3 to 0.6)
  year_8_position <- 5.5    # Position of year 8
  year_9_position <- 6.8    # Position of year 9
  dodge_1_7 <- 0.3        # Dodge width for years 1-7
  dodge_8_9 <- 0.9        # Dodge width for years 8-9
  
  
  # Combine the two prediction dataframes into one
  # Make sure these dataframes already have: year, treatment, predicted, LCI, UCI, x_custom
  both_preds <- bind_rows(
    year_treatment_preds_o %>% mutate(species = "B. wrighti"),
    year_treatment_preds_e %>% mutate(species = "E. eschscholtzii")
  )
  head(both_preds)
  
  # separate years 1-7 and 8-9
  data_1_7 <- both_preds %>% filter(year <= 7)
  data_8_9 <- both_preds %>% filter(year >= 8)
  
  
  # order trt levels
  both_preds$treatment <- factor(both_preds$treatment,
                                 levels = c("UU", "BU", "HU", "HB", "BS"))
  
  # Shared color scale / labels
  treatment_colors <- c(
    "BS" = "#EA7317",
    "HB" = "#FEC601",
    "HU" = "#62B6CB",
    "UU" = "#2364AA",
    "BU" = "#69995D"
  )
  
  treatment_labels <- c(
    "BS" = "Salvage",
    "BU" = "Burn-only",
    "HB" = "Harvest-Burn",
    "HU" = "Harvest-only",
    "UU" = "Control"
  )
  
  
  p_facet <- ggplot(both_preds,
                    aes(x = x_custom, y = predicted,
                        color = treatment, group = treatment)) +
    # Years 1-7
      geom_pointrange(
        data = data_1_7,
        aes(x = x_custom, y = predicted, ymin = LCI, ymax = UCI, color = treatment, group = treatment),
        position = position_dodge(width = dodge_1_7),
        size = 0.55
      ) +
    # Years 8-9
      geom_pointrange(
        data = data_8_9,
        aes(x = x_custom, y = predicted, ymin = LCI, ymax = UCI, color = treatment, group = treatment),
        position = position_dodge(width = dodge_8_9),
        size = 0.55
      ) +
    facet_wrap(~ species, ncol = 1, scales = "fixed") +
    scale_color_manual(name = "Treatment",
                       values = treatment_colors,
                       labels = treatment_labels) +
    scale_x_continuous(
      breaks = c(1:7 * compress_factor, year_8_position, year_9_position),
      labels = c(2013, 2014, 2015, 2016, 2017, 2018, 2019, 2023, 2024),
      limits = c(0.6, 7.5)
    ) +
    # the same vertical-dotted markers you had before (they will line up across facets)
    geom_vline(xintercept = (7 * compress_factor + year_8_position) / 2 - 0.2,
               linetype = "dotted", alpha = 0.6) +
    geom_vline(xintercept = (2 * compress_factor + 3 * compress_factor) / 2,
               linetype = "dotted", alpha = 0.6) +
    labs(
      title = "Occupancy Estimates by Year and Treatment",
      y = expression("Predicted " * psi),
      x = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      #strip.text = element_text(size = 14, face = "italic"), # species as subtitle-like labels
      strip.text = element_blank(),
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 16, margin = margin(r = 10)),
      legend.position = "right",
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 14),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      #axis.line.x = element_line(color = "gray40"),
      #axis.line.y = element_line(color = "gray40"),
      panel.border = element_rect(color = "gray40", fill = NA, linewidth = 0.5),
      axis.ticks = element_line(color = "gray40", size = 0.5),
      axis.ticks.length = unit(0.1, "cm")
    ) 
  
  # Plot it
  p_facet
  
  
  ggsave("figures/facet-yearly-preds-v2.png",
         plot = p_facet,
         width = 10,
         height = 5,
         dpi = 300,
         bg = "white")



