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


##### Setup --------------------------------------------
  library(ggplot2)
  library(dplyr)

# Load data
  e_covs <- read.csv("data/enes.prepost.multiscale.occu.csv") 
  o_covs <- read.csv("data/oss.prepost.multiscale.occu.csv")
  load("data/msc-enes-data-workspace.RData")
  load("data/multiscale_output_and_data_072125_enes_small.RData")
  E = a2
  load("data/multiscale_output_and_data_072525_oss_small.RData")
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
    year = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9),
    treatment = c("UU", "HU", "UU", "HU", "UU", "HU", "UU", "HU", 
                  "UU", "HU", "UU", "HU", "UU", "HU", # 2013-2019, just UU, HU
                  "UU", "HU", "BU", "HB", "BS", # 2023, all trts
                  "UU", "HU", "BU", "HB", "BS")) # 2024, all trts
  
  
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
  BU_psi_CIs <- apply(BU_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)
  
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
  HB_psi_CIs <- apply(HB_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)
  
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
  HU_psi_CIs <- apply(HU_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)
  
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
  BS_psi_CIs <- apply(BS_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)
  
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
  UU_psi_CIs <- apply(UU_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)
  
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
  
  
# Plot with points and lines, automatically handling missing combinations
  p2 <- ggplot(year_treatment_preds, aes(x = year, y = predicted, color = treatment)) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.1, alpha = 0.7) +
    labs(x = "Year", 
         y = "Predicted Occupancy Probability",
         title = "Occupancy Estimates by Year and Treatment - OSS") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks = unique(year_treatment_preds$year))
  
  
  ggsave("figures/o-yearly-trt-preds.png", plot = p2, dpi = 300)




