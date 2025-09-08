## =================================================
##
## Title: msc-occu-posterior-probs-for-pub
## Author: Jasmine Williamson 
## Date Created: 09/08/2025
##
## Description: Assess estimates for publication results section
##
## =================================================
setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/multi-scale-occu-oss")


## Setup 
  library(ggplot2)
  library(dplyr)

# Load data
  e_covs <- read.csv("data/enes.prepost.multiscale.occu.csv") 
  o_covs <- read.csv("data/oss.prepost.multiscale.occu.csv")
  #load("data/msc-enes-data-workspace.RData")
  load("data/multiscale_output_and_data_082525_enes_small.RData")
  E = a2
  load("data/multiscale_output_and_data_082525_oss_full.RData")
  O = a

# Combine
  summary(E)
  summary(O)
  E2 = runjags::combine.mcmc(E)
  O2 = runjags::combine.mcmc(O)  
  
  
#### OSS - calculate posterior probabilities from MCMC samples ------------------
  
  
  # Extract posterior samples for each treatment effect from combined MCMC object
  samples_HU <- O2[, "beta3.psi.HU"]  # harvest only
  samples_BS <- O2[, "beta4.psi.BS"]  # burn-salvage  
  samples_HB <- O2[, "beta2.psi.HB"]  # harvest-burn
  samples_BU <- O2[, "beta1.psi.BU"]  # burn only
  
  
  # Calculate probabilities that treatments are LOWER than controls
  # (negative coefficients mean lower occupancy than baseline)
  
  # Probability that harvest-only is lower than controls
  prob_HU_lower <- mean(samples_HU < 0)
  cat("Probability HU < Control:", prob_HU_lower, "\n")
  
  mean(samples_HU)
  table(samples_HU<0)
  
  # Probability that burn-salvage is lower than controls  
  prob_BS_lower <- mean(samples_BS < 0)
  cat("Probability BS < Control:", prob_BS_lower, "\n")
  
  # Probability that harvest-burn is lower than controls
  prob_HB_lower <- mean(samples_HB < 0)
  cat("Probability HB < Control:", prob_HB_lower, "\n")
  
  # Probability that burn-only is lower than controls
  prob_BU_lower <- mean(samples_BU < 0)
  cat("Probability BU < Control:", prob_BU_lower, "\n")
  
  
  # Compare treatments to each other
  
  # Probability that burn-salvage is worse than harvest-only
  prob_BS_worse_than_HU <- mean(samples_BS < samples_HU)
  cat("Probability BS < HU:", prob_BS_worse_than_HU, "\n")
  
  # Probability that burn-salvage is worse than harvest-burn
  prob_BS_worse_than_HB <- mean(samples_BS < samples_HB)
  cat("Probability BS < HB:", prob_BS_worse_than_HB, "\n")
  
  
  # Convert to percentages for reporting
  cat("\nFor reporting:\n")
  cat("Harvest-only lower than controls:", round(prob_HU_lower * 100, 1), "%\n")
  cat("Burn-salvage lower than controls:", round(prob_BS_lower * 100, 1), "%\n")
  cat("Burn-salvage worse than harvest-only:", round(prob_BS_worse_than_HU * 100, 1), "%\n")
  
  
  
  
  
#### ENES - calculate posterior probabilities from MCMC samples ------------------
  
  
  # Extract posterior samples for each treatment effect from combined MCMC object
  samples_HU_e <- E2[, "beta3.psi.HU"]  # harvest only
  samples_BS_e <- E2[, "beta4.psi.BS"]  # burn-salvage  
  samples_HB_e <- E2[, "beta2.psi.HB"]  # harvest-burn
  samples_BU_e <- E2[, "beta1.psi.BU"]  # burn only
  
  
  # Calculate probabilities that treatments are LOWER than controls
  # (negative coefficients mean lower occupancy than baseline)
  
  # Probability that harvest-only is lower than controls
  prob_HU_lower <- mean(samples_HU_e < 0)
  cat("Probability HU < Control:", prob_HU_lower, "\n")
  
  # Probability that burn-salvage is lower than controls  
  prob_BS_lower <- mean(samples_BS_e < 0)
  cat("Probability BS < Control:", prob_BS_lower, "\n")
  
  # Probability that harvest-burn is lower than controls
  prob_HB_lower <- mean(samples_HB_e < 0)
  cat("Probability HB < Control:", prob_HB_lower, "\n")
  
  # Probability that burn-only is lower than controls
  prob_BU_lower <- mean(samples_BU_e < 0)
  cat("Probability BU < Control:", prob_BU_lower, "\n")
  
  
  # Compare treatments to each other
  
  # Probability that burn-salvage is worse than harvest-only
  prob_BS_worse_than_HU <- mean(samples_BS_e < samples_HU)
  cat("Probability BS < HU:", prob_BS_worse_than_HU, "\n")
  
  # Probability that burn-salvage is worse than harvest-burn
  prob_BS_worse_than_HB <- mean(samples_BS_e < samples_HB)
  cat("Probability BS < HB:", prob_BS_worse_than_HB, "\n")
  
  
  # Convert to percentages for reporting
  cat("\nFor reporting:\n")
  cat("Harvest-only lower than controls:", round(prob_HU_lower * 100, 1), "%\n")
  cat("Burn-salvage lower than controls:", round(prob_BS_lower * 100, 1), "%\n")
  cat("Burn-salvage worse than harvest-only:", round(prob_BS_worse_than_HU * 100, 1), "%\n")
  
  # posterior density plots 
  
  
  