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

# Estimates
  summary(E)
  summary(O)

# Combine
  E2 = runjags::combine.mcmc(E)
  O2 = runjags::combine.mcmc(O)   


##### Plot predicted occupancy just for each year - ignore this plot its not useful ---------------------
  
  dat <- as.matrix(E)  # # choose species # #
  
  # Extract all samples for each year intercept (beta0.psi.year[t])
  years <- 1:9
  year_psi_preds <- lapply(years, function(t) plogis(dat[, paste0("beta0.psi.year[", t, "]")]))
  year_means <- sapply(year_psi_preds, mean)
  year_CIs <- sapply(year_psi_preds, function(x) quantile(x, probs = c(0.025, 0.975)))
  
  # Put into data frame for ggplot
  e_year_preds_df <- data.frame(
    year = years,
    predicted = year_means,
    LCI = year_CIs[1, ],
    UCI = year_CIs[2, ]
  )
  
  # Plot
  ggplot(e_year_preds_df, aes(x = year, y = predicted)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2) +
    ylab(expression("Predicted "*psi*"")) +
    xlab("Year") +
    ggtitle("Predicted Occupancy by Year - ENES") +
    theme_classic()



##### Plot predicted occupancy for each year and treatment ---------------------
  
##### Method 1 ---------------------  
  
  dat <- as.matrix(E)
  
  treatments <- c("UU", "BU", "HB", "HU", "BS")
  years <- 1:9
  
  # Create a function to compute logit_psi for a given treatment and year
  get_psi_preds <- function(treatment, year_index) {
    # Extract posterior samples
    beta0 <- dat[, "beta0.psi"]
    psi_year <- dat[, paste0("beta0.psi.year[", year_index, "]")]
    
    # Set treatment coefficients (0 for control group)
    beta_BU <- ifelse(treatment == "BU", dat[, "beta1.psi.BU"], 0)
    beta_HB <- ifelse(treatment == "HB", dat[, "beta2.psi.HB"], 0)
    beta_HU <- ifelse(treatment == "HU", dat[, "beta3.psi.HU"], 0)
    beta_BS <- ifelse(treatment == "BS", dat[, "beta4.psi.BS"], 0)
    
    # Set covs to mean 
    beta_lat <- 0
    beta_lon <- 0
    beta_elev <- 0
    
    # Compute linear predictor and transform to probability
    logit_psi <- beta0 + psi_year + beta_BU + beta_HB + beta_HU + beta_BS
    psi <- plogis(logit_psi)
    
    # combine
    data.frame(
      treatment = treatment,
      year = year_index,
      mean_psi = mean(psi),
      lci = quantile(psi, 0.025),
      uci = quantile(psi, 0.975)
    )
  }
  
  
  # dataframe with only valid year-treatment combos
  valid_combos <- data.frame(
    year = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9),
    treatment = c("UU", "HU", "UU", "HU", "UU", "HU", "UU", "HU", 
                  "UU", "HU", "UU", "HU", "UU", "HU", 
                  "UU", "HU", "BU", "HB", "BS", "UU", "HU", "BU", "HB", "BS"))
  
  # predictions for each trt each year
  all_preds1 <- do.call(rbind, lapply(1:max(years), function(t) {
    # get treatments for this year
    trts_for_year <- valid_combos$treatment[valid_combos$year == t]
    
    do.call(rbind, lapply(trts_for_year, function(trt) get_psi_preds(trt, t)))
  }))
  
  row.names(all_preds1) <- NULL
  head(all_preds1)
  
  
  # plot
  ggplot(all_preds1, aes(x = factor(year), y = mean_psi, color = treatment, group = treatment)) +
    geom_point(position = position_dodge(0.3), size = 2.5) +
    geom_errorbar(aes(ymin = lci, ymax = uci), position = position_dodge(0.3), width = 0.2) +
    labs(x = "Year", y = "Predicted Occupancy (ψ)", title = "Occupancy by Treatment and Year - E") +
    theme_minimal()
  
  

  
  
### Method 2 --------------------  
  
  dat <- as.matrix(E)
  
  get_psi_preds <- function(treatment, year) {
    # Intercept for year
    beta0 <- dat[, "beta0.psi"] + dat[, paste0("beta0.psi.year[", year, "]")]
    
    # Treatment coefficients (UU is reference, so no added term)
    treatment_term <- switch(treatment,
                             "BU" = dat[, "beta1.psi.BU"],
                             "HB" = dat[, "beta2.psi.HB"],
                             "HU" = dat[, "beta3.psi.HU"],
                             "BS" = dat[, "beta4.psi.BS"],
                             0  # for "UU"
    )
    
    # Add other covariates, set to mean (0)
    cov_effect <- 0
    cov_effect <- cov_effect +
      0 * dat[, "beta8.psi.elev"] +
      0 * dat[, "beta5.psi.lat"] +
      0 * dat[, "beta6.psi.lon"]
    
    # transform
    logit_psi <- beta0 + treatment_term + cov_effect
    psi <- plogis(logit_psi)
    
    # combine
    data.frame(
      treatment = treatment,
      year = year,
      mean_psi = mean(psi),
      lci = quantile(psi, 0.025),
      uci = quantile(psi, 0.975)
    )
  }
  
  
  # dataframe with only valid year-treatment combos
  valid_combos <- data.frame(
    year = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9),
    treatment = c("UU", "HU", "UU", "HU", "UU", "HU", "UU", "HU", 
                  "UU", "HU", "UU", "HU", "UU", "HU", 
                  "UU", "HU", "BU", "HB", "BS", "UU", "HU", "BU", "HB", "BS"))
  
  
  # preds for each trt each year
  all_preds2 <- do.call(rbind, lapply(1:nrow(valid_combos), function(i) {
    get_psi_preds(valid_combos$treatment[i], valid_combos$year[i])
  }))
  
  
  # plot
  ggplot(all_preds2, aes(x = factor(year), y = mean_psi, color = treatment, group = treatment)) +
    geom_point(position = position_dodge(0.3), size = 3) +
    geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.1, position = position_dodge(0.3)) +
    labs(x = "Year", y = "Predicted occupancy (ψ)", color = "Treatment") +
    theme_bw(base_size = 14)
  
  
  
  
  
