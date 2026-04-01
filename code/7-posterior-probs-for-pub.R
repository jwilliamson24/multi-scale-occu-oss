## =================================================
##
## Title: posterior-estimates-faster-code
## Author: Jasmine Williamson 
## Date Created: 09/08/2025
##
## Description: Calculates both (1) posterior probabilities with CI and 
## (2) beta coefficient estimates with CI with table for publication
##
## These (1) are the same numbers calculated in script 4 for occupancy plot,
## but a faster version of the code. Help from Claude AI.
##
## =================================================
setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/multi-scale-occu-oss")


## posterior occupancy predictions and credible intervals go in-text to discuss
## direct changes in occupancy across treatments

## beta values go in a table
## can also be used in text to illustrate effect sizes; use wisely


#### --------------------------------------------------------------------------

# Setup 
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(nimble)
  library(flextable)

# Load data
  source("data/attach.nimble_v2 copy.R")
  load("data/msc-data-workspace-V3.RData")
  load("data/multiscale_output_and_data_121225_enes_full.RData")
  E = a
  load("data/multiscale_output_and_data_121225_oss_full.RData")
  O = a

# Combine
  summary(E)
  summary(O)
  E2 = runjags::combine.mcmc(E)
  O2 = runjags::combine.mcmc(O)  
  

#### Predicted occupancy for trts with uncertainty - ENES ----------------------

## back transforming parameters of interest
  
  # Get number of posterior samples
  n.samples <- nrow(E2)
  
  # Extract trt parameter columns
  beta0_samples <- E2[, "beta0.psi"]
  beta1_samples <- E2[, "beta1.psi.BU"]
  beta2_samples <- E2[, "beta2.psi.HB"]
  beta3_samples <- E2[, "beta3.psi.HU"]
  beta4_samples <- E2[, "beta4.psi.BS"]
  
  # Create matrix to stick estimates in: rows = posterior samples, columns = treatments
  logit_psi_preds_e <- matrix(NA, nrow = n.samples, ncol = 5)
  
  # Calculate linear predictors for each treatment across all posterior samples
  for(i in 1:n.samples){
    logit_psi_preds_e[i, 1] <- beta0_samples[i]                     # control
    logit_psi_preds_e[i, 2] <- beta0_samples[i] + beta1_samples[i]  # BU
    logit_psi_preds_e[i, 3] <- beta0_samples[i] + beta2_samples[i]  # HB
    logit_psi_preds_e[i, 4] <- beta0_samples[i] + beta3_samples[i]  # HU
    logit_psi_preds_e[i, 5] <- beta0_samples[i] + beta4_samples[i]  # BS
  }
  
  # Transform to probability scale
  psi_preds_e <- plogis(logit_psi_preds_e)
  
  # Get means and CIs
  psi_means_e <- colMeans(psi_preds_e)
  psi_CIs_e <- apply(psi_preds_e, 2, quantile, c(0.025, 0.5, 0.975))
  
  # Package into dataframe
  occupancy_summary_e <- data.frame(
    Treatment = c("Reference", "BU", "HB", "HU", "BS"),
    Mean = psi_means_e,
    Median = psi_CIs_e[2, ],
    CI_2.5 = psi_CIs_e[1, ],
    CI_97.5 = psi_CIs_e[3, ]
  )
  
  
  occupancy_summary_e[, 2:5] <- round(occupancy_summary_e[, 2:5], 3)
  print(occupancy_summary_e)


  # save
  #write.csv(occupancy_summary_e, "posterior-preds-e-V3.csv", row.names = FALSE)
  
  
  
#### Predicted occupancy for trts with uncertainty - OSS -----------------------  

  n.samples <- nrow(O2)
  
  beta0_samples <- O2[, "beta0.psi"]
  beta1_samples <- O2[, "beta1.psi.BU"]
  beta2_samples <- O2[, "beta2.psi.HB"]
  beta3_samples <- O2[, "beta3.psi.HU"]
  beta4_samples <- O2[, "beta4.psi.BS"]
  
  logit_psi_preds_o <- matrix(NA, nrow = n.samples, ncol = 5)
  
  for(i in 1:n.samples){
    logit_psi_preds_o[i, 1] <- beta0_samples[i]
    logit_psi_preds_o[i, 2] <- beta0_samples[i] + beta1_samples[i]
    logit_psi_preds_o[i, 3] <- beta0_samples[i] + beta2_samples[i]
    logit_psi_preds_o[i, 4] <- beta0_samples[i] + beta3_samples[i]
    logit_psi_preds_o[i, 5] <- beta0_samples[i] + beta4_samples[i]
  }
  
  psi_preds_o <- plogis(logit_psi_preds_o)
  
  psi_means_o <- colMeans(psi_preds_o)
  psi_CIs_o <- apply(psi_preds_o, 2, quantile, c(0.025, 0.5, 0.975))
  
  occupancy_summary_o <- data.frame(
    Treatment = c("Reference", "BU", "HB", "HU", "BS"),
    Mean = psi_means_o,
    Median = psi_CIs_o[2, ],
    CI_2.5 = psi_CIs_o[1, ],
    CI_97.5 = psi_CIs_o[3, ]
  )
  

  occupancy_summary_o[, 2:5] <- round(occupancy_summary_o[, 2:5], 3)
  print(occupancy_summary_o)
  
  # save
  #write.csv(occupancy_summary_o, "posterior-preds-o-V3.csv", row.names = FALSE)
  

#### Beta values (cov effects) for pub with CI -----------------
  
## just the non-back-transformed version of the above parameter values
## useful for reporting, not much else
  
## OSS 
  
  # Key parameters to report
  params <- c("beta0.psi", "beta1.psi.BU", "beta2.psi.HB", 
              "beta3.psi.HU", "beta4.psi.BS", "beta5.psi.lat", 
              "beta6.psi.lon", "beta8.psi.elev", "beta1.theta.DW",
              "alpha1", "alpha2")
  
  coef_summary_o <- data.frame(
    Parameter = params,
    Median = NA,
    CI_2.5 = NA,
    CI_97.5 = NA
  )
  
  for(i in 1:length(params)) {
    param_data <- O2[, params[i]]
    coef_summary_o[i, "Median"] <- median(param_data)
    coef_summary_o[i, "CI_2.5"] <- quantile(param_data, 0.025)
    coef_summary_o[i, "CI_97.5"] <- quantile(param_data, 0.975)
  }
  
  coef_summary_o[, 2:4] <- round(coef_summary_o[, 2:4], 3)
  print(coef_summary_o) 

  
## ENES 
  
  # Key parameters to report
  params <- c("beta0.psi", "beta1.psi.BU", "beta2.psi.HB", 
              "beta3.psi.HU", "beta4.psi.BS", "beta5.psi.lat", 
              "beta6.psi.lon", "beta8.psi.elev")
  
  coef_summary_e <- data.frame(
    Parameter = params,
    Median = NA,
    CI_2.5 = NA,
    CI_97.5 = NA
  )
  
  for(i in 1:length(params)) {
    param_data <- E2[, params[i]]
    coef_summary_e[i, "Median"] <- median(param_data)
    coef_summary_e[i, "CI_2.5"] <- quantile(param_data, 0.025)
    coef_summary_e[i, "CI_97.5"] <- quantile(param_data, 0.975)
  }
  
  coef_summary_e[, 2:4] <- round(coef_summary_e[, 2:4], 3)
  print(coef_summary_e) 
  
  
  

#### Create beta coeff flextable for pub ---------------------------------------

  
## OSS
  
  # Key parameters to report
  params <- c("beta0.psi", "beta1.psi.BU", "beta2.psi.HB", 
              "beta3.psi.HU", "beta4.psi.BS", "beta5.psi.lat", 
              "beta6.psi.lon", "beta8.psi.elev", "beta1.theta.DW",
              "alpha1", "alpha2")
  
  coef_summary_o <- data.frame(
    Parameter = params,
    Median = NA,
    CI_2.5 = NA,
    CI_97.5 = NA
  )
  
  for(i in 1:length(params)) {
    param_data <- O2[, params[i]]
    coef_summary_o[i, "Median"] <- median(param_data)
    coef_summary_o[i, "CI_2.5"] <- quantile(param_data, 0.025)
    coef_summary_o[i, "CI_97.5"] <- quantile(param_data, 0.975)
    coef_summary_o[i, "Percent_Pos"] <- mean(param_data > 0) * 100
  }
  
  coef_summary_o[, 2:4] <- round(coef_summary_o[, 2:4], 2)
  coef_summary_o$Percent_Pos <- round(coef_summary_o$Percent_Pos, 1)
  
  # Add submodel and covariate names
  coef_summary_o$Submodel <- c("ψ", "ψ", "ψ", "ψ", "ψ", "ψ", "ψ", "ψ", "θ", "p", "p")
  
  coef_summary_o$Covariate <- c("Intercept", "Burn-only", "Harvest-burn",
                                "Harvest-only", "Burn-salvage", "Latitude",
                                "Longitude", "Elevation", "Downed wood", 
                                "Temperature", "Temperature²")
  
## ENES
  
  # Key parameters to report
  params <- c("beta0.psi", "beta1.psi.BU", "beta2.psi.HB", 
              "beta3.psi.HU", "beta4.psi.BS", "beta5.psi.lat", 
              "beta6.psi.lon", "beta8.psi.elev", "beta1.theta.DW",
              "alpha1", "alpha2")
  
  coef_summary_e <- data.frame(
    Parameter = params,
    Median = NA,
    CI_2.5 = NA,
    CI_97.5 = NA
  )
  
  for(i in 1:length(params)) {
    param_data <- E2[, params[i]]
    coef_summary_e[i, "Median"] <- median(param_data)
    coef_summary_e[i, "CI_2.5"] <- quantile(param_data, 0.025)
    coef_summary_e[i, "CI_97.5"] <- quantile(param_data, 0.975)
    coef_summary_e[i, "Percent_Pos"] <- mean(param_data > 0) * 100
  }
  
  coef_summary_e[, 2:4] <- round(coef_summary_e[, 2:4], 2)
  coef_summary_e$Percent_Pos <- round(coef_summary_e$Percent_Pos, 1)
  
  # Add submodel and covariate names
  coef_summary_e$Submodel <- c("ψ", "ψ", "ψ", "ψ", "ψ", "ψ", "ψ", "ψ", "θ", "p", "p")
  
  coef_summary_e$Covariate <- c("Intercept", "Burn-only", "Harvest-burn",
                                "Harvest-only", "Burn-salvage", "Latitude",
                                "Longitude", "Elevation", "Downed wood", 
                                "Temperature", "Temperature²")
  
## COMBINE
  
  # Prepare OSS
  oss_table <- coef_summary_o %>%
    select(Submodel, Covariate, Median, CI_2.5, CI_97.5, Percent_Pos) %>%
    rename(OSS_Median = Median, OSS_CI_2.5 = CI_2.5, OSS_CI_97.5 = CI_97.5,
           OSS_Percent_Pos = Percent_Pos)
  
  # Prepare ENES
  enes_table <- coef_summary_e %>%
    select(Submodel, Covariate, Median, CI_2.5, CI_97.5, Percent_Pos) %>%
    rename(ENES_Median = Median, ENES_CI_2.5 = CI_2.5, ENES_CI_97.5 = CI_97.5,
           ENES_Percent_Pos = Percent_Pos)
  
  # Combine both species
  combined_table <- full_join(oss_table, enes_table, by = c("Submodel", "Covariate"))
  
  # Create nice CI columns with unique names
  combined_table <- combined_table %>%
    mutate(
      `OSS_Estimate` = ifelse(!is.na(OSS_Median), 
                              paste0(OSS_Median, " (", OSS_CI_2.5, ", ", OSS_CI_97.5, ")"),
                              "—"),
      `OSS_P` = ifelse(!is.na(OSS_Percent_Pos), 
                       paste0(OSS_Percent_Pos, "%"),
                       "—"),
      `ENES_Estimate` = ifelse(!is.na(ENES_Median),
                               paste0(ENES_Median, " (", ENES_CI_2.5, ", ", ENES_CI_97.5, ")"),
                               "—"),
      `ENES_P` = ifelse(!is.na(ENES_Percent_Pos),
                        paste0(ENES_Percent_Pos, "%"),
                        "—")
    ) %>%
    select(Submodel, Covariate, OSS_Estimate, OSS_P, ENES_Estimate, ENES_P)
  
  # Create flextable and rename cols
  ft <- flextable(combined_table) %>%
    set_header_labels(
      OSS_Estimate = "Estimate (95% CI)",
      OSS_P = "P(β>0)",
      ENES_Estimate = "Estimate (95% CI)",
      ENES_P = "P(β>0)"
    ) %>%
    merge_v(j = "Submodel") %>%
    theme_booktabs() %>%
    autofit() %>%
    align(j = 3:6, align = "center", part = "all") %>%
    align(j = 1, align = "left", part = "all") %>%
    bold(part = "header") %>%
    valign(j = 1, valign = "top") %>%
    add_header_row(values = c("", "", "Oregon Slender Salamander", "Ensatina"), 
                   colwidths = c(1, 1, 2, 2)) %>%
    align(i = 1, align = "center", part = "header") %>%
    bold(i = 1, part = "header")
  
  # View the table
  ft
  
  # Save to Word
  save_as_docx(ft, path = "beta_coefficients_table.docx")
  
    
#### Attach.nimble method ----------------------------------------------
  
  # Using attach.nimble function from Josh Stewarts Bayes class
  # Summary table: mean, med, CI, % pos/neg for each parameter 
  # Predicted occu prob point estimates (means), and percent change from control
  # NOT calculating posterior distributions here
  
  ## ENES
  
  # Convert mcmc.list to the format attach.nimble expects
  mcmc.output.e <- list()
  mcmc.output.e$chain1 <- as.matrix(E[[1]])
  mcmc.output.e$chain2 <- as.matrix(E[[2]])
  mcmc.output.e$chain3 <- as.matrix(E[[3]])
  
  attach.nimble(mcmc.output.e)
  
  hist(beta1.psi.BU)  
  median(beta1.psi.BU)
  mean(beta1.psi.BU<0)*100  # probability that BU is less than 0, aka negative effect
  apply(beta1.psi.BU,2,quantile, c(0.025, 0.5, 0.975)) # looking at CI
  
  
  # Table with estimates and CI
  
  params <- c("beta0.psi", "beta0.theta", "alpha0", "beta1.psi.BU", "beta2.psi.HB", 
              "beta3.psi.HU", "beta4.psi.BS", "beta5.psi.lat", "beta6.psi.lon", 
              "beta8.psi.elev", "beta1.theta.DW", "alpha1", "alpha2")
  
  # Create empty dataframe
  summary_table_e <- data.frame(
    Parameter = params,
    Mean = NA,
    Median = NA,
    CI_2.5 = NA,
    CI_97.5 = NA,
    Percent_Pos = NA,
    Percent_Neg = NA
  )
  
  # Use a for loop to fill in the values
  for(i in 1:length(params)) {
    param_name <- params[i]
    
    if(exists(param_name)) {     # Check if parameter exists in environment
      param_data <- get(param_name)
      
      if(is.matrix(param_data)) {       # Convert to vector if it's a matrix
        param_data <- as.vector(param_data)
      }
      
      # Calculate statistics
      summary_table_e[i, "Mean"] <- mean(param_data)
      summary_table_e[i, "Median"] <- median(param_data)
      summary_table_e[i, "CI_2.5"] <- quantile(param_data, 0.025)
      summary_table_e[i, "CI_97.5"] <- quantile(param_data, 0.975)
      summary_table_e[i, "Percent_Pos"] <- mean(param_data > 0) * 100
      summary_table_e[i, "Percent_Neg"] <- mean(param_data < 0) * 100
      
    } else {
      cat("Parameter", param_name, "not found\n")
    }
  }
  
  # Round the results
  summary_table_e$Mean <- round(summary_table_e$Mean, 3)
  summary_table_e$Median <- round(summary_table_e$Median, 3)
  summary_table_e$CI_2.5 <- round(summary_table_e$CI_2.5, 3)
  summary_table_e$CI_97.5 <- round(summary_table_e$CI_97.5, 3)
  summary_table_e$Percent_Pos <- round(summary_table_e$Percent_Pos, 2)
  summary_table_e$Percent_Neg <- round(summary_table_e$Percent_Neg, 2)
  
  print(summary_table_e)
  
  #write.csv(summary_table_e, "e_parameter_summary.csv", row.names = FALSE)
  
  
  # detection probability = 0.160
  plogis(-1.658) # back transform the mean of alpha 0
  
  
  
  
  ## OSS
  
  mcmc.output.o <- list()
  mcmc.output.o$chain1 <- as.matrix(O[[1]])
  mcmc.output.o$chain2 <- as.matrix(O[[2]])
  mcmc.output.o$chain3 <- as.matrix(O[[3]])
  
  attach.nimble(mcmc.output.o)
  
  # Table with estimates and CI
  
  params <- c("beta0.psi", "beta0.theta", "alpha0", "beta1.psi.BU", "beta2.psi.HB", 
              "beta3.psi.HU", "beta4.psi.BS", "beta5.psi.lat", "beta6.psi.lon", 
              "beta8.psi.elev", "beta1.theta.DW", "alpha1", "alpha2")
  
  summary_table_o <- data.frame(
    Parameter = params,
    Mean = NA,
    Median = NA,
    CI_2.5 = NA,
    CI_97.5 = NA,
    Percent_Pos = NA,
    Percent_Neg = NA
  )
  
  # Use a for loop to fill in the values
  for(i in 1:length(params)) {
    param_name <- params[i]
    
    if(exists(param_name)) {     # Check if parameter exists in environment
      param_data <- get(param_name)
      
      if(is.matrix(param_data)) {       # Convert to vector if it's a matrix
        param_data <- as.vector(param_data)
      }
      
      # Calculate statistics
      summary_table_o[i, "Mean"] <- mean(param_data)
      summary_table_o[i, "Median"] <- median(param_data)
      summary_table_o[i, "CI_2.5"] <- quantile(param_data, 0.025)
      summary_table_o[i, "CI_97.5"] <- quantile(param_data, 0.975)
      summary_table_o[i, "Percent_Pos"] <- mean(param_data > 0) * 100
      summary_table_o[i, "Percent_Neg"] <- mean(param_data < 0) * 100
      
    } else {
      cat("Parameter", param_name, "not found\n")
    }
  }
  
  # Round the results
  summary_table_o$Mean <- round(summary_table_o$Mean, 3)
  summary_table_o$Median <- round(summary_table_o$Median, 3)
  summary_table_o$CI_2.5 <- round(summary_table_o$CI_2.5, 3)
  summary_table_o$CI_97.5 <- round(summary_table_o$CI_97.5, 3)
  summary_table_o$Percent_Pos <- round(summary_table_o$Percent_Pos, 2)
  summary_table_o$Percent_Neg <- round(summary_table_o$Percent_Neg, 2)
  
  # View the table
  print(summary_table_o)
  
  #write.csv(summary_table_o, "o_parameter_summary.csv", row.names = FALSE)
  
  
  # detection probability = 0.248
  plogis(-1.107) # back transform the mean of alpha 0
  
  
  
  
  
  
  
  
  