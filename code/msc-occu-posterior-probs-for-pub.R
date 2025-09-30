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
  library(tibble)
  library(nimble)

# Load data
  source("data/attach.nimble_v2 copy.R")
  load("data/msc-enes-data-workspace-V2.RData")
  load("data/multiscale_output_and_data_082525_enes_full.RData")
  E = a
  load("data/multiscale_output_and_data_082525_oss_full.RData")
  O = a

# Combine
  summary(E)
  summary(O)
  E2 = runjags::combine.mcmc(E)
  O2 = runjags::combine.mcmc(O)  
  

#### Attach.nimble: mean, median, CI, percent pos/neg for each parameter ---------------

# Using attach.nimble function from Josh Stewarts Bayes class
  
## ENES
  
  # Convert mcmc.list to the format attach.nimble expects
  mcmc.output.e <- list()
  mcmc.output.e$chain1 <- as.matrix(E[[1]])
  mcmc.output.e$chain2 <- as.matrix(E[[2]])
  mcmc.output.e$chain3 <- as.matrix(E[[3]])
  
  attach.nimble(mcmc.output.e)

  hist(beta9.psi.dwd)  
  median(beta9.psi.dwd)
  mean(beta9.psi.dwd<0)*100
  apply(beta9.psi.dwd,2,quantile, c(0.025, 0.5, 0.975))
  
  # There is a significant negative effect of dwd on ENES occupancy
  # CI includes zero, so it isn't an incredibly strong effect
  

# Table with estimates and CI
  
  params <- c("beta0.psi", "beta0.theta", "alpha0", "beta1.psi.BU", "beta2.psi.HB", 
              "beta3.psi.HU", "beta4.psi.BS", "beta5.psi.lat", "beta6.psi.lon", 
              "beta8.psi.elev", "beta9.psi.dwd", "alpha1", "alpha2")
  
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
  
  
  # detection probability = 0.163
  plogis(-1.636) # back transform the mean of alpha 0
  
  
  
## OSS

  mcmc.output.o <- list()
  mcmc.output.o$chain1 <- as.matrix(O[[1]])
  mcmc.output.o$chain2 <- as.matrix(O[[2]])
  mcmc.output.o$chain3 <- as.matrix(O[[3]])

  attach.nimble(mcmc.output.o)
  
# Table with estimates and CI
  
  params <- c("beta0.psi", "beta0.theta", "alpha0", "beta1.psi.BU", "beta2.psi.HB", 
              "beta3.psi.HU", "beta4.psi.BS", "beta5.psi.lat", "beta6.psi.lon", 
              "beta8.psi.elev", "beta9.psi.dwd", "alpha1", "alpha2")
  
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
  
  write.csv(summary_table_o, "o_parameter_summary.csv", row.names = FALSE)
  
  
  # detection probability = 0.256
  plogis(-1.068) # back transform the mean of alpha 0

  
  
#### Calculate predicted occupancy probabilities for each treatment ---------------- 
  
  
## ENES
  
  if("beta0.psi" %in% summary_table_e$Parameter) {
    # get the intercept (control)
    intercept <- summary_table_e$Mean[summary_table_e$Parameter == "beta0.psi"]
    
    # reference treatment occupancy - convert logit to probability scale
    ref_occupancy <- plogis(intercept)
    
    # add the treatment effect to the control estimate, then convert to prob scale
    occupancy_results_e <- data.frame(
      Parameter = treatment_effects$Parameter,
      Occupancy = plogis(intercept + treatment_effects$Mean),
      Percent_Change = ((plogis(intercept + treatment_effects$Mean) - ref_occupancy) / ref_occupancy) * 100
    )
    
    # Add reference row at the top
    occupancy_results_e <- rbind(
      data.frame(Parameter = "Reference", Occupancy = ref_occupancy, Percent_Change = 0),
      occupancy_results_e
    )
    
    print(occupancy_results_e)
  }

  
  
## OSS
  
  if("beta0.psi" %in% summary_table_o$Parameter) {
    # get the intercept (control)
    intercept <- summary_table_o$Mean[summary_table_o$Parameter == "beta0.psi"]
    
    # reference treatment occupancy - convert logit to probability scale
    ref_occupancy <- plogis(intercept)
    
    # add the treatment effect to the control estimate, then convert to prob scale
    occupancy_results_o <- data.frame(
      Parameter = treatment_effects$Parameter,
      Occupancy = plogis(intercept + treatment_effects$Mean),
      Percent_Change = ((plogis(intercept + treatment_effects$Mean) - ref_occupancy) / ref_occupancy) * 100
    )
    
    # Add reference row at the top
    occupancy_results_o <- rbind(
      data.frame(Parameter = "Reference", Occupancy = ref_occupancy, Percent_Change = 0),
      occupancy_results_o
    )
    
    print(occupancy_results_o)
  }
  
  
  
  
#### CI, percent neg, using old methods from Claude - not using ----------------------------------
  
  # OSS 
  
  # Calculate 95% credible intervals
  # significant effect if CI does not include zero
  quantile(O2[,"beta1.psi.BU"], c(0.025, 0.975))    # Burned-Unharvested
  quantile(O2[,"beta2.psi.HB"], c(0.025, 0.975))    # Harvested-Burned  
  quantile(O2[,"beta3.psi.HU"], c(0.025, 0.975))    # Harvested-Unburned      ## excludes zero
  quantile(O2[,"beta4.psi.BS"], c(0.025, 0.975))    # Burned-Salvage          ## excludes zero
  # these are the same as in summary(O)
  
  
  # Check proportion of samples < 0 
  # significant neg effect if 95% samples <0
  mean(O2[,"beta1.psi.BU"] < 0)  # BU           ## no effect
  mean(O2[,"beta2.psi.HB"] < 0)  # HB           ## moderate neg effect
  mean(O2[,"beta3.psi.HU"] < 0)  # HU           ## significant neg effect
  mean(O2[,"beta4.psi.BS"] < 0)  # BS           ## significant neg effect
  
  
  # Probability that burn-salvage is worse than harvest-only
  c1 <- mean(O2[,"beta4.psi.BS"] < O2[,"beta3.psi.HU"])
  
  
  # Probability that burn-salvage is worse than harvest-burn
  c2 <- mean(O2[,"beta4.psi.BS"] < O2[,"beta2.psi.HB"])
  
  # Probability that harvest-burn is higher than harvest-only
  c5 <- mean(O2[,"beta2.psi.HB"] > O2[,"beta3.psi.HU"])
  
  # Convert to percentages for reporting
  cat("\nFor reporting:\n")
  cat("Harvest-only lower than controls:", round(mean(O2[,"beta3.psi.HU"] < 0) * 100, 1), "%\n")
  cat("Burn-salvage lower than controls:", round(mean(O2[,"beta4.psi.BS"] < 0) * 100, 1), "%\n")
  cat("Burn-salvage worse than harvest-only:", round(c1 * 100, 1), "%\n")
  cat("Burn-salvage worse than harvest-burn:", round(c2 * 100, 1), "%\n")
  cat("Harvest-burn higher than harvest-only:", round(c5 * 100, 1), "%\n")
  
  
  
  # ENES 
  
  # Calculate 95% credible intervals
  # significant effect if CI does not include zero
  quantile(E2[,"beta1.psi.BU"], c(0.025, 0.975))    # Burned-Unharvested
  quantile(E2[,"beta2.psi.HB"], c(0.025, 0.975))    # Harvested-Burned        ## excludes zero
  quantile(E2[,"beta3.psi.HU"], c(0.025, 0.975))    # Harvested-Unburned      ## excludes zero
  quantile(E2[,"beta4.psi.BS"], c(0.025, 0.975))    # Burned-Salvage          ## excludes zero
  # these are the same as in summary(E)
  
  
  # Check proportion of samples < 0 
  # significant neg effect if 95% samples <0
  mean(E2[,"beta1.psi.BU"] < 0)  # BU           ## moderate negative effect
  mean(E2[,"beta2.psi.HB"] < 0)  # HB           ## significant negative effect
  mean(E2[,"beta3.psi.HU"] < 0)  # HU           ## significant negative effect
  mean(E2[,"beta4.psi.BS"] < 0)  # BS           ## significant negative effect
  
  
  # Probability that burn-salvage is worse than harvest-only
  c3 <- mean(E2[,"beta4.psi.BS"] < E2[,"beta3.psi.HU"])
  
  
  # Probability that burn-salvage is worse than harvest-burn
  c4 <- mean(E2[,"beta4.psi.BS"] < E2[,"beta2.psi.HB"])
  
  
  # Convert to percentages for reporting
  cat("\nFor reporting:\n")
  cat("Harvest-only lower than controls:", round(mean(E2[,"beta3.psi.HU"] < 0) * 100, 1), "%\n")
  cat("Burn-salvage lower than controls:", round(mean(E2[,"beta4.psi.BS"] < 0) * 100, 1), "%\n")
  cat("Burn-salvage worse than harvest-only:", round(c3 * 100, 1), "%\n")
  cat("Burn-salvage worse than harvest-burn:", round(c4 * 100, 1), "%\n")
  
  
  