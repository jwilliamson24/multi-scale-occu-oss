## =================================================
##
## Title: posterior-distribution-plots
## Author: Jasmine Williamson 
## Date Created: 09/09/2025
##
## Description: Creating posterior distribution plots for covariates from
##              multi-scale hierarchical occupancy model output
##
## =================================================
setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/multi-scale-occu-oss")


## Setup 
  library(ggplot2)
  library(dplyr)
  library(tidyr)

# Load data
  load("data/multiscale_output_and_data_082525_enes_full.RData")
  E = a
  load("data/multiscale_output_and_data_082525_oss_full.RData")
  O = a

# Combine
  summary(E)
  summary(O)
  E2 = runjags::combine.mcmc(E)
  O2 = runjags::combine.mcmc(O)  


#### OSS -----------------------------------------------------------------------
  

  # Convert model output and prepare data
  posterior_data_O <- as.data.frame(O2) %>%
    select(beta1.psi.BU, beta2.psi.HB, beta3.psi.HU, beta4.psi.BS, 
           beta5.psi.lat, beta6.psi.lon, beta8.psi.elev, beta9.psi.dwd) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(
      parameter_clean = case_when(
        parameter == "beta1.psi.BU" ~ "Burned-Unharvested",
        parameter == "beta2.psi.HB" ~ "Harvested-Burned", 
        parameter == "beta3.psi.HU" ~ "Harvested-Unburned",
        parameter == "beta4.psi.BS" ~ "Burned-Salvage",
        parameter == "beta5.psi.lat" ~ "Latitude",
        parameter == "beta6.psi.lon" ~ "Longitude",
        parameter == "beta8.psi.elev" ~ "Elevation", 
        parameter == "beta9.psi.dwd" ~ "Downed Wood"
      )
    ) %>%
    # Order by mean effect size
    group_by(parameter_clean) %>%
    mutate(mean_val = mean(value)) %>%
    ungroup() %>%
    mutate(parameter_clean = reorder(parameter_clean, mean_val))
  
  
  # Plot
  plot_densities_O <- posterior_data_O %>%
    # Calculate density for each parameter
    group_by(parameter_clean) %>%
    do(data.frame(x = density(.$value, n = 512)$x, 
                  y = density(.$value, n = 512)$y)) %>%
    # Scale densities
    mutate(
      y_pos = as.numeric(parameter_clean),
      y_final = y_pos + (y / max(y)) * 0.8  # Scale and offset in one step
    ) %>%
    
    # dashed line at 0
    ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    
    # Add horizontal baseline for each parameter
    geom_hline(yintercept = 1:length(levels(posterior_data_O$parameter_clean)), 
               color = "gray20", alpha = 0.7, size = 0.5) +
    
    # Draw density curves
    geom_ribbon(aes(x = x, ymin = y_pos, ymax = y_final, group = parameter_clean),
                fill = "steelblue", alpha = 0.7, size = 0.3) +
    
    # Add parameter labels
    scale_y_continuous(
      breaks = 1:length(levels(posterior_data_O$parameter_clean)),
      labels = levels(posterior_data_O$parameter_clean)
    ) +
    
    labs(
      x = "Coefficient Value",
      title = "OSS: Posterior Distributions",
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(vjust = 0, size = 12)
    )
  
  print(plot_densities_O)
  
  ggsave("figures/o-posterior-dist.png", plot = plot_densities_O, dpi = 300)
  
  
  
#### ENES -----------------------------------------------------------------------
  
  
  # Convert model output and prepare data
  posterior_data_E <- as.data.frame(E2) %>%
    select(beta1.psi.BU, beta2.psi.HB, beta3.psi.HU, beta4.psi.BS, 
           beta5.psi.lat, beta6.psi.lon, beta8.psi.elev, beta9.psi.dwd) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(
      parameter_clean = case_when(
        parameter == "beta1.psi.BU" ~ "Burned-Unharvested",
        parameter == "beta2.psi.HB" ~ "Harvested-Burned", 
        parameter == "beta3.psi.HU" ~ "Harvested-Unburned",
        parameter == "beta4.psi.BS" ~ "Burned-Salvage",
        parameter == "beta5.psi.lat" ~ "Latitude",
        parameter == "beta6.psi.lon" ~ "Longitude",
        parameter == "beta8.psi.elev" ~ "Elevation", 
        parameter == "beta9.psi.dwd" ~ "Downed Wood"
      )
    ) %>%
    # Order by mean effect size
    group_by(parameter_clean) %>%
    mutate(mean_val = mean(value)) %>%
    ungroup() %>%
    mutate(parameter_clean = reorder(parameter_clean, mean_val))
  
  
  # Plot
  plot_densities_E <- posterior_data_E %>%
    # Calculate density for each parameter
    group_by(parameter_clean) %>%
    do(data.frame(x = density(.$value, n = 512)$x, 
                  y = density(.$value, n = 512)$y)) %>%
    # Scale densities
    mutate(
      y_pos = as.numeric(parameter_clean),
      y_final = y_pos + (y / max(y)) * 0.8  # Scale and offset in one step
    ) %>%
    
    # dashed line at 0
    ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    
    # Add horizontal baseline for each parameter
    geom_hline(yintercept = 1:length(levels(posterior_data_O$parameter_clean)), 
               color = "gray20", alpha = 0.7, size = 0.5) +
    
    # Draw density curves
    geom_ribbon(aes(x = x, ymin = y_pos, ymax = y_final, group = parameter_clean),
                fill = "steelblue", alpha = 0.7, size = 0.3) +
    
    # Add parameter labels
    scale_y_continuous(
      breaks = 1:length(levels(posterior_data_O$parameter_clean)),
      labels = levels(posterior_data_O$parameter_clean)
    ) +
    
    labs(
      x = "Coefficient Value",
      title = "ENES: Posterior Distributions",
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(vjust = 0, size = 12)
    )
  
  print(plot_densities_E)
  
  ggsave("figures/e-posterior-dist.png", plot = plot_densities_E, dpi = 300)
  
  
  

