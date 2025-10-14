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
  library(ggtext)
  library(patchwork)


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
        parameter == "beta1.psi.BU" ~ "Burn-only",
        parameter == "beta2.psi.HB" ~ "Harvest-Burn", 
        parameter == "beta3.psi.HU" ~ "Harvest-only",
        parameter == "beta4.psi.BS" ~ "Burn-Salvage",
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
  
  
  # Define colors for each parameter
  param_colors <- c(
    "Downed Wood" = "steelblue",
    "Elevation" = "steelblue",
    "Latitude" = "steelblue",
    "Burn-only" = "#69995D",
    "Longitude" = "steelblue",
    "Harvest-Burn" = "#FEC601",
    "Harvest-only" = "#62B6CB",
    "Burn-Salvage" = "#EA7317"
  )
  
  
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
    
    ggplot() + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    
    # horizontal baseline for each parameter
    geom_hline(yintercept = 1:length(levels(posterior_data_O$parameter_clean)), 
               color = "gray20", alpha = 0.7, linewidth = 0.5) +
    
    # density curves
    geom_ribbon(aes(x = x, ymin = y_pos, ymax = y_final, group = parameter_clean, fill = parameter_clean),
                alpha = 0.7, linewidth = 0.3) +
    scale_fill_manual(values = param_colors, guide = "none") +
    
    # parameter labels
    scale_y_continuous(
      breaks = 1:length(levels(posterior_data_O$parameter_clean)),
      labels = levels(posterior_data_O$parameter_clean)
    ) +
    # Narrow the x-axis to focus on the distributions
    coord_cartesian(xlim = c(-5, 3)) +
    labs(
      #x = expression(paste("Effect Size: ", italic("B. wrighti")))
      x = "Effect Size"
      ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(vjust = 0, size = 16),      # bigger y-axis labels
      axis.text.x = element_text(size = 14),                 # bigger x-axis labels
      axis.title.x = element_text(size = 16),                # bigger x-axis title
    )
  
  
  print(plot_densities_O)
  
  ggsave("figures/o-posterior-dist.png", 
         plot = plot_densities_O, 
         width = 6,      
         height = 10,  
         dpi = 300)
  
  
  
#### ENES -----------------------------------------------------------------------
  
  # Convert model output and prepare data
  posterior_data_E <- as.data.frame(E2) %>%
    select(beta1.psi.BU, beta2.psi.HB, beta3.psi.HU, beta4.psi.BS, 
           beta5.psi.lat, beta6.psi.lon, beta8.psi.elev, beta9.psi.dwd) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(
      parameter_clean = case_when(
        parameter == "beta1.psi.BU" ~ "Burn-only",
        parameter == "beta2.psi.HB" ~ "Harvest-Burn", 
        parameter == "beta3.psi.HU" ~ "Harvest-only",
        parameter == "beta4.psi.BS" ~ "Burn-Salvage",
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
  
  
  # Define colors for each parameter
  param_colors <- c(
    "Downed Wood" = "steelblue",
    "Elevation" = "steelblue",
    "Latitude" = "steelblue",
    "Burn-only" = "#69995D",
    "Longitude" = "steelblue",
    "Harvest-Burn" = "#FEC601",
    "Harvest-only" = "#62B6CB",
    "Burn-Salvage" = "#EA7317"
  )
  
  
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
    
    ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    
    # horizontal baseline for each parameter
    geom_hline(yintercept = 1:length(levels(posterior_data_E$parameter_clean)), 
               color = "gray20", alpha = 0.7, size = 0.5) +
    
    # density curves
    geom_ribbon(aes(x = x, ymin = y_pos, ymax = y_final, group = parameter_clean, fill = parameter_clean),
                alpha = 0.7, linewidth = 0.3) +
    scale_fill_manual(values = param_colors, guide = "none") +
    
    # parameter labels
    scale_y_continuous(
      breaks = 1:length(levels(posterior_data_E$parameter_clean)),
      labels = levels(posterior_data_E$parameter_clean)
    ) +
    coord_cartesian(xlim = c(-5, 3)) +
    labs(
      #x = expression(paste("Effect Size: ", italic("E. eschscholtzii")))
      x = "Effect Size"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(vjust = 0, size = 16),      # bigger y-axis labels
      axis.text.x = element_text(size = 14),                 # bigger x-axis labels
      axis.title.x = element_text(size = 16),                # bigger x-axis title
    )
  
  
  print(plot_densities_E)
  
  ggsave("figures/e-posterior-dist.png", 
         plot = plot_densities_E, 
         width = 6,      
         height = 10,  
         dpi = 300)
  
  
  
  
## COMBINE
  
  p <- plot_densities_O + plot_densities_E + 
    plot_annotation(
      tag_levels = 'A',
      title = "Posterior Distributions of Covariate Effects on Occupancy"
    ) &
    theme(
      plot.tag = element_text(size = 18, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  
  
  ggsave("figures/both-posterior-dist-facet.png", 
         plot = p, 
         width = 12,      
         height = 10,  
         dpi = 300)
  
  
