## =================================================
##
## Title: occu-summary-stats-for-pub
## Author: Jasmine Williamson 
## Date Created: 09/08/2025
##
## Description: Summary stats on occupancy for paper results section
##
## =================================================
setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/multi-scale-occu-oss")


## Setup 
  library(dplyr)

# Load data
  e_covs <- read.csv("data/enes.prepost.multiscale.occu.csv") 
  o_covs <- read.csv("data/oss.prepost.multiscale.occu.csv")
  
  
#### ENES ----------------------------------------------------------------------
  
  # remove NA
  enes_clean <- enes %>% 
    filter(!is.na(V1) | !is.na(V2) | !is.na(V3))
  
  # total individual detections across all surveys
  total_dets_e <- sum(c(enes_clean$V1, enes_clean$V2, enes_clean$V3), na.rm = TRUE)
  
  # detection status by site (detected = 1 if ANY plot within site had detection)
  enes_dets <- enes_clean %>%
    mutate(
      plot_detected = pmax(V1, V2, V3, na.rm = TRUE),  # 1 if detected in any survey, 0 if not
      plot_detected = ifelse(is.na(plot_detected), 0, plot_detected)
    ) %>%
    group_by(site_id, year, trt) %>%
    summarise(
      detected = max(plot_detected, na.rm = TRUE),  # 1 if any plot in site had detection
      .groups = 'drop'
    )
  
  total_sites_surveyed_e <- nrow(enes_dets)
  
  
  # detections before and after 2020 wildfire
  pre_2020_e <- sum(enes_dets$detected[enes_dets$year < 2020], na.rm = TRUE)
  post_2020_e <- sum(enes_dets$detected[enes_dets$year >= 2020], na.rm = TRUE)
  total_site_dets_e <- sum(enes_dets$detected, na.rm = TRUE)
  
  # detections across all surveys (V1-V3)
  total_survey_dets_e <- sum(c(enes_clean$V1, enes_clean$V2, enes_clean$V3), na.rm = TRUE)
  
  
  # Percentage of sites with detections
  pct_sites_detected_e <- round((total_site_dets_e / total_sites_surveyed_e) * 100, 1)
  
  # Overall treatment summary (at site level)
  treatment_summary <- enes_dets %>%
    group_by(trt) %>%
    summarise(
      mean_occupancy = round(mean(detected, na.rm = TRUE), 3),
      sites_detected = sum(detected, na.rm = TRUE),
      total_sites = n(),
      detection_rate = round(sum(detected, na.rm = TRUE) / n(), 3),
      .groups = 'drop'
    )
  
  
  # Detection frequency across surveys
  detection_pattern <- enes_clean %>%
    rowwise() %>%
    mutate(num_detections = sum(c(V1, V2, V3), na.rm = TRUE)) %>%
    ungroup() %>%
    count(num_detections, name = "plots") %>%
    mutate(percentage = round(plots/sum(plots)*100, 1))

  
  
  total_dets_e # total individual detections across all surveys
  total_sites_surveyed_e
  pre_2020_e # site detections before and after 2020 wildfire
  post_2020_e
  total_site_dets_e # total sites where enes was detected
  total_survey_dets_e # total detections across all surveys (V1-V3)
  pct_sites_detected_e # percentage of sites with detections

  
  
#### OSS ----------------------------------------------------------------------
  
  # remove NA
  oss_clean <- oss %>% 
    filter(!is.na(V1) | !is.na(V2) | !is.na(V3))
  
  # total individual detections across all surveys
  total_dets_o <- sum(c(oss_clean$V1, oss_clean$V2, oss_clean$V3), na.rm = TRUE)
  
  # detection status by site (detected = 1 if ANY plot within site had detection)
  oss_dets <- oss_clean %>%
    mutate(
      plot_detected = pmax(V1, V2, V3, na.rm = TRUE),  # 1 if detected in any survey, 0 if not
      plot_detected = ifelse(is.na(plot_detected), 0, plot_detected)
    ) %>%
    group_by(site_id, year, trt) %>%
    summarise(
      detected = max(plot_detected, na.rm = TRUE),  # 1 if any plot in site had detection
      .groups = 'drop'
    )
  
  total_sites_surveyed_o <- nrow(oss_dets)
  
  
  # detections before and after 2020 wildfire
  pre_2020_o <- sum(oss_dets$detected[oss_dets$year < 2020], na.rm = TRUE)
  post_2020_o <- sum(oss_dets$detected[oss_dets$year >= 2020], na.rm = TRUE)
  total_site_dets_o <- sum(oss_dets$detected, na.rm = TRUE)
  
  # detections across all surveys (V1-V3)
  total_survey_dets_o <- sum(c(oss_clean$V1, oss_clean$V2, oss_clean$V3), na.rm = TRUE)
  
  
  # Percentage of sites with detections
  pct_sites_detected_o <- round((total_site_dets_o / total_sites_surveyed_o) * 100, 1)
  
  # Overall treatment summary (at site level)
  treatment_summary <- oss_dets %>%
    group_by(trt) %>%
    summarise(
      mean_occupancy = round(mean(detected, na.rm = TRUE), 3),
      sites_detected = sum(detected, na.rm = TRUE),
      total_sites = n(),
      detection_rate = round(sum(detected, na.rm = TRUE) / n(), 3),
      .groups = 'drop'
    )
  
  
  # Detection frequency across surveys
  detection_pattern <- oss_clean %>%
    rowwise() %>%
    mutate(num_detections = sum(c(V1, V2, V3), na.rm = TRUE)) %>%
    ungroup() %>%
    count(num_detections, name = "plots") %>%
    mutate(percentage = round(plots/sum(plots)*100, 1))
  
  
  
  
    total_dets_o # total individual detections across all surveys
    total_sites_surveyed_o
    pre_2020_o # site detections before and after 2020 wildfire
    post_2020_o
    total_site_dets_o # total sites where enes was detected
    total_survey_dets_o # total detections across all surveys (V1-V3)
    pct_sites_detected_o # percentage of sites with detections
  
  
  
  
  