## =================================================
##
## Title: raw-data-summary-for-pub
## Author: Jasmine Williamson 
## Date Created: 09/08/2025
##
## Description: Summary stats on RAW DATA for paper results section
##
## =================================================
setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/multi-scale-occu-oss")


## Setup 
  library(dplyr)

# Load data
  enes <- read.csv("data/enes.prepost.multiscale.occu.csv") 
  oss <- read.csv("data/oss.prepost.multiscale.occu.csv")
  
  
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
  total_sites_surveyed_e_2020 <- nrow(enes_dets[enes_dets$year >= 2020, ])
  
  # detections before and after 2020 wildfire
  pre_2020_e <- sum(enes_dets$detected[enes_dets$year < 2020], na.rm = TRUE)
  post_2020_e <- sum(enes_dets$detected[enes_dets$year >= 2020], na.rm = TRUE)
  total_site_dets_e <- sum(enes_dets$detected, na.rm = TRUE)
  
  
  # detections across all surveys (V1-V3) for post-2020 only
  total_survey_dets_e_post2020 <- sum(c(enes_clean$V1[enes_clean$year >= 2020], 
                                        enes_clean$V2[enes_clean$year >= 2020], 
                                        enes_clean$V3[enes_clean$year >= 2020]), na.rm = TRUE)
  
  # Percentage of sites with detections
  pct_sites_detected_e <- round((total_site_dets_e / total_sites_surveyed_e) * 100, 1)
  
  # Percentage of sites with detections - post-2020
  pct_sites_detected_e_2020 <- round((post_2020_e / total_sites_surveyed_e_2020) * 100, 1)
  
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

  
  ## All yrs data ##
  total_sites_surveyed_e # number of total sites surveyed
  total_dets_e # total individual enes detections across all surveys
  pre_2020_e # site level detections before 2020 wildfire
  total_site_dets_e # total sites where enes was detected
  pct_sites_detected_e # percentage of sites with detections
  
  ## 2023-2024 data ##
  post_2020_e # site detections post-fire
  total_sites_surveyed_e_2020 # total surveyed plots post fire
  total_survey_dets_e_post2020 # total enes detections post-fire
  pct_sites_detected_e_2020 # percent of sites with enes detections
  
  
  
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
  total_sites_surveyed_o_2020 <- nrow(oss_dets[oss_dets$year >= 2020, ])
  
  
  # detections before and after 2020 wildfire
  pre_2020_o <- sum(oss_dets$detected[oss_dets$year < 2020], na.rm = TRUE)
  post_2020_o <- sum(oss_dets$detected[oss_dets$year >= 2020], na.rm = TRUE)
  total_site_dets_o <- sum(oss_dets$detected, na.rm = TRUE)
  
  # detections across all surveys (V1-V3)
  total_survey_dets_o <- sum(c(oss_clean$V1, oss_clean$V2, oss_clean$V3), na.rm = TRUE)
  
  # detections across all surveys (V1-V3) for post-2020 only - OSS
  total_survey_dets_o_post2020 <- sum(c(oss_clean$V1[oss_clean$year >= 2020], 
                                        oss_clean$V2[oss_clean$year >= 2020], 
                                        oss_clean$V3[oss_clean$year >= 2020]), na.rm = TRUE)
  
  # Percentage of sites with detections - all
  pct_sites_detected_o <- round((total_site_dets_o / total_sites_surveyed_o) * 100, 1)
  
  # Percentage of sites with detections - post-2020
  pct_sites_detected_o_2020 <- round((post_2020_o / total_sites_surveyed_o_2020) * 100, 1)
  
  
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
  
  
  
  ## All yrs data ##
  total_sites_surveyed_o # number of total sites surveyed
  total_dets_o # total individual enes detections across all surveys
  pre_2020_o # site level detections before 2020 wildfire
  total_site_dets_o # total sites where enes was detected
  pct_sites_detected_o # percentage of sites with detections
  
  ## 2023-2024 data ##
  post_2020_o # site detections post-fire
  total_sites_surveyed_o_2020 # total surveyed plots post fire
  total_survey_dets_o_post2020 # total enes detections post-fire
  pct_sites_detected_o_2020 # percent of sites with enes detections
    
    
  
#### Combines Species ---------------------------------------------------------
    
    # Combine species at site level
    combined_site_dets <- enes_dets %>%
      full_join(oss_dets, by = c("site_id", "year", "trt")) %>%
      mutate(
        enes_detected = ifelse(is.na(detected.x), 0, detected.x),
        oss_detected = ifelse(is.na(detected.y), 0, detected.y),
        any_species_detected = pmax(enes_detected, oss_detected, na.rm = TRUE)
      )
  
    # Subset df for only 2023-2024
    combined_site_dets <- combined_site_dets[combined_site_dets$year %in% c(2023, 2024), ]
    
    # Calculate percentage of sites with either species
    total_sites_combined <- nrow(combined_site_dets)
    sites_with_either_species <- sum(combined_site_dets$any_species_detected, na.rm = TRUE)
    pct_sites_either_species <- round((sites_with_either_species / total_sites_combined) * 100, 1)
    
    pct_sites_either_species  
  
    
    
    
#### Number of stands ---------------------------------------------------------
    
    # unique stands from each period
    stands_2013_2019 <- unique(enes$stand[enes$year >= 2013 & enes$year <= 2019])
    length(stands_2013_2019)
    
    stands_2023_2024 <- unique(enes$stand[enes$year >= 2023 & enes$year <= 2024])
    length(stands_2023_2024)
    
    # stands that were surveyed in BOTH periods (repeated stands)
    repeated_stands <- intersect(stands_2013_2019, stands_2023_2024)
    length(repeated_stands)
    
    # new stands in 2023-2024 only
    new_stands <- setdiff(stands_2023_2024, stands_2013_2019)
    length(new_stands)
 
    
    
    # Filter to unique site-year combos
    df_new <- enes %>%
      select(site_id, stand, year, trt) %>%
      distinct()
    
    # Split into old study and new study
    df_old <- filter(df_new, year < 2023)
    df_new <- filter(df_new, year >= 2023)
    
    # Check whether stand was resurveyed or whether it is new
    df_new$resurvey <- ifelse(df_new$stand %in% df_old$stand, "resurvey", "new")
    
    # Check counts of resurveys/new by treatment
    table(df_new$trt, df_new$resurvey)
    
    
    
    