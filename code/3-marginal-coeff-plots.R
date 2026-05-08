## =================================================
##
## Title: marginal-coeff-plots
## Author: Jasmine Williamson 
## Date Created: 07/22/2025
##
## Description: Visualize and interpret outputs for 
## multi-scale occupancy model, both spp, all-years dataset
## added detection prob for each species at mean temp
##
## =================================================
setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/multi-scale-occu-oss")

##### Setup --------------------------------------------
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggh4x)


# Load data
  e_covs <- read.csv("data/enes.prepost.multiscale.occu.csv") 
  o_covs <- read.csv("data/oss.prepost.multiscale.occu.csv")
  load("data/msc-data-workspace-V3.RData")
  load("data/multiscale_output_121225_enes_small.RData")
  E = a2
  load("data/multiscale_output_121225_oss_small.RData")
  O = a2
  
# Estimates
  summary(E)
  summary(O)

# Combine
  E2 = runjags::combine.mcmc(E)
  O2 = runjags::combine.mcmc(O)   


##### Posterior Coefficient Estimates Plot -------------------------------
  
  summary_table <- summary(O) # # # # # choose species # # # # #
  
  params_to_plot <- c("beta0.psi", "beta0.theta", "alpha0", 
                      "beta1.psi.BU", "beta2.psi.HB", "beta3.psi.HU", "beta4.psi.BS", 
                      "beta5.psi.lat", "beta6.psi.lon", "beta8.psi.elev", "beta1.theta.DW",
                      "alpha1", "alpha2")
  
  # summarize
  coef_df <- data.frame(
    Parameter = rownames(summary_table$statistics),
    Mean = summary_table$statistics[, "Mean"],
    LCI = summary_table$quantiles[, "2.5%"],
    UCI = summary_table$quantiles[, "97.5%"]
  )
  
  # subset parameters
  coef_df <- subset(coef_df, Parameter %in% params_to_plot)
  
  # add submodel column to group by
  coef_df$Submodel <- case_when(
    coef_df$Parameter %in% c("beta0.psi", "beta1.psi.BU", "beta2.psi.HB", "beta3.psi.HU", 
                             "beta4.psi.BS", "beta5.psi.lat", "beta6.psi.lon", 
                             "beta8.psi.elev") ~ "Site Occupancy (ψ)",
    
    coef_df$Parameter %in% c("beta0.theta", "beta1.theta.DW") ~ "Subplot Use (θ)",
    
    coef_df$Parameter %in% c("alpha0", "alpha1", "alpha2") ~ "Detection (p)",
    
    TRUE ~ "Other"
  )
  
  # order submodel groupings in plot
  coef_df$Submodel <- factor(coef_df$Submodel, levels = c(
    "Site Occupancy (ψ)",
    "Subplot Use (θ)", 
    "Detection (p)"
  ))
  
  # rename - optional
  coef_df$Parameter <- dplyr::recode(coef_df$Parameter,
                              "beta1.psi.BU" = "Burned (BU)",
                              "beta2.psi.HB" = "Harvest+Burn (HB)",
                              "beta3.psi.HU" = "Harvested (HU)",
                              "beta4.psi.BS" = "Burn+Salvage (BS)",
                              "beta0.psi" = "Site Occu Intercept",
                              "beta0.theta" = "Plot Use Intercept",
                              "alpha0" = "Detection Intercept",
                              "beta5.psi.lat" = "Latitude",
                              "beta6.psi.lon" = "Longitude",
                              "beta8.psi.elev" = "Elevation",
                              "beta1.theta.DW" = "Downed Wood",
                              "alpha1" = "Linear Temp",
                              "alpha2" = "Quadratic Temp")
  
  
  p <- ggplot(coef_df, aes(x = Mean, y = reorder(Parameter, Mean))) +
    geom_point() +
    geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_wrap(~ Submodel, scales = "free_y", ncol=1) +
    coord_cartesian(xlim = c(min(coef_df$LCI) - 0.2, max(coef_df$UCI) + 0.2)) +
    labs(title = "Posterior Coefficient Estimates - ENES",
         x = "Estimate (logit scale)",
         y = "Parameter") +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13),
      legend.position = "none",
      strip.text = element_text(size = 13, face = "bold"),
      plot.title = element_text(face = "bold")
    )
  
  p
  
  ggsave("figures/o-coeff-plot.png", plot = p, dpi = 300)
  
  
##### Marginal Psi Predictions Dataframe - ENES -------------------------------------------  
  
# Isolate the effect of a single covariate by varying it across a range,
# holding everything else constant
  
  b <- E2    # # # # # # choose species here # # # # # #
  dat <- e_covs    # # # # # # choose species here # # # # # #
  n.samples = nrow(b) # number of posterior samples
  
# Latitude on psi
  
  # range of lat
  r<- range(lat.2D)
  dim(lat.2D)
  
  # create sequence along range
  lat_data <- seq(r[1], r[2], length.out=50) 
  
  # set all other covs at mean
  HB = 0
  HU = 0
  BS = 0
  BU = 0
  long = 0
  elev = 0
  
  # create matrices to stick estimates in
  logit_lat_psi = matrix(NA, n.samples, length(lat_data))
  
  # psi predictions for enes
  # Sample from posterior for the sequence of values for cov - raw linear predictor for each MCMC sample 
  for (i in 1:n.samples){
    for (j in 1:length(lat_data)){
      # create the linear predictors for dominant species
      logit_lat_psi[i,j] = 
        b[,'beta0.psi'][[i]] + 
        b[,'beta1.psi.BU'][[i]] * BU + 
        b[,'beta2.psi.HB'][[i]] * HB + 
        b[,'beta3.psi.HU'][[i]] * HU + 
        b[,'beta4.psi.BS'][[i]] * BS + 
        b[,'beta5.psi.lat'][[i]] * lat_data[j] +  
        b[,'beta6.psi.lon'][[i]] * long +  
        b[,'beta8.psi.elev'][[i]] * elev 
    }}
  
  # create array 
  lat_psi = matrix(NA, n.samples, length(lat_data))
  
  # transform psi off logit-scale back to probability scale - for each sample
  lat_psi <- plogis(logit_lat_psi)
  
  # calculate means and credible intervals - mean predicted occupancy probability across all posterior samples
  lat_psi_means = colMeans(lat_psi) # - taking the mean of all samples (transformed occu prob)
  lat_psi_CIs <- apply(lat_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)
  
  # stuff into df
  lat_psi_preds <- data.frame(predicted = lat_psi_means, 
                             cov_zsc = lat_data,
                             LCI = lat_psi_CIs[1,],
                             UCI = lat_psi_CIs[2,])
  
  # add back-transformed values to df (just in case i want them later)
  lat_mean  <- mean(dat$lat, na.rm = TRUE)
  lat_sd    <- sd(dat$lat, na.rm = TRUE)
  lat_psi_preds$cov_value  <- lat_psi_preds$cov_zsc  * lat_sd  + lat_mean
  
  # add species column
  lat_psi_preds_e = lat_psi_preds
  lat_psi_preds_e$species = "ENES"
  
  # e.lat <- ggplot(lat_psi_preds_e, aes(x = cov_value, y = predicted)) +
  #       geom_line() +
  #       geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.2) +
  #       ylab(bquote("Predicted "*psi~"")) +
  #       xlab("Latitude") +
  #       #labs(title = "Marginal Effect of Covariates on ENES Occupancy") +
  #       theme_classic() +
  #       theme(legend.position = "none",
  #             strip.text = element_text(size = 12, face = "bold"),
  #             plot.title = element_text(hjust = 0.5, face = "bold"),
  #             axis.text = element_text(size = 18),      # Axis numbers
  #             axis.title = element_text(size = 20)) +    # Axis labels  
  #       scale_x_continuous(n.breaks = 3)  # Reduce number of x-axis labels
  # 
  # #ggsave("figures/e-lat-effect.png", plot = e.lat, dpi = 300)
  
  
# Longitude on psi
  
  # range
  r<- range(lon.2D)
  # create sequence along range
  lon_data <- seq(r[1], r[2], length.out=50) 
  
  # set all other covs at mean
  HB = 0
  HU = 0
  BS = 0
  BU = 0
  lat = 0
  elev = 0
  
  # create matrices to stick estimates in
  logit_lon_psi = matrix(NA, n.samples, length(lon_data))
  
  # psi predictions for enes
  # Sample from posterior for the sequence of values for cov
  for (i in 1:n.samples){
    for (j in 1:length(lon_data)){
      logit_lon_psi[i,j] = 
        b[,'beta0.psi'][[i]] + 
        b[,'beta1.psi.BU'][[i]] * BU + 
        b[,'beta2.psi.HB'][[i]] * HB + 
        b[,'beta3.psi.HU'][[i]] * HU + 
        b[,'beta4.psi.BS'][[i]] * BS + 
        b[,'beta5.psi.lat'][[i]] * lat +  
        b[,'beta6.psi.lon'][[i]] * lon_data[j] +  
        b[,'beta8.psi.elev'][[i]] * elev 
    }}
  
  # create array 
  lon_psi = matrix(NA, n.samples, length(lon_data))
  
  # transform psi off logit-scale back to probability scale
  lon_psi <- plogis(logit_lon_psi)
  
  # calculate means and credible intervals 
  lon_psi_means = colMeans(lon_psi) 
  lon_psi_CIs <- apply(lon_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)
  
  # stuff into df
  lon_psi_preds <- data.frame(predicted = lon_psi_means, 
                              cov_zsc = lon_data,
                              LCI = lon_psi_CIs[1,],
                              UCI = lon_psi_CIs[2,])
  
  # add back-transformed values to df (just in case i want them later)
  long_mean  <- mean(dat$long, na.rm = TRUE)
  long_sd    <- sd(dat$long, na.rm = TRUE)
  lon_psi_preds$cov_value  <- lon_psi_preds$cov_zsc  * long_sd  + long_mean
  
  # add species column
  lon_psi_preds_e = lon_psi_preds
  lon_psi_preds_e$species = "ENES"
  
  
# Elevation on psi
  
  # range
  r<- range(elev.2D)
  # create sequence along range
  elev_data <- seq(r[1], r[2], length.out=50) 
  
  # set all other covs at mean
  HB = 0
  HU = 0
  BS = 0
  BU = 0
  lat = 0
  long = 0
  
  # create matrices to stick estimates in
  logit_elev_psi = matrix(NA, n.samples, length(elev_data))
  
  # psi predictions for enes
  # Sample from posterior for the sequence of values for cov 
  for (i in 1:n.samples){
    for (j in 1:length(elev_data)){
      logit_elev_psi[i,j] = 
        b[,'beta0.psi'][[i]] + 
        b[,'beta1.psi.BU'][[i]] * BU + 
        b[,'beta2.psi.HB'][[i]] * HB + 
        b[,'beta3.psi.HU'][[i]] * HU + 
        b[,'beta4.psi.BS'][[i]] * BS + 
        b[,'beta5.psi.lat'][[i]] * lat +  
        b[,'beta6.psi.lon'][[i]] * long +  
        b[,'beta8.psi.elev'][[i]] * elev_data[j] 
    }}
  
  # create array 
  elev_psi = matrix(NA, n.samples, length(elev_data))
  
  # transform psi off logit-scale back to probability scale
  elev_psi <- plogis(logit_elev_psi)
  
  # calculate means and credible intervals
  elev_psi_means = colMeans(elev_psi) 
  elev_psi_CIs <- apply(elev_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)
  
  # stuff into df
  elev_psi_preds <- data.frame(predicted = elev_psi_means, 
                               cov_zsc = elev_data,
                              LCI = elev_psi_CIs[1,],
                              UCI = elev_psi_CIs[2,])
  
  # add back-transformed values to df (just in case i want them later)
  elev_mean  <- mean(dat$elev, na.rm = TRUE)
  elev_sd    <- sd(dat$elev, na.rm = TRUE)
  elev_psi_preds$cov_value  <- round(elev_psi_preds$cov_zsc  * elev_sd  + elev_mean)
  
  # add species column
  elev_psi_preds_e = elev_psi_preds
  elev_psi_preds_e$species = "ENES"
  

# # DWD on psi
#   
#   # range
#   r<- range(downedwood.new)
#           scale_x_continuous(n.breaks = 4)  # Reduce number of x-axis labels
# # create sequence along range
#   dwd_data <- seq(r[1], r[2], length.out=50) 
#   
#   # set all other covs at mean
#   HB = 0
#   HU = 0
#   BS = 0
#   BU = 0
#   lat = 0
#   long = 0
#   elev = 0
#   
#   # create matrices to stick estimates in
#   logit_dwd_psi = matrix(NA, n.samples, length(dwd_data))
#   
#   # psi predictions for enes
#   # Sample from posterior for the sequence of values for cov 
#   for (i in 1:n.samples){
#     for (j in 1:length(dwd_data)){
#       logit_dwd_psi[i,j] = 
#         b[,'beta0.psi'][[i]] + 
#         b[,'beta1.psi.BU'][[i]] * BU + 
#         b[,'beta2.psi.HB'][[i]] * HB + 
#         b[,'beta3.psi.HU'][[i]] * HU + 
#         b[,'beta4.psi.BS'][[i]] * BS + 
#         b[,'beta5.psi.lat'][[i]] * lat +  
#         b[,'beta6.psi.lon'][[i]] * long +  
#         b[,'beta8.psi.elev'][[i]] * elev +
#         b[,'beta9.psi.dwd'][[i]] * dwd_data[j] 
#     }}
#   
#   # create array 
#   dwd_psi = matrix(NA, n.samples, length(dwd_data))
#   
#   # transform psi off logit-scale back to probability scale
#   dwd_psi <- plogis(logit_dwd_psi)
#   
#   # calculate means and credible intervals
#   dwd_psi_means = colMeans(dwd_psi) 
#   dwd_psi_CIs <- apply(dwd_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)
#   
#   # stuff into df
#   dwd_psi_preds <- data.frame(predicted = dwd_psi_means, 
#                               cov_zsc = dwd_data,
#                               LCI = dwd_psi_CIs[1,],
#                               UCI = dwd_psi_CIs[2,])
#   
#   # add back-transformed values to df (just in case i want them later)
#   dwd_mean  <- mean(dat$DW, na.rm = TRUE)
#   dwd_sd    <- sd(dat$DW, na.rm = TRUE)
#   dwd_psi_preds$cov_value  <- dwd_psi_preds$cov_zsc  * dwd_sd  + dwd_mean
#   
#   # add species column
#   dwd_psi_preds_e = dwd_psi_preds
#   dwd_psi_preds_e$species = "ENES"
#   
#   e.dwd <- ggplot(dwd_psi_preds_e, aes(x = cov_value, y = predicted)) +
#     geom_line(size = 1) +
#     geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.2) +
#     ylab(bquote("Predicted "*psi~"")) +
#     xlab("Downed Wood Count") +
#     #labs(title = "Marginal Effect of DWD on Occupancy") +
#     theme_classic() +
#     theme(legend.position = "none",
#           strip.text = element_text(size = 12, face = "bold"),
#           plot.title = element_text(hjust = 0.5, face = "bold"),
#           axis.text = element_text(size = 18),      # Axis numbers
#           axis.title = element_text(size = 20),
#           axis.title.y = element_blank(),    
#           axis.text.y  = element_blank()) +    
#     scale_x_continuous(n.breaks = 4)  # Reduce number of x-axis labels
  
    
# Combine
  # p.enes <- (e.lat | e.long | e.elev)
  
  # ggsave("figures/e-covs-effect.png",          
  #        width = 14,   # Wider
  #        height = 3,   # Taller to reduce the "skinny" look
  #        units = "in",
  #        plot = p.enes, dpi = 300)
  
  
##### Theta Prediction Dataframe - ENES --------------------------------------  
  
  b <- E2
  dat <- e_covs

# Downed Wood effect on Plot Use

  r<- range(downedwood.3D)
  # create sequence along range
  DW_data <- seq(r[1], r[2], length.out=50)

  # create matrices to stick estimates in
  logit_DW_theta = matrix(NA, n.samples, length(DW_data))

  # theta predictions for enes
  # Sample from posterior for the sequence of values for cov
  for (i in 1:n.samples){
    for (j in 1:length(DW_data)){
      logit_DW_theta[i,j] =
        b[,'beta0.theta'][[i]] +
        b[,'beta1.theta.DW'][[i]] * DW_data[j]
    }}

  # create array
  DW_theta = matrix(NA, n.samples, length(DW_data))

  # transform psi off logit-scale back to probability scale
  DW_theta <- plogis(logit_DW_theta)

  # calculate means and credible intervals
  DW_means = colMeans(DW_theta)
  DW_CIs <- apply(DW_theta,2,quantile, c(0.025,0.975), na.rm=TRUE)

  # stuff into df
  DW_preds <- data.frame(predicted = DW_means,
                             cov_zsc = DW_data,
                               LCI = DW_CIs[1,],
                               UCI = DW_CIs[2,])

  # add back-transformed values to df (just in case i want them later)
  dw_mean  <- mean(dat$DW, na.rm = TRUE)
  dw_sd    <- sd(dat$DW, na.rm = TRUE)
  DW_preds$cov_value  <- round(DW_preds$cov_zsc  * dw_sd  + dw_mean, 2)

  # add species
  DW_preds_e <- DW_preds
  DW_preds_e$species <- "ENES"


##### Detection Prediction Dataframe - ENES -------------------------------------
  
  b <- E2
  dat <- e_covs
  n.samples = nrow(b) # number of posterior samples
  
  r<- range(temp.3D)
  # create sequence along range
  temp_data <- seq(r[1], r[2], length.out=50)

  # create matrices to stick estimates in
  logit_temp = matrix(NA, n.samples, length(temp_data))

  # alpha predictions for enes
  # Sample from posterior for the sequence of values for cov
  for (i in 1:n.samples){
    for (j in 1:length(temp_data)){
      logit_temp[i,j] =
        # alpha0_mean[i] +
        b[,'alpha0'][[i]] +
        b[,'alpha1'][[i]] * temp_data[j] +
        b[,'alpha2'][[i]] * temp_data[j]^2
    }}

  # create array
  temp_p = matrix(NA, n.samples, length(temp_data))

  # transform psi off logit-scale back to probability scale
  temp_p <- plogis(logit_temp)

  # calculate means and credible intervals
  temp_means = colMeans(temp_p )
  temp_CIs <- apply(temp_p ,2,quantile, c(0.025,0.975), na.rm=TRUE)

  # stuff into df
  temp_preds <- data.frame(predicted = temp_means,
                             cov_zsc = temp_data,
                             LCI = temp_CIs[1,],
                             UCI = temp_CIs[2,])

  # add back-transformed values to df (just in case i want them later)
  temp_mean  <- mean(dat$temp, na.rm = TRUE)
  temp_sd    <- sd(dat$temp, na.rm = TRUE)
  temp_preds$cov_value  <- round(temp_preds$cov_zsc  * temp_sd  + temp_mean, 2)

  # add species
  temp_preds_e <- temp_preds
  temp_preds_e$species <- "ENES"
  
  
  
##### Detection probability ENES ------------------------------------------
  
  # posterior samples of intercept
  alpha0_e <- E2[, "alpha0"]
  
  # convert from logit scale to probability
  p_avg_e <- plogis(alpha0_e)
  
  mean(p_avg_e) # posterior mean detection probability
  quantile(p_avg_e, c(0.025, 0.975)) # 95% credible interval

  # 0.1610015
  # 2.5%     97.5% 
  # 0.1277433 0.1980939
  
##### Marginal Psi Predictions Dataframe - OSS -------------------------------------------  
  
  b <- O2 
  dat <- o_covs 
  n.samples = nrow(b) # number of posterior samples
  
  # Latitude on psi
  
  # range of lat
  r<- range(lat.2D)
  dim(lat.2D)
  
  # create sequence along range
  lat_data <- seq(r[1], r[2], length.out=50) 
  
  # set all other covs at mean
  HB = 0
  HU = 0
  BS = 0
  BU = 0
  long = 0
  elev = 0
  
  # create matrices to stick estimates in
  logit_lat_psi = matrix(NA, n.samples, length(lat_data))
  
  # psi predictions for enes
  # Sample from posterior for the sequence of values for cov - raw linear predictor for each MCMC sample 
  for (i in 1:n.samples){
    for (j in 1:length(lat_data)){
      # create the linear predictors for dominant species
      logit_lat_psi[i,j] = 
        b[,'beta0.psi'][[i]] + 
        b[,'beta1.psi.BU'][[i]] * BU + 
        b[,'beta2.psi.HB'][[i]] * HB + 
        b[,'beta3.psi.HU'][[i]] * HU + 
        b[,'beta4.psi.BS'][[i]] * BS + 
        b[,'beta5.psi.lat'][[i]] * lat_data[j] +  
        b[,'beta6.psi.lon'][[i]] * long +  
        b[,'beta8.psi.elev'][[i]] * elev 
    }}
  
  # create array 
  lat_psi = matrix(NA, n.samples, length(lat_data))
  
  # transform psi off logit-scale back to probability scale - for each sample
  lat_psi <- plogis(logit_lat_psi)
  
  # calculate means and credible intervals - mean predicted occupancy probability across all posterior samples
  lat_psi_means = colMeans(lat_psi) # - taking the mean of all samples (transformed occu prob)
  lat_psi_CIs <- apply(lat_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)
  
  # stuff into df
  lat_psi_preds <- data.frame(predicted = lat_psi_means, 
                              cov_zsc = lat_data,
                              LCI = lat_psi_CIs[1,],
                              UCI = lat_psi_CIs[2,])
  
  # add back-transformed values to df (just in case i want them later)
  lat_mean  <- mean(dat$lat, na.rm = TRUE)
  lat_sd    <- sd(dat$lat, na.rm = TRUE)
  lat_psi_preds$cov_value  <- lat_psi_preds$cov_zsc  * lat_sd  + lat_mean
  
  # add species column
  lat_psi_preds_o = lat_psi_preds
  lat_psi_preds_o$species = "OSS"
  
  # Longitude on psi
  
  # range
  r<- range(lon.2D)
  # create sequence along range
  lon_data <- seq(r[1], r[2], length.out=50) 
  
  # set all other covs at mean
  HB = 0
  HU = 0
  BS = 0
  BU = 0
  lat = 0
  elev = 0
  
  # create matrices to stick estimates in
  logit_lon_psi = matrix(NA, n.samples, length(lon_data))
  
  # psi predictions for enes
  # Sample from posterior for the sequence of values for cov
  for (i in 1:n.samples){
    for (j in 1:length(lon_data)){
      logit_lon_psi[i,j] = 
        b[,'beta0.psi'][[i]] + 
        b[,'beta1.psi.BU'][[i]] * BU + 
        b[,'beta2.psi.HB'][[i]] * HB + 
        b[,'beta3.psi.HU'][[i]] * HU + 
        b[,'beta4.psi.BS'][[i]] * BS + 
        b[,'beta5.psi.lat'][[i]] * lat +  
        b[,'beta6.psi.lon'][[i]] * lon_data[j] +  
        b[,'beta8.psi.elev'][[i]] * elev 
    }}
  
  # create array 
  lon_psi = matrix(NA, n.samples, length(lon_data))
  
  # transform psi off logit-scale back to probability scale
  lon_psi <- plogis(logit_lon_psi)
  
  # calculate means and credible intervals 
  lon_psi_means = colMeans(lon_psi) 
  lon_psi_CIs <- apply(lon_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)
  
  # stuff into df
  lon_psi_preds <- data.frame(predicted = lon_psi_means, 
                              cov_zsc = lon_data,
                              LCI = lon_psi_CIs[1,],
                              UCI = lon_psi_CIs[2,])
  
  # add back-transformed values to df (just in case i want them later)
  long_mean  <- mean(dat$long, na.rm = TRUE)
  long_sd    <- sd(dat$long, na.rm = TRUE)
  lon_psi_preds$cov_value  <- lon_psi_preds$cov_zsc  * long_sd  + long_mean
  
  # add species column
  lon_psi_preds_o = lon_psi_preds
  lon_psi_preds_o$species = "OSS"
  
  
# Elevation on psi
  
  # range
  r<- range(elev.2D)
  # create sequence along range
  elev_data <- seq(r[1], r[2], length.out=50) 
  
  # set all other covs at mean
  HB = 0
  HU = 0
  BS = 0
  BU = 0
  lat = 0
  long = 0
  
  # create matrices to stick estimates in
  logit_elev_psi = matrix(NA, n.samples, length(elev_data))
  
  # psi predictions for enes
  # Sample from posterior for the sequence of values for cov 
  for (i in 1:n.samples){
    for (j in 1:length(elev_data)){
      logit_elev_psi[i,j] = 
        b[,'beta0.psi'][[i]] + 
        b[,'beta1.psi.BU'][[i]] * BU + 
        b[,'beta2.psi.HB'][[i]] * HB + 
        b[,'beta3.psi.HU'][[i]] * HU + 
        b[,'beta4.psi.BS'][[i]] * BS + 
        b[,'beta5.psi.lat'][[i]] * lat +  
        b[,'beta6.psi.lon'][[i]] * long +  
        b[,'beta8.psi.elev'][[i]] * elev_data[j] 
    }}
  
  # create array 
  elev_psi = matrix(NA, n.samples, length(elev_data))
  
  # transform psi off logit-scale back to probability scale
  elev_psi <- plogis(logit_elev_psi)
  
  # calculate means and credible intervals
  elev_psi_means = colMeans(elev_psi) 
  elev_psi_CIs <- apply(elev_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)
  
  # stuff into df
  elev_psi_preds <- data.frame(predicted = elev_psi_means, 
                               cov_zsc = elev_data,
                               LCI = elev_psi_CIs[1,],
                               UCI = elev_psi_CIs[2,])
  
  # add back-transformed values to df (just in case i want them later)
  elev_mean  <- mean(dat$elev, na.rm = TRUE)
  elev_sd    <- sd(dat$elev, na.rm = TRUE)
  elev_psi_preds$cov_value  <- round(elev_psi_preds$cov_zsc  * elev_sd  + elev_mean)
  
  # add species column
  elev_psi_preds_o = elev_psi_preds
  elev_psi_preds_o$species = "OSS"

  
##### Theta Prediction Dataframe - OSS ------------------------------------------
  
  b <- O2    
  dat <- o_covs 
  
# Downed Wood effect on Plot Use
  
  r<- range(downedwood.3D)
  # create sequence along range
  DW_data <- seq(r[1], r[2], length.out=50) 
  
  # create matrices to stick estimates in
  logit_DW_theta = matrix(NA, n.samples, length(DW_data))
  
  # theta predictions for enes
  # Sample from posterior for the sequence of values for cov 
  for (i in 1:n.samples){
    for (j in 1:length(DW_data)){
      logit_DW_theta[i,j] = 
        b[,'beta0.theta'][[i]] + 
        b[,'beta1.theta.DW'][[i]] * DW_data[j] 
    }}
  
  # create array 
  DW_theta = matrix(NA, n.samples, length(DW_data))
  
  # transform psi off logit-scale back to probability scale
  DW_theta <- plogis(logit_DW_theta)
  
  # calculate means and credible intervals
  DW_means = colMeans(DW_theta) 
  DW_CIs <- apply(DW_theta,2,quantile, c(0.025,0.975), na.rm=TRUE)
  
  # stuff into df
  DW_preds <- data.frame(predicted = DW_means, 
                             cov_zsc = DW_data,
                             LCI = DW_CIs[1,],
                             UCI = DW_CIs[2,])
  
  # add back-transformed values to df (just in case i want them later)
  dw_mean  <- mean(dat$DW, na.rm = TRUE)
  dw_sd    <- sd(dat$DW, na.rm = TRUE)
  DW_preds$cov_value  <- round(DW_preds$cov_zsc  * dw_sd  + dw_mean, 2)
  
  # add species
  DW_preds_o <- DW_preds
  DW_preds_o$species <- "OSS"
  
  
##### Detection Prediction Dataframe - OSS --------------------------------------
  
  b <- O2    
  dat <- o_covs 
  
  r<- range(temp.3D)
  # create sequence along range
  temp_data <- seq(r[1], r[2], length.out=50) 
  
  # create matrices to stick estimates in
  logit_temp = matrix(NA, n.samples, length(temp_data))
  
  # theta predictions for enes
  # Sample from posterior for the sequence of values for cov 
  for (i in 1:n.samples){
    for (j in 1:length(temp_data)){
      logit_temp[i,j] = 
        # alpha0_mean[i] +
        b[,'alpha0'][[i]] + 
        b[,'alpha1'][[i]] * temp_data[j] + 
        b[,'alpha2'][[i]] * temp_data[j]^2 
    }}
  
  # create array 
  temp_p = matrix(NA, n.samples, length(temp_data))
  
  # transform psi off logit-scale back to probability scale
  temp_p <- plogis(logit_temp)
  
  # calculate means and credible intervals
  temp_means = colMeans(temp_p ) 
  temp_CIs <- apply(temp_p ,2,quantile, c(0.025,0.975), na.rm=TRUE)
  
  # stuff into df
  temp_preds <- data.frame(predicted = temp_means, 
                           cov_zsc = temp_data,
                           LCI = temp_CIs[1,],
                           UCI = temp_CIs[2,])
  
  # add back-transformed values to df (just in case i want them later)
  temp_mean  <- mean(dat$temp, na.rm = TRUE)
  temp_sd    <- sd(dat$temp, na.rm = TRUE)
  temp_preds$cov_value  <- round(temp_preds$cov_zsc  * temp_sd  + temp_mean, 2)
  
  # add species
  temp_preds_o <- temp_preds
  temp_preds_o$species <- "OSS"
  
##### Detection probability OSS ------------------------------------------
  
  # posterior samples of intercept
  alpha0_o <- O2[, "alpha0"]
  
  # convert from logit scale to probability
  p_avg_o <- plogis(alpha0_o)
  
  mean(p_avg_o) # posterior mean detection probability
  quantile(p_avg_o, c(0.025, 0.975)) # 95% credible interval
  
  # 0.2517736
  # 2.5%     97.5% 
  # 0.1917336 0.3315254 
  
  
##### Plot ----------------------------------------------------------------------
  
# After creating all prediction dataframes for both species, combine them:

  # Latitude predictions
  lat_all <- rbind(
    lat_psi_preds_e %>% mutate(covariate = "Latitude", parameter = "ψ"),
    lat_psi_preds_o %>% mutate(covariate = "Latitude", parameter = "ψ")
  )
  
  # Longitude predictions
  lon_all <- rbind(
    lon_psi_preds_e %>% mutate(covariate = "Longitude", parameter = "ψ"),
    lon_psi_preds_o %>% mutate(covariate = "Longitude", parameter = "ψ")
  )
  
  # Elevation predictions
  elev_all <- rbind(
    elev_psi_preds_e %>% mutate(covariate = "Elevation", parameter = "ψ"),
    elev_psi_preds_o %>% mutate(covariate = "Elevation", parameter = "ψ")
  )
  
  # Downed wood predictions
  dwd_all <- rbind(
    DW_preds_e %>% mutate(covariate = "Downed Wood", parameter = "θ"),
    DW_preds_o %>% mutate(covariate = "Downed Wood", parameter = "θ")
  )
  
  # Temperature predictions
  temp_all <- rbind(
    temp_preds_e %>% mutate(covariate = "Temperature", parameter = "p"),
    temp_preds_o %>% mutate(covariate = "Temperature", parameter = "p")
  )
  
  # Combine into one dataframe
  all_preds <- bind_rows(lat_all, lon_all, elev_all, dwd_all, temp_all)
  
  # Set factor levels to control ordering
  all_preds$covariate <- factor(all_preds$covariate, 
                                levels = c("Latitude", "Longitude", "Elevation", 
                                           "Downed Wood", "Temperature"))
  
  # Rename species labels
  all_preds$species <- factor(all_preds$species,
                              levels = c("ENES", "OSS"),
                              labels = c("Ensatina", "Oregon Slender"))
  
  # Create y-axis label that varies by parameter
  all_preds <- all_preds %>%
    mutate(y_label = case_when(
      parameter == "ψ" ~ "Occupancy (ψ)",
      parameter == "θ" ~ "Subplot Use (θ)",
      parameter == "p" ~ "Detection (p)"
    ))

  
  
#### Plot
  
  # Define color palette
  species_colors <- c("Ensatina" = "#6091C2", 
                      "Oregon Slender" = "gray20")
  
  # Reorder species factor
  all_preds$species <- factor(all_preds$species, levels = c("Oregon Slender", "Ensatina"))
  
  # Separate the data by parameter type
  psi_data <- all_preds %>% filter(parameter == "ψ")
  theta_data <- all_preds %>% filter(parameter == "θ")
  p_data <- all_preds %>% filter(parameter == "p")
  
  # Create occupancy (ψ) plot - 3 columns
  p_psi <- ggplot(psi_data, aes(x = cov_value, y = predicted,
                                color = species, fill = species)) +
    geom_line(linewidth = 0.8) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.2, color = NA) +
    facet_grid(species ~ covariate, scales = "free_x",
               labeller = labeller(covariate = c("Latitude" = "Latitude (°N)",
                                                 "Longitude" = "Longitude (°W)",
                                                 "Elevation" = "Elevation (m)"))) +
    scale_color_manual(values = species_colors) +
    scale_fill_manual(values = species_colors) +
    facetted_pos_scales(
      x = list(
        covariate == "Latitude" ~ scale_x_continuous(breaks = c(44.6, 45.1)),
        covariate == "Longitude" ~ scale_x_continuous(breaks = c(-122.6, -122.3)),
        covariate == "Elevation" ~ scale_x_continuous(breaks = c(750, 2000, 3200))
      )
    ) +    labs(y = expression(bold("Occupancy ("*psi*")")), x = "") +
    theme_classic() +
    theme(
      strip.text.x = element_text(size = 10, face = "bold"),
      strip.text.y = element_blank(),  # Remove species labels
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 10, face = "bold"),
      panel.spacing.x = unit(0.3, "lines"),
      plot.margin = margin(t = 5, r = 5, b = 2, l = 5),
      legend.position = "none"
    )
  
  # Create subplot use (θ) plot - 1 column
  p_theta <- ggplot(theta_data, aes(x = cov_value, y = predicted,
                                    color = species, fill = species)) +
    geom_line(linewidth = 0.8) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.2, color = NA) +
    facet_grid(species ~ covariate, scales = "free",
           labeller = labeller(covariate = c("Downed Wood" = "Downed Wood (count)"))) +
    scale_color_manual(values = species_colors) +
    scale_fill_manual(values = species_colors) +
    scale_x_continuous(n.breaks = 5) +  
    labs(y = expression(bold("Subplot Use ("*theta*")")), x = "") +
    theme_classic() +
    theme(
      strip.text.x = element_text(size = 10, face = "bold"),
      strip.text.y = element_blank(),  # Remove species labels 
      strip.background.x = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.background.y = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 10, face = "bold"),
      panel.spacing.x = unit(0.3, "lines"),
      plot.margin = margin(t = 5, r = 5, b = 2, l = 5),
      legend.position = "none"
    )
  
  # Create detection (p) plot - 1 column
  p_detection <- ggplot(p_data, aes(x = cov_value, y = predicted,
                                    color = species, fill = species)) +
    geom_line(linewidth = 0.8) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.2, color = NA) +
    facet_grid(species ~ covariate, scales = "free",
               labeller = labeller(covariate = c("Temperature" = "Temperature (°C)"))) +
    scale_color_manual(values = species_colors) +
    scale_fill_manual(values = species_colors) +
    scale_x_continuous(n.breaks = 5) +  
    #scale_x_continuous(breaks = c(0, 15, 30)) +  
    labs(y = expression(bold("Detection ("*italic(p)*")")), x = "") +
    theme_classic() +
    theme(
      strip.text.x = element_text(size = 10, face = "bold"),
      strip.text.y = element_blank(),
      strip.background.x = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.background.y = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 10, face = "bold"),
      panel.spacing.x = unit(0.3, "lines"),
      plot.margin = margin(t = 5, r = 5, b = 2, l = 5),
      legend.position = "bottom",
      legend.title = element_blank(),  
      legend.text = element_text(size = 10),
      legend.key.spacing.x = unit(3, "cm")
    )
  
  
  
  # Combine with patchwork
  
  p_final <- p_psi + p_theta + p_detection + 
    #plot_layout(widths = c(3, 1, 1), nrow = 1)
    plot_layout(ncol=1)
  
  p_final
  
  # Save
  ggsave("figures/marginal-effects-combined-2.png", 
         plot = p_final,
         width = 14, 
         height = 22,
         units = "cm",
         dpi = 600)
  

  
  
  
  