## =================================================
##
## Title: both-occu-preds-treatment-plot
## Author: Jasmine Williamson 
## Date Created: 07/22/2025
##
## Description: Makes plot that shows occupancy predictions by treatment for
## both species, on avg using all data (NOT by year)
##
## =================================================
setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/multi-scale-occu-oss")

##### Setup --------------------------------------------
library(ggplot2)
library(dplyr)
library(patchwork)

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


### the below code is really slow ###
### I made an updated version of it to get estimates in script # 7 in this folder ###


##### Predicted occupancy for each treatment - ENES ----------------------------
# This produces marginal estimates of occupancy for different treatment types


b <- E2 
# number of posterior samples
n.samples = nrow(b) # why is this such a weird number?

# BU    
# range of BU
r<- range(BU.new)

# create sequence along range of BU
BU_data <- seq(r[1], r[2], length.out=2) #is this not identical to r created above?

# set all other covs at mean
HB = 0
HU = 0
BS = 0
lat = 0
long = 0
elev = 0

# create matrices to stick estimates in
logit_psi = matrix(NA, n.samples, length(BU_data))

# psi predictions for enes
# Sample from posterior for the sequence of values for BU - raw linear predictor for each MCMC sample 
for (i in 1:n.samples){
  for (j in 1:length(BU_data)){
    # create the linear predictors for dominant species
    logit_psi[i,j] = b[,'beta0.psi'][[i]] + b[,'beta1.psi.BU'][[i]] * BU_data[j] + b[,'beta2.psi.HB'][[i]] * HB + b[,'beta3.psi.HU'][[i]] * HU
    + b[,'beta4.psi.BS'][[i]] * BS + b[,'beta5.psi.lat'][[i]] * lat +  b[,'beta6.psi.lon'][[i]] * long +  b[,'beta8.psi.elev'][[i]] * elev 
  }}

# create array 
BU_psi = matrix(NA, n.samples, length(BU_data))

# transform psi off logit-scale back to probability scale - for each sample
BU_psi <- plogis(logit_psi)

# calculate means and credible intervals - mean predicted occupancy probability across all posterior samples
BU_psi_means = colMeans(BU_psi) # - taking the mean of all samples (transformed occu prob)
BU_psi_CIs <- apply(BU_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)

# stuff into df
BU_psi_preds <- data.frame(yesno = BU_data, 
                           predicted = BU_psi_means, 
                           treatment = "BU",
                           LCI = BU_psi_CIs[1,],
                           UCI = BU_psi_CIs[2,])

# HB
r<- range(HB.new)
HB_data <- seq(r[1], r[2], length.out=2)

# set all other covs at mean
BU = 0
HU = 0
BS = 0
lat = 0
long = 0
elev = 0

# create matrices to stick estimates in
logit_HB_psi = matrix(NA, n.samples, length(HB_data))

# psi predictions for enes
# Sample from posterior for the sequence of values for HB
for (i in 1:n.samples){
  for (j in 1:length(HB_data)){
    # create the linear predictors for dominant species
    logit_psi[i,j] = b[,'beta0.psi'][[i]] + b[,'beta1.psi.BU'][[i]] * BU + b[,'beta2.psi.HB'][[i]] * HB_data[j] + b[,'beta3.psi.HU'][[i]] * HU
    + b[,'beta4.psi.BS'][[i]] * BS + b[,'beta5.psi.lat'][[i]] * lat +  b[,'beta6.psi.lon'][[i]] * long +  b[,'beta8.psi.elev'][[i]] * elev 
  }}

HB_psi = matrix(NA, n.samples, length(HB_data)) # create array 
HB_psi <- plogis(logit_psi) # transform psi off logit-scale back to probability scale
HB_psi_means = colMeans(HB_psi) # calculate means and credible intervals
HB_psi_CIs <- apply(HB_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)

HB_psi_preds <- data.frame(yesno = HB_data, 
                           predicted = HB_psi_means, 
                           treatment = "HB",
                           LCI = HB_psi_CIs[1,],
                           UCI = HB_psi_CIs[2,])

# HU   
r<- range(HU.new)
HU_data <- seq(r[1], r[2], length.out=2)

# set all other covs at mean
BU = 0
HB = 0
BS = 0
lat = 0
long = 0
elev = 0

# create matrices to stick estimates in
logit_HU_psi = matrix(NA, n.samples, length(HU_data))

# psi predictions for enes
# Sample from posterior for the sequence of values for HU 
for (i in 1:n.samples){
  for (j in 1:length(HU_data)){
    # create the linear predictors for dominant species
    logit_psi[i,j] = b[,'beta0.psi'][[i]] + b[,'beta1.psi.BU'][[i]] * BU + b[,'beta2.psi.HB'][[i]] * HB + b[,'beta3.psi.HU'][[i]] * HU_data[j]
    + b[,'beta4.psi.BS'][[i]] * BS + b[,'beta5.psi.lat'][[i]] * lat +  b[,'beta6.psi.lon'][[i]] * long +  b[,'beta8.psi.elev'][[i]] * elev 
  }}

HU_psi = matrix(NA, n.samples, length(HU_data))   # create array 
HU_psi <- plogis(logit_psi)   # transform psi off logit-scale back to probability scale
HU_psi_means = colMeans(HU_psi)   # calculate means and credible intervals
HU_psi_CIs <- apply(HU_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)

HU_psi_preds <- data.frame(yesno = HU_data, 
                           predicted = HU_psi_means, 
                           treatment = "HU",
                           LCI = HU_psi_CIs[1,],
                           UCI = HU_psi_CIs[2,])

# BS
r<- range(BS.new)
BS_data <- seq(r[1], r[2], length.out=2)

# set all other covs at mean
BU = 0
HB = 0
HU = 0
lat = 0
long = 0
elev = 0

# create matrices to stick estimates in
logit_BS_psi = matrix(NA, n.samples, length(BS_data))

# psi predictions for enes
# Sample from posterior for the sequence of values for BS 
for (i in 1:n.samples){
  for (j in 1:length(BS_data)){
    # create the linear predictors for dominant species
    logit_psi[i,j] = b[,'beta0.psi'][[i]] + b[,'beta1.psi.BU'][[i]] * BU + b[,'beta2.psi.HB'][[i]] * HB + b[,'beta3.psi.HU'][[i]] * HU +
      b[,'beta4.psi.BS'][[i]] * BS_data[j] + b[,'beta5.psi.lat'][[i]] * lat +  b[,'beta6.psi.lon'][[i]] * long +  b[,'beta8.psi.elev'][[i]] * elev 
  }}

BS_psi = matrix(NA, n.samples, length(BS_data))   # create array 
BS_psi <- plogis(logit_psi)   # transform psi off logit-scale back to probability scale
BS_psi_means = colMeans(BS_psi)  # calculate means and credible intervals
BS_psi_CIs <- apply(BS_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)

BS_psi_preds <- data.frame(yesno = BS_data, 
                           predicted = BS_psi_means, 
                           treatment = "BS",
                           LCI = BS_psi_CIs[1,],
                           UCI = BS_psi_CIs[2,])

# UU  
UU_psi_preds <- subset(BS_psi_preds, yesno == 0)
UU_psi_preds$treatment <- "UU"

# Combine
alltreatment_preds <- rbind(BS_psi_preds, HU_psi_preds, BU_psi_preds, HB_psi_preds)
alltreatment_preds_noUU <- subset(alltreatment_preds, yesno==1)
alltreatment_preds <- rbind(alltreatment_preds_noUU, UU_psi_preds)

# Plot

p <-  ggplot(alltreatment_preds, aes(x = treatment, y = predicted)) +
  geom_point(position = position_dodge(0.5), size = 1.5)+ 
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.1, position = position_dodge(0.5))+
  ylab(bquote("Predicted "*psi~""))+ xlab("Treatment Group")+
  labs(title="Predicted Occupancy Estimates by Treatment - ENES") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

#ggsave("figures/e-trt-psi-preds.png", plot = p, dpi = 500)  


# add species column for combined plot
alltreatment_preds$species <- "ENES"


##### Predicted occupancy for each treatment - OSS -----------------------------

b <- O2
n.samples = nrow(b)

# BU    
# range of BU
r<- range(BU.new)

# create sequence along range of BU
BU_data <- seq(r[1], r[2], length.out=2)

# set all other covs at mean
HB = 0
HU = 0
BS = 0
lat = 0
long = 0
elev = 0

# create matrices to stick estimates in
logit_psi = matrix(NA, n.samples, length(BU_data))

# psi predictions for enes
# Sample from posterior for the sequence of values for BU - raw linear predictor for each MCMC sample 
for (i in 1:n.samples){
  for (j in 1:length(BU_data)){
    # create the linear predictors for dominant species
    logit_psi[i,j] = b[,'beta0.psi'][[i]] + b[,'beta1.psi.BU'][[i]] * BU_data[j] + b[,'beta2.psi.HB'][[i]] * HB + b[,'beta3.psi.HU'][[i]] * HU
    + b[,'beta4.psi.BS'][[i]] * BS + b[,'beta5.psi.lat'][[i]] * lat +  b[,'beta6.psi.lon'][[i]] * long +  b[,'beta8.psi.elev'][[i]] * elev 
  }}

# create array 
BU_psi = matrix(NA, n.samples, length(BU_data))

# transform psi off logit-scale back to probability scale - for each sample
BU_psi <- plogis(logit_psi)

# calculate means and credible intervals - mean predicted occupancy probability across all posterior samples
BU_psi_means = colMeans(BU_psi) # - taking the mean of all samples (transformed occu prob)
BU_psi_CIs <- apply(BU_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)

# stuff into df
BU_psi_preds <- data.frame(yesno = BU_data, 
                           predicted = BU_psi_means, 
                           treatment = "BU",
                           LCI = BU_psi_CIs[1,],
                           UCI = BU_psi_CIs[2,])

# HB
r<- range(HB.new)
HB_data <- seq(r[1], r[2], length.out=2)

# set all other covs at mean
BU = 0
HU = 0
BS = 0
lat = 0
long = 0
elev = 0

# create matrices to stick estimates in
logit_HB_psi = matrix(NA, n.samples, length(HB_data))

# psi predictions for enes
# Sample from posterior for the sequence of values for HB
for (i in 1:n.samples){
  for (j in 1:length(HB_data)){
    # create the linear predictors for dominant species
    logit_psi[i,j] = b[,'beta0.psi'][[i]] + b[,'beta1.psi.BU'][[i]] * BU + b[,'beta2.psi.HB'][[i]] * HB_data[j] + b[,'beta3.psi.HU'][[i]] * HU
    + b[,'beta4.psi.BS'][[i]] * BS + b[,'beta5.psi.lat'][[i]] * lat +  b[,'beta6.psi.lon'][[i]] * long +  b[,'beta8.psi.elev'][[i]] * elev 
  }}

HB_psi = matrix(NA, n.samples, length(HB_data)) # create array 
HB_psi <- plogis(logit_psi) # transform psi off logit-scale back to probability scale
HB_psi_means = colMeans(HB_psi) # calculate means and credible intervals
HB_psi_CIs <- apply(HB_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)

HB_psi_preds <- data.frame(yesno = HB_data, 
                           predicted = HB_psi_means, 
                           treatment = "HB",
                           LCI = HB_psi_CIs[1,],
                           UCI = HB_psi_CIs[2,])

# HU   
r<- range(HU.new)
HU_data <- seq(r[1], r[2], length.out=2)

# set all other covs at mean
BU = 0
HB = 0
BS = 0
lat = 0
long = 0
elev = 0

# create matrices to stick estimates in
logit_HU_psi = matrix(NA, n.samples, length(HU_data))

# psi predictions for enes
# Sample from posterior for the sequence of values for HU 
for (i in 1:n.samples){
  for (j in 1:length(HU_data)){
    # create the linear predictors for dominant species
    logit_psi[i,j] = b[,'beta0.psi'][[i]] + b[,'beta1.psi.BU'][[i]] * BU + b[,'beta2.psi.HB'][[i]] * HB + b[,'beta3.psi.HU'][[i]] * HU_data[j]
    + b[,'beta4.psi.BS'][[i]] * BS + b[,'beta5.psi.lat'][[i]] * lat +  b[,'beta6.psi.lon'][[i]] * long +  b[,'beta8.psi.elev'][[i]] * elev 
  }}

HU_psi = matrix(NA, n.samples, length(HU_data))   # create array 
HU_psi <- plogis(logit_psi)   # transform psi off logit-scale back to probability scale
HU_psi_means = colMeans(HU_psi)   # calculate means and credible intervals
HU_psi_CIs <- apply(HU_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)

HU_psi_preds <- data.frame(yesno = HU_data, 
                           predicted = HU_psi_means, 
                           treatment = "HU",
                           LCI = HU_psi_CIs[1,],
                           UCI = HU_psi_CIs[2,])

# BS
r<- range(BS.new)
BS_data <- seq(r[1], r[2], length.out=2)

# set all other covs at mean
BU = 0
HB = 0
HU = 0
lat = 0
long = 0
elev = 0

# create matrices to stick estimates in
logit_BS_psi = matrix(NA, n.samples, length(BS_data))

# psi predictions for enes
# Sample from posterior for the sequence of values for BS 
for (i in 1:n.samples){
  for (j in 1:length(BS_data)){
    # create the linear predictors for dominant species
    logit_psi[i,j] = b[,'beta0.psi'][[i]] + b[,'beta1.psi.BU'][[i]] * BU + b[,'beta2.psi.HB'][[i]] * HB + b[,'beta3.psi.HU'][[i]] * HU +
      b[,'beta4.psi.BS'][[i]] * BS_data[j] + b[,'beta5.psi.lat'][[i]] * lat +  b[,'beta6.psi.lon'][[i]] * long +  b[,'beta8.psi.elev'][[i]] * elev 
  }}

BS_psi = matrix(NA, n.samples, length(BS_data))   # create array 
BS_psi <- plogis(logit_psi)   # transform psi off logit-scale back to probability scale
BS_psi_means = colMeans(BS_psi)  # calculate means and credible intervals
BS_psi_CIs <- apply(BS_psi,2,quantile, c(0.025,0.975), na.rm=TRUE)

BS_psi_preds <- data.frame(yesno = BS_data, 
                           predicted = BS_psi_means, 
                           treatment = "BS",
                           LCI = BS_psi_CIs[1,],
                           UCI = BS_psi_CIs[2,])

# UU  
UU_psi_preds <- subset(BS_psi_preds, yesno == 0)
UU_psi_preds$treatment <- "UU"

# Combine
alltreatment_preds2 <- rbind(BS_psi_preds, HU_psi_preds, BU_psi_preds, HB_psi_preds)
alltreatment_preds_noUU <- subset(alltreatment_preds2, yesno==1)
alltreatment_preds2 <- rbind(alltreatment_preds_noUU, UU_psi_preds)

# Plot

p2 <-  ggplot(alltreatment_preds2, aes(x = treatment, y = predicted)) +
  geom_point(position = position_dodge(0.5), size = 1.5)+ 
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.1, position = position_dodge(0.5))+
  ylab(bquote("Predicted "*psi~""))+ xlab("Treatment Group")+
  labs(title="Predicted Occupancy Estimates by Treatment - OSS") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

#ggsave("figures/o-trt-psi-preds.png", plot = p2, dpi = 500)  


# add species column for combined plot
alltreatment_preds2$species <- "OSS"



both_preds <- rbind(alltreatment_preds, alltreatment_preds2)
both_preds[, c(2,4:5)] <- round(both_preds[, c(2,4:5)], 3)
print(both_preds)


### SAVE ###
write.csv(both_preds, "posterior-preds-both-V3.csv", row.names = FALSE)



##### Combined Plot - predicted occu by trt for both spp -----------------------


# Set the order of treatments
both_preds$treatment <- factor(both_preds$treatment, levels = c("UU", "BU", "HU", "HB", "BS"))

p3 <- ggplot(both_preds, aes(x = treatment, y = predicted, color = species)) +
  geom_point(position = position_dodge(0.5), size = 4) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.1, size = 0.8, position = position_dodge(0.5)) +
  # scale_color_manual(values = c("BS" = "#EA7317", 
  #                               "HB" = "#FEC601",
  #                               "HU" = "#62B6CB", 
  #                               "UU" = "#2364AA",
  #                               "BU" = "#69995D")) + 
  scale_color_manual(values = c("ENES" = "#6091C2",
                                "OSS" = "gray20"),
                     labels = c("ENES" = "Ensatina",
                                "OSS" = "Oregon Slender salamander"))+ 
  scale_x_discrete(labels = c("UU" = "Control",
                              "BU" = "Burn-only",
                              "HU" = "Harvest",
                              "HB" = "Harvest-Burn",
                              "BS" = "Burn-Salvage")) +
  # scale_linetype_manual(values = c("ENES" = "solid", "OSS" = "42"),
  #                       labels = c("OSS" = expression(italic("B. wrighti")), 
  #                                  "ENES" = expression(italic("E. eschscholtzii"))
  #                                  )) +
  ylab(bquote("Predicted "*psi~"")) +
  #xlab("Treatment Group") +
  labs(title = "Predicted Occupancy by Treatment and Species") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "inside",
    legend.position.inside = c(0.25, 0.13),
    #legend.text = element_text(size = 16),
    #legend.title = element_text(size = 16),
    strip.text = element_text(size = 15, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(), 
    panel.border = element_rect(color = "gray40", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(color = "gray40", linewidth = 0.5),
    axis.ticks.length = unit(0.1, "cm")
  ) #+
#guides(color = "none") # removes trt legend

p3

ggsave("figures/both-spp-trt-psi-preds.png", 
       plot = p3, 
       width = 8,
       height = 5,
       bg = "white",
       dpi = 500)

