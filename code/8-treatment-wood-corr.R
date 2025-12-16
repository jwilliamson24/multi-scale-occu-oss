## =================================================
##
## Title: treatment-wood-corr
## Author: Jasmine Williamson 
## Date Created: 10/01/2025
##
## Description: Checking if DWD is correlated to treatment
##
## =================================================
setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/multi-scale-occu-oss")


## Setup 
library(dplyr)

# Load data
source("data/attach.nimble_v2 copy.R")
load("data/msc-enes-data-workspace-V2.RData")



# Flatten data to vectors
  dw_vec <- c(downedwood.new)
  bs_vec <- c(BS.new)
  hu_vec <- c(HU.new)
  hb_vec <- c(HB.new)
  bu_vec <- c(BU.new)

# Create treatment category variable
  treatment <- ifelse(bs_vec == 1, "Salvage",
                      ifelse(hb_vec == 1, "Harvest-Burn",
                             ifelse(hu_vec == 1, "Harvest",
                                    ifelse(bu_vec == 1, "Burn", "Control"))))

# Convert to factor for ordering
  treatment <- factor(treatment, levels=c("Control", "Burn", "Harvest", 
                                          "Harvest-Burn", "Salvage"))




# 1. Summary statistics by treatment
  aggregate(dw_vec ~ treatment, FUN=function(x) c(mean=mean(x), 
                                                  sd=sd(x), 
                                                  n=length(x)))


# 2. Boxplot across all treatments
  boxplot(dw_vec ~ treatment, 
          xlab="Treatment", 
          ylab="Downed Wood (standardized)",
          main="Downed Wood by Treatment Category",
          las=2)  # rotate x-axis labels

# 3. Statistical test - does DW differ among treatments?
  kruskal.test(dw_vec ~ treatment)  # non-parametric
  anova(lm(dw_vec ~ treatment))  # parametric

# 4. Post-hoc pairwise comparisons (if overall test significant)
  pairwise.wilcox.test(dw_vec, treatment, p.adjust.method = "bonferroni")

# 5. Table of means
  tapply(dw_vec, treatment, mean)








