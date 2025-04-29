################################################################################
## This code creates summary plots of model of tree radial growth

## Code by Kyle C. Rodman, Ecological Restoration Institute.

# 4/24/2023; updated on 5/10/2024
################################################################################

### Bring in necessary packages
package.list <- c("here", "tidyverse", "pals", "ggpubr", "FNN", "ggnewscale")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}

################################################################################

### Read in files

## Read in coefficient values

growthCoefs <- readRDS(here("ModelOutputs", "growthv3_bayesianFinal", "Coefs_full.rds"))

## Get growth data that was used for analysis

growthData <- read.csv(here("Data", "AnalysisReady", "growthData_merged.csv"))


################################################################################

## Creating the antecedent climate terms for plotting and prediction. Weighted

# sums based on the fitted weight parameters

## Soil moisture

wtsCoef <- growthCoefs$quantiles[str_detect(rownames(growthCoefs$quantiles), "wASW"),3]

temp <- growthData %>%
  dplyr::select(contains("SWA")) %>%
  as.matrix()

growthData$aswVals <- apply(temp, MARGIN = 1, FUN = function(x){sum(x*wtsCoef)})

## VPD

wtsCoef <- growthCoefs$quantiles[str_detect(rownames(growthCoefs$quantiles), "wVPD"),3]

temp <- growthData %>%
  dplyr::select(contains("VPD")) %>%
  as.matrix()

growthData$vpdVals <- apply(temp, MARGIN = 1, FUN = function(x){sum(x*wtsCoef)})

################################################################################

### Make plots of main effects, interactions, and antecedent weights

## Format coefficient estimates as DF

coefs <- growthCoefs$quantiles %>%
  as.data.frame() %>%
  dplyr::select(c('2.5%', "50%", "97.5%")) %>%
  rownames_to_column() %>%
  filter(!str_detect(as.character(rowname), "bSiteTrt") &
           !str_detect(as.character(rowname), "bSite") &
           !str_detect(as.character(rowname), "bTrt") &
           !str_detect(as.character(rowname), "b0") &
           !str_detect(as.character(rowname), "bAR1") &
           !str_detect(as.character(rowname), "sig") &
           !str_detect(as.character(rowname), "r2") &
           !str_detect(as.character(rowname), "wASW") &
           !str_detect(as.character(rowname), "wVPD") &
           !str_detect(as.character(rowname), "diff") &
           !str_detect(as.character(rowname), "deviance")) %>%
  rename(parameter = rowname, lower = "2.5%", median = "50%", upper = "97.5%") %>%
  mutate(Treatment = ifelse(str_sub(parameter, start = -2, end = -2) == "1", "Untreated", "Treated"),
         parameter = str_sub(parameter, start = 2, end = -4)) %>%
  mutate(parameter = factor(case_when(
    parameter == "BA" ~ "Basal Area\nW/In Treatment",
    parameter == "DBH" ~ "Diam. at Breast\nHeight (DBH)",
    parameter == "ASW" ~ "Available Soil\nWater (ASW)",
    parameter == "VPD" ~ "Vapor Pressure\nDeficit (VPD)",
    parameter == "ASWbyVPD" ~ "ASW by\nVPD",
    parameter == "DBHbyASW" ~ "DBH by\nASW",
    parameter == "DBHbyVPD" ~ "DBH by\nVPD"
  ), levels = rev(c("Basal Area\nW/In Treatment", "Diam. at Breast\nHeight (DBH)",
                    "Available Soil\nWater (ASW)", "Vapor Pressure\nDeficit (VPD)",
                    "ASW by\nVPD", "DBH by\nASW", "DBH by\nVPD"))))

## Getting coefficients with any antecedent term

coefs2 <- coefs[str_detect(coefs$parameter, "VPD")|str_detect(coefs$parameter, "ASW"),c(1,3,5)]

coefs2[1:2,2] <- coefs2[1:2,2]*sd(growthData$aswVals) # Adjust by SD of ASW values

coefs2[3:4,2] <- coefs2[3:4,2]*sd(growthData$vpdVals) # Adjust by SD of VPD values

coefs2[5:6,2] <- coefs2[5:6,2]*sd(growthData$vpdVals)*sd(growthData$aswVals)  # Adjust by SD of both values since it's an interaction

coefs2[7:8,2] <- coefs2[7:8,2]*sd(growthData$aswVals) # Adjust by SD of ASW values

coefs2[9:10,2] <- coefs2[9:10,2]*sd(growthData$vpdVals) # Adjust by SD of VPD values

## Create coefficient plot

dodge <- position_dodge(width=0.5)

(a <- ggplot(coefs, aes(x = median, y = parameter, color = Treatment, group = Treatment)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5) +
    geom_point(size = 2.5, position = dodge) +
    theme_bw() +
    xlab("Covariate Effect\n(posterior median and 95% CI)") + xlim(c(-1.1, 1.1)) +
    geom_linerange(aes(xmin = lower, xmax = upper), position=dodge, linewidth = 1) +
    scale_color_manual(values = c("#7fbf7b", "#af8dc3")) +
    geom_point(data = coefs2, size = 2.5, color = "grey80", shape = 21,
               alpha = 0.7, position = dodge, stroke = 1) +
    theme(legend.position = "top",
          axis.title.y=element_blank()))

