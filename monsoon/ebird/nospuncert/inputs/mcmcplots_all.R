# Load packages -----------------------------------------------------------

(start.time <- Sys.time())

package.list <- c("jagsUI", "coda",
                  'dplyr', 'stringr',
                  'magrittr', 'tidyr',
                  'tibble', 'purrr',
                  'mcmcplots','ggplot2') 

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load models -------------------------------------------------------------

check1 <- readRDS("/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model_checkN.RDS")

check2 <- readRDS("/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model2_checkN.RDS")

blob1 <- readRDS("/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model_blobN.RDS")

blob2 <- readRDS("/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model2_blobN.RDS")

# Diagnose models ----------------------------------------------------------

mcmcplot(check1$samples,
         dir = "/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/mcmcplots/checkN")

mcmcplot(check2$samples,
         dir = "/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/mcmcplots/checkN_run2")

mcmcplot(blob1$samples,
         dir = "/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/mcmcplots/blobN")

mcmcplot(blob2$samples,
         dir = "/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/mcmcplots/blobN_run2")




