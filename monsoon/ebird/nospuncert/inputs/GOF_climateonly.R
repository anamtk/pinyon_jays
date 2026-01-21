# Load packages -----------------------------------------------------------

(start.time <- Sys.time())

package.list <- c("jagsUI", "coda",
                  'dplyr', 'stringr',
                  'magrittr', 'tidyr',
                  'mcmcplots','ggplot2') 

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load model --------------------------------------------------------------

mod <- readRDS("/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model2_checkN_climateonly.RDS")

# samples of yrep ---------------------------------------------------------

mod4 <- update(mod,
               parameters.to.save = parms2,
               n.iter = 350)

saveRDS(mod4$samples, file ="/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model_yrepsamples_climateonly.RDS")

(end.time <- Sys.time())

# Replicated data ---------------------------------------------------------

parms2 <- c('ebird.count.rep')

mod3 <- update(mod,
               parameters.to.save = parms2,
               n.iter = 3000)

sum3 <- summary(mod3$samples)

saveRDS(sum3, file ="/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model_yrep_climateonly.RDS")

(end.time <- Sys.time())


