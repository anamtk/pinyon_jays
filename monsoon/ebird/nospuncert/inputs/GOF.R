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

mod <- readRDS("/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model2.RDS")

# Update for residuals ----------------------------------------------------

parms <- c("resid")

mod2 <- update(mod,
               parameters.to.save = parms,
               n.iter = 3000)

sum2 <- summary(mod2$samples)

saveRDS(sum2, file ="/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model2_residuals.RDS")

# Replicated data ---------------------------------------------------------

parms2 <- c('ebird.count.rep')

mod3 <- update(mod,
               parameters.to.save = parms2,
               n.iter = 3000)

sum3 <- summary(mod3$samples)

saveRDS(sum3, file ="/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model2_yrep.RDS")

(end.time <- Sys.time())

