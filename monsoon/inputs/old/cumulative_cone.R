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

mod <- readRDS("/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model2_checkN.RDS")

# update for cumulative cone weight ---------------------------------------

mod2 <- update(mod,
               parameters.to.save = 'cum.cone.wt',
               n.iter = 3000)

sum <- summary(mod2$samples)

saveRDS(sum, file ="/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model_cumulative_cone_wt.RDS")

(end.time <- Sys.time())
