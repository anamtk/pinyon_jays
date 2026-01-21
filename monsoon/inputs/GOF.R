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

# Update for residuals ----------------------------------------------------

parms <- c("resid")

mod2 <- update(mod,
               parameters.to.save = parms,
               n.iter = 3000)

sum2 <- summary(mod2$samples)

saveRDS(sum2, file ="/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model_residuals.RDS")

# Replicated data ---------------------------------------------------------

parms2 <- c('ebird.count.rep')

mod3 <- update(mod,
               parameters.to.save = parms2,
               n.iter = 3000)

sum3 <- summary(mod3$samples)

saveRDS(sum3, file ="/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model_yrep.RDS")

(end.time <- Sys.time())


# samples of yrep ---------------------------------------------------------

mod4 <- update(mod,
               parameters.to.save = parms2,
               n.iter = 350)

saveRDS(mod4$samples, file ="/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model_yrepsamples.RDS")

(end.time <- Sys.time())


# Samples of covariate effects --------------------------------------------

parms3 <- c('a0',
            'a',
            'wA',
            'wB',
            'wC',
            'c0',
            'c')

mod5 <- update(mod,
               parameters.to.save = parms3,
               n.iter = 350)

saveRDS(mod5$samples, file ="/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model_covariate_effect_samples.RDS")

(end.time <- Sys.time())


# RMSE --------------------------------------------------------------------


mod6 <- update(mod,
               parameters.to.save = 'RMSE',
               n.iter = 3000)

saveRDS(mod6$samples, file ="/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model_RMSE_samples.RDS")

(end.time <- Sys.time())
