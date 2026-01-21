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

# Load data ---------------------------------------------------------------

data <- readRDS("/scratch/atm234/pinyon_jays/ebird/nospuncert/inputs/ebird_data_list_nospuncert.RDS")

# Define data list --------------------------------------------------------

data_list <- list(#latent N loop:
  n.blobs = data$n.blobs,
  n.years = data$n.years,
  n.lag = data$n.lag,
  n.clag = data$n.clag,
  Cone = data$Cone,
  Temp = data$Temp,
  PPT = data$PPT,
  Monsoon = data$Monsoon,
  PinyonBA = data$PinyonBA,
  blobArea = data$blobArea,
  #ebird loop
  n.ebird.check = data$n.ebird.check,
  ebird.count = data$ebird.count,
  SurveyType = data$SurveyType,
  StartTime = data$StartTime,
  Duration = data$Duration,
  Distance = data$Distance,
  Speed = data$Speed,
  NumObservers = data$NumObservers,
  listArea = data$listArea,
  n.checklists = data$n.checklists,
  deltaA = rep(1, data$n.lag),
  deltaB = rep(1, data$n.clag),
  deltaC = rep(1, data$n.clag))

# Parameters to save ------------------------------------------------------

parameters <- c('a0',
                'a',
                'wA',
                'wB',
                'wC',
                'b0',
                'b',
                'c0',
                'c1',
                'c')

# Initials ----------------------------------------------------------------

inits_list <- readRDS('/scratch/atm234/pinyon_jays/ebird/nospuncert/inputs/ebird_init_list_nospuncert.RDS')

# Run model ---------------------------------------------------------------

model <- jagsUI::jags(data = data_list,
                      parameters.to.save = parameters,
                      inits = inits_list,
                      model.file = '/scratch/atm234/pinyon_jays/ebird/nospuncert/inputs/ebird_abund_JAGS_checkN_climateonly.R',
                      parallel = TRUE,
                      n.chains = 3,
                      n.iter = 25000,
                      n.burnin = 5000,
                      n.thin = 10,
                      DIC = TRUE)

#save as an R data object
saveRDS(model, 
        file ="/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model_checkN_climateonly.RDS")

(end.time <- Sys.time())

(tot.time <- end.time - start.time)

# Diagnose model ----------------------------------------------------------

mcmcplot(model$samples,
         dir = "/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/mcmcplots/checkN_climateonly")