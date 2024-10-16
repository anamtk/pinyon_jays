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

data <- readRDS("/scratch/atm234/pinyon_jays/nospuncert/inputs/bbs_ebird_joint_data_list_nospuncert.RDS")

# Define data list --------------------------------------------------------

data_list <- list(n.grids = data$n.grids,
                  n.years = data$n.years,
                  n.lag = data$n.lag,
                  n.clag = data$n.clag,
                  Cone = data$Cone,
                  Temp = data$Temp,
                  PPT = data$PPT,
                  Monsoon = data$Monsoon,
                  PinyonBA = data$PinyonBA,
                  #BBS loop
                  n.bbs.years = data$n.bbs.years,
                  n.bbs.trans = data$n.bbs.trans,
                  ObserverExp = data$ObserverExp,
                  n.bbs.points = data$n.bbs.points,
                  bbs.count = data$bbs.count,
                  bbs.grid = data$bbs.grid,
                  bbs.year = data$bbs.year,
                  #ebird loop
                  n.ebird.pairs = data$n.ebird.pairs,
                  n.ebird.check = data$n.ebird.check,
                  ebird.count = data$ebird.count,
                  ebird.grid = data$ebird.grid,
                  SurveyType = data$SurveyType,
                  StartTime = data$StartTime,
                  Duration = data$Duration,
                  Distance = data$Distance,
                  NumObservers = data$NumObservers)

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

inits_list <- readRDS('/scratch/atm234/pinyon_jays/nospuncert/inputs/bbs_ebird_joint_init_list_nospuncert.RDS')

# Run model ---------------------------------------------------------------

model <- jagsUI::jags(data = data_list,
                      parameters.to.save = parameters,
                      inits = inits_list,
                      model.file = '/scratch/atm234/pinyon_jays/nospuncert/inputs/ebird_bbs_joint_abund_JAGS_nospuncert.R',
                      parallel = TRUE,
                      n.chains = 3,
                      n.iter = 1,
                      DIC = TRUE)

#save as an R data object
saveRDS(model, 
        file ="/scratch/atm234/pinyon_jays/nospuncert/outputs/ebird_bbs_joint_abund_model.RDS")

(end.time <- Sys.time())

(tot.time <- end.time - start.time)

# Diagnose model ----------------------------------------------------------

mcmcplot(model$samples,
         dir = "/scratch/atm234/pinyon_jays/nospuncert/outputs/mcmcplots/test")