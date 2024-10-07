# Testing the covariate paraemter sampler
# Ana Miller-ter Kuile
# March 2, 2023

# this script tests the sampler Kiona generated that lets me
# sample from the posterior samples of each covariate parameter
# from statistical models for the full IPM

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", 
                  "tidyverse", 
                  "jagsUI")


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Import data -------------------------------------------------------------

data_list <- readRDS(here('data',
                          'jags_input_data',
                          'bbs_ebird_joint_data_list.RDS'))

# Model path --------------------------------------------------------------


model_file <- (here('code',
                    'jags_models',
                    'grid_sampler_ebirdbbs.R'))


# Parameters to save ------------------------------------------------------


parameters <- c('bbs.grid.index',
                'bbs.grid',
                'ebird.grid.index',
                'ebird.grid')

# Run model ---------------------------------------------------------------


model <- jagsUI::jags(data = data_list,
                      parameters.to.save = parameters,
                      inits = NULL,
                      model.file = model_file,
                      parallel = TRUE,
                      n.chains = 3,
                      n.iter = 10,
                      DIC = TRUE)



