# Load packages -----------------------------------------------------------

package.list <- c("tidyverse", 'here', #general packages
                  'jagsUI', #jags wrapper
                  'coda', #gelman.diag() function
                  'mcmcplots') #trace plot function

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

#theme_set(theme_classic())


# Load data ---------------------------------------------------------------

data_list <- readRDS(here('data',
                          'jags_input_data',
                          'test_ebird_data_list.RDS'))


# Model path --------------------------------------------------------------


model_file <- (here('code',
               'jags_models',
               'ebird_abund_JAGS_spuncert2.R'))


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

inits_list <- readRDS(here('data',
                           'jags_input_data',
                           'test_ebird_init_list.RDS'))

# Run model ---------------------------------------------------------------

start <- Sys.time()

model <- jagsUI::jags(data = data_list,
                      parameters.to.save = parameters,
                      inits = inits_list,
                      #inits = NULL,
                      model.file = model_file,
                      parallel = TRUE,
                      n.chains = 3,
                      n.iter = 1,
                      DIC = TRUE)

end <- Sys.time()
end-start

