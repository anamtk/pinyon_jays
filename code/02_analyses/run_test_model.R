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
                          'bbs_ebird_joint_data_list_nospuncert.RDS'))


# Model path --------------------------------------------------------------


model_file <- (here('code',
               'jags_models',
               'ebird_bbs_joint_abund_JAGS_nospuncert.R'))


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
                           'bbs_ebird_joint_init_list_nospuncert.RDS'))

# Run model ---------------------------------------------------------------

start <- Sys.time()
#this takes 5 minutes hahahahah
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

