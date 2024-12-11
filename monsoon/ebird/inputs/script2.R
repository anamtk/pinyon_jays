# Load packages -----------------------------------------------------------

(start.time <- Sys.time())

package.list <- c("jagsUI", "coda",
                  'dplyr', 'stringr',
                  'magrittr', 'tidyr',
                  'mcmcplots','ggplot2',
                  'tibble') 

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load model --------------------------------------------------------------

mod <- readRDS("/scratch/atm234/pinyon_jays/ebird/outputs/ebird_abund_model.RDS")
      

# Get initials from previous model ----------------------------------------

#get the MCMC chains
samples <- mod$samples

#function to make each chain a dataframe
df_fun <- function(chain){
  df <- as.data.frame(chain) %>%
    rownames_to_column(var = "iteration")
  return(df)
}

#use that function on all list elements
samp_dfs <- lapply(samples, df_fun)

#make into one dataframe
samp_df <- bind_rows(samp_dfs, .id = "chain")

#get values for all parameters from the last iteration of the
#chain with the lowest deviance
samp_df2 <- samp_df %>%
  group_by(chain) %>%
  #get mean deviance by chain
  mutate(mean_dev = mean(deviance, na.rm = T)) %>%
  ungroup() %>%
  #get only the chain with the minimum average deviance
  filter(mean_dev == min(mean_dev)) %>%
  #pull out the final iteration from that chain
  filter(iteration == max(iteration)) %>%
  dplyr::select(-chain, -iteration,
                -deviance, -mean_dev) 

#root nodes
#a0
a0 <- as.vector(samp_df2$a0)
#a
a <- as.vector(samp_df2$a)

#cs
c0 <- as.vector(samp_df2$c0)
c1 <- as.vector(samp_df2$c1)
c <- as.vector(samp_df2$c)

# Load data ---------------------------------------------------------------

data <- readRDS("/scratch/atm234/pinyon_jays/ebird/inputs/ebird_data_list.RDS")

# Define data list --------------------------------------------------------

data_list <- list(#latent N loop:
                  n.grids = data$n.grids,
                  n.years = data$n.years,
                  n.lag = data$n.lag,
                  n.clag = data$n.clag,
                  Cone = data$Cone,
                  Temp = data$Temp,
                  PPT = data$PPT,
                  Monsoon = data$Monsoon,
                  PinyonBA = data$PinyonBA,
                  #ebird loop
                  n.ebird.grids = data$n.ebird.grids,
                  n.ebird.check = data$n.ebird.check,
                  ebird.count = data$ebird.count,
                  ebird.pi = data$ebird.pi,
                  ebird.grid.array = data$ebird.grid.array,
                  n.ebird.cells = data$n.ebird.cells,
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

inits_list <- readRDS('/scratch/atm234/pinyon_jays/ebird/inputs/ebird_init_list.RDS')

N <- inits_list[[1]]$N

inits_list <- list(list(N = N,
                        a0 = a0,
                        a = a,
                        c0 = c0,
                        c1 = c1,
                        c = c),
                   list(N = N,
                        a0 = a0 + 0.5,
                        a = a+ 0.25,
                        c0 = c0 +0.05 ,
                        c1 = c1 + 0.05,
                        c = c + 0.05),
                   list(N = N,
                        a0 = a0 - 0.5,
                        a = a- 0.25,
                        c0 = c0 -0.05 ,
                        c1 = c1 - 0.05,
                        c = c - 0.05))
# Run model ---------------------------------------------------------------

model <- jagsUI::jags(data = data_list,
                      parameters.to.save = parameters,
                      inits = inits_list,
                      model.file = '/scratch/atm234/pinyon_jays/ebird/inputs/ebird_abund_JAGS_spuncert2.R',
                      parallel = TRUE,
                      n.chains = 3,
                      n.iter = 15000,
                      n.thin = 5,
                      n.burnin = 5000,
                      DIC = TRUE)

#save as an R data object
saveRDS(model, 
        file ="/scratch/atm234/pinyon_jays/ebird/outputs/ebird_abund_model2.RDS")

(end.time <- Sys.time())

(tot.time <- end.time - start.time)

# Diagnose model ----------------------------------------------------------

mcmcplot(model$samples,
         dir = "/scratch/atm234/pinyon_jays/ebird/outputs/mcmcplots/run2")