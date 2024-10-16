#Getting initials for next model run
#Ana Miller-ter Kuile

#this script pulls out the last values from the chain with the lowest deviance
#to use as initials in a follow-up model run

# Load packages ---------------------------------------------------------------
(start.time <- Sys.time())


# Load packages,
package.list <- c("jagsUI", "coda",
                  'dplyr', 'stringr',
                  'magrittr', 'tidyr',
                  'mcmcplots','ggplot2',
                  'tibble') 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load model --------------------------------------------------------------

mod <- readRDS(file ="/scratch/atm234/pinyon_jays/outputs/ebird_bbs_joint_abund_model.RDS")


# Load initials for n -----------------------------------------------------

n_inits <- readRDS('/scratch/atm234/pinyon_jays/inputs/bbs_ebird_joint_init_list.RDS')

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


# Set initials ------------------------------------------------------------

#root nodes:
a0 <- as.vector(samp_df2$a0)
a <- as.vector(samp_df2$a)
b0 <- as.vector(samp_df2$b0)
b <- as.vector(samp_df2$b)
c0 <- as.vector(samp_df2$c0)
c <- as.vector(samp_df2$c)
N <- n_inits[[1]]$N

inits <- list(list(a0 = a0,
                   a = a,
                   b0 = b0,
                   b = b,
                   c0 = c0,
                   c = c,
                   N = N),
              list(a0 = a0 + 0.1,
                   a = a + 0.1,
                   b0 = b0 + 0.1,
                   b = b + 0.01,
                   c0 = c0 + 0.1,
                   c = c + 0.01,
                   N = N),
              list(a0 = a0 - 0.1,
                   a = a - 0.1,
                   b0 = b0 - 0.1,
                   b = b - 0.01,
                   c0 = c0 - 0.1,
                   c = c - 0.01,
                   N = N))

# Load data ---------------------------------------------------------------

data <- readRDS("/scratch/atm234/pinyon_jays/inputs/bbs_ebird_joint_data_list.RDS")

# Define data list --------------------------------------------------------

data_list <- list(n.grids = data$n.grids,
                  n.years = data$n.years,
                  n.lag = data$n.lag,
                  Cone = data$Cone,
                  Temp = data$Temp,
                  PPT = data$PPT,
                  VPD = data$VPD,
                  #BBS loop
                  n.bbs.years = data$n.bbs.years,
                  n.bbs.trans = data$n.bbs.trans,
                  ObserverExp = data$ObserverExp,
                  n.bbs.points = data$n.bbs.points,
                  bbs.count = data$bbs.count,
                  bbs.grid = data$bbs.grid,
                  bbs.year = data$bbs.year,
                  #ebird loop
                  n.ebird.grids = data$n.ebird.grids,
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

# Run model ---------------------------------------------------------------

model <- jagsUI::jags(data = data_list,
                      parameters.to.save = parameters,
                      inits = inits,
                      model.file = '/scratch/atm234/pinyon_jays/inputs/ebird_bbs_joint_abund_JAGS.R',
                      parallel = TRUE,
                      n.chains = 3,
                      n.iter = 9000,
                      n.burnin = 5000,
                      DIC = TRUE)

#save as an R data object
saveRDS(model, 
        file ="/scratch/atm234/pinyon_jays/outputs/ebird_bbs_joint_abund_model2.RDS")

(end.time <- Sys.time())

(tot.time <- end.time - start.time)

# Diagnose model ----------------------------------------------------------

mcmcplot(model$samples,
         dir = "/scratch/atm234/pinyon_jays/outputs/mcmcplots/test2")