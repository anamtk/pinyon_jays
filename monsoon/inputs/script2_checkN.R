# Load packages -----------------------------------------------------------

(start.time <- Sys.time())

package.list <- c("jagsUI", "coda",
                  'dplyr', 'stringr',
                  'magrittr', 'tidyr',
                  'tibble', 'purrr',
                  'mcmcplots','ggplot2') 

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load data ---------------------------------------------------------------

mod <- readRDS("/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model_checkN.RDS")

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
a <- samp_df2 %>%
  dplyr::select('a[1]':'a[9]') %>%
  pivot_longer('a[1]':'a[9]',
               names_to = 'a',
               values_to = "value") %>%
  dplyr::select(value) %>%
  as_vector()

#cs
c0 <- as.vector(samp_df2$c0)

#c1 <- samp_df2 %>%
#  dplyr::select("c1[2]") %>%
#  as_vector()

#c1 <- c(NA, c1)

c <- samp_df2 %>%
  dplyr::select('c[1]':'c[4]') %>%
  pivot_longer('c[1]':'c[4]',
               names_to = 'c',
               values_to = "value") %>%
  dplyr::select(value) %>%
  as_vector()

#c <- c(NA, c)

# Get data ----------------------------------------------------------------

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
  n.checklists = data$n.checklists)

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

N <- inits_list[[1]]$N

inits_list2 <- list(list(N = N,
                         a0 = a0,
                         a = a,
                         c0 = c0,
                         #c1 = c1,
                         c = c),
                    list(N = N,
                         a0 = a0 + 0.5,
                         a = a+ 0.25,
                         c0 = c0 +0.05 ,
                         #c1 = c1 + 0.05,
                         c = c + 0.05),
                    list(N = N,
                         a0 = a0 - 0.5,
                         a = a- 0.25,
                         c0 = c0 -0.05 ,
                         #c1 = c1 - 0.05,
                         c = c - 0.05))
# Run model ---------------------------------------------------------------

model <- jagsUI::jags(data = data_list,
                      parameters.to.save = parameters,
                      inits = inits_list2,
                      model.file = '/scratch/atm234/pinyon_jays/ebird/nospuncert/inputs/ebird_abund_JAGS_checkN_KO.R',
                      parallel = TRUE,
                      n.chains = 3,
                      n.iter = 25000,
                      n.burnin = 5000,
                      n.thin = 10,
                      DIC = TRUE)

#save as an R data object
saveRDS(model, 
        file ="/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model2_checkN.RDS")

(end.time <- Sys.time())

(tot.time <- end.time - start.time)

# Diagnose model ----------------------------------------------------------

mcmcplot(model$samples,
         dir = "/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/mcmcplots/checkN_run2")


