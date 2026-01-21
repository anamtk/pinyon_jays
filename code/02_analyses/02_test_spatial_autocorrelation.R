#Checking for spatial autocorrelation in residuals
#Ana Miller-ter Kuile and Kyle Rodman
#March 27, 2025

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "sf",  
                  "terra",
                  'readxl',
                  'sf',
                  'jagsUI',
                  'rjags',
                  'mcmcplots',
                  'coda',
                  # "DHARMa',
                  # 'BayesPostEst",
                  'ncf')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load data ---------------------------------------------------------------

ebird <- readRDS(here('data',
                      '01_ebird_data',
                      'cleaned_data',
                      '04_JAGS_indexIDs',
                      'ebird_check_blob_yr_ids.RDS')) %>%
  st_as_sf()

residuals_sum <- readRDS(here('monsoon',
                       'outputs',
                       'ebird_abund_model_residuals.RDS'))

# yrep_samps <- readRDS(here('monsoon',
#                          'ebird',
#                          'nospuncert',
#                          'outputs',
#                          'ebird_abund_model2_yrepsamples.RDS'))


# Get dataframes in order -------------------------------------------------

#ebird dataframe just needs numID, yrID, checkID, geometry

#residuals_sum needs to just be mean residuals and then prepped to 
#be tidy with those columns of numID, yrID, checkID, residuals

#yrep samples needs to be a matrix of observations x simulations

#prep ebird dataframe
ebird2 <- ebird %>%
  dplyr::select(numID, yrID, checkID, geometry, obsrvtn_c)

#prep residuals dataframe
resids <- as.data.frame(residuals_sum$statistics) %>%
  dplyr::select(Mean) %>%
  rownames_to_column(var = 'par') %>%
  filter(par != "deviance") %>%
  separate(par, 
           into = c('yrID', 'numID', 'checkID'),
           sep = ",") %>%
  mutate(yrID = str_sub(yrID, 7, nchar(yrID)),
         checkID = str_sub(checkID, 1, (nchar(checkID)-1))) %>%
  mutate(yrID = as.numeric(yrID),
         numID = as.numeric(numID),
         checkID = as.numeric(checkID)) %>%
  rename(residuals = Mean)

#yrep <- yrep_samps

# Combine dataframes ------------------------------------------------------


resids2 <- resids %>%
  left_join(ebird2, by = c("yrID", "numID", "checkID"))

par(mfrow = c(1,1))

spatial_fun <- function(yearID){
  
  df <- resids2 %>%
    filter(yrID == {{yearID}})
  
  coords <- as.data.frame(st_coordinates(df$geometry))
  
  plotspline <- plot(spline.correlog(coords$X, coords$Y,
                            df$residuals,
                            latlon = T,
                            resamp = 200,
                            xmax = 500),
                     main = yearID)
  
  return(plotspline)

}

years <- 1:13
lapply(years, spatial_fun)

## Now looking at residuals
#create DHARMa requires:
#simulatedResponse (yrep) (as a matrix of observations x simulations)
#observedResponse (y)
#optional:
#fittedPredictedResponse: this would be mu of a model (or lambda)
#but I'm not really sure how to convert this in this case because
#these are at the "blob" because of the detection component

# sim <- createDHARMa(simulatedResponse = yrep, 
#                     observedResponse = resids2$obsrvtn_c,
#                     integerResponse = F)
# 
# coords <- as.data.frame(st_coordinates(resids2$geometry))
# 
# # Spline correlogram for all points
# s <- spline.correlog(coords$X, coords$Y, resids2$residuals, 
#                      latlon = T, resamp = 200) ## 
# testSpatialAutocorrelation(refit, refit$X, refit$Y) ## Significant, but not particularly problematic autocorrelation here, especially given above test
# 

