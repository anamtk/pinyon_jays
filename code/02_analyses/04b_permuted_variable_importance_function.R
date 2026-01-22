
# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 'data.table',
                  'corrplot',
                  'sf', 'coda')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# What this needs to do ---------------------------------------------------

#TO UPDATE
#Get weighted AntCone, AntTmax, AntPPT into format in the function
#and then run with those shuffled instead of raw data
#Loop through posterior samples??? 

#track RMSE in a predict function

#read in data that is the originally fitted data and the mean/median 
##results for parameters from the model

#Shuffle each variable a total of 10-20 times and track RMSE for each
##reshuffling.

#calculate importance (Importance = baselineRMSE - shuffledRMSE) for 
##each of these runs

#get an average importance for each variable

#relativize that importance to sum to 1 (100%)

# Load data objects -------------------------------------------------------

data <- readRDS(here('data',
                     '03_jags_input_data',
                     'ebird_data_list_nospuncert.RDS'))

betas <- readRDS(here('monsoon',
                      'outputs',
                      'ebird_abund_model2_summary.RDS'))

beta_df <- as.data.frame(betas$statistics) %>%
  rownames_to_column(var = "parm")

beta_samples <- readRDS(here('monsoon',
                             'outputs',
                             'ebird_abund_model_covariate_effect_samples.RDS'))

beta_samps <- bind_rows(as.data.frame(beta_samples[[1]]),
               as.data.frame(beta_samples[[2]]),
               as.data.frame(beta_samples[[3]])) %>%
  mutate(sample = 1:n()) %>%
  pivot_longer(-sample,
               names_to = "parm",
               values_to = "value")

rmse_samples <- readRDS(here('monsoon',
                             'outputs',
                             'ebird_abund_model_RMSE_samples.RDS'))

baselineRMSE <- bind_rows(as.data.frame(rmse_samples[[1]]),
                     as.data.frame(rmse_samples[[2]]),
                     as.data.frame(rmse_samples[[3]])) %>%
  dplyr::select(RMSE) %>%
  summarise(mean = mean(RMSE)) %>%
  as_vector()

# Prep data for the loop --------------------------------------------------

#need:
#data objects for loop:

#indexing values:
#n.years (one value)
n.years <- data$n.years
#n.blobs[t] (vector length n.years)
n.blobs <- data$n.blobs
#n.ebird.check[t,i] (matrix rows of years, columns of blobs)
n.ebird.check <- data$n.ebird.check
#n.lag (cones)
n.lag <- data$n.lag
#n.clag (ppt and tmax)
n.clag <- data$n.clag

#other variables: 
#listArea[t,i,r] (year, blob, checklist)
listArea <- data$listArea
#ebird.count[t,i,r] (y data for the count for each checklist)
ebird.count <- data$ebird.count
#n.checklists - total number of checklists in the dataset
n.checklists <- data$n.checklists

#biological covariate data:
#Cone[t,i,l] (year, blob, lag)
Cone <- data$Cone
#Temp[t,i,l] (year, blob, lag)
Temp <- data$Temp
#PPT[t,i,l] (year, blob, lag)
PPT <- data$PPT
#Monsoon[t,i]
Monsoon <- data$Monsoon
#PinyonBA[t,i]
PinyonBA <- data$PinyonBA
#checklist covariates:
#SurveyType[t,i,r]
SurveyType <- data$SurveyType
#StartTime[t,i,r]
StartTime <- data$StartTime
#Duration[t,i,r]
Duration <- data$Duration
#Speed[t,i,r]
Speed <- data$Speed
#NumObservers[t,i,r] 
NumObservers <- data$NumObservers

#UPDATE THISSSSSSS
cov_list <- list(
  "Cone" = Cone,
  "Temp" = Temp,
  "PPT" = PPT,
  "Monsoon" = Monsoon,
  "PinyonBA" = PinyonBA
)

predict_fun <- function(sample, covariate){
  

# Weight SAM variables to get an AntX for each ----------------------------

  #create empty variables to fill in the for loop for this
  AntCone  <- matrix(NA, nrow = n.years, ncol = max(n.blobs))
  AntTmax <- matrix(NA, nrow = n.years, ncol = max(n.blobs))
  AntPPT <- matrix(NA, nrow = n.years, ncol = max(n.blobs))
  ConeTemp <- array(NA, dim = c(n.years, max(n.blobs), n.lag))
  TmaxTemp<- array(NA, dim = c(n.years, max(n.blobs), n.clag))
  PPTTemp<- array(NA, dim = c(n.years, max(n.blobs), n.clag))
  
  #get weights from model
  #wA
  wA <- beta_samps %>%
    filter(sample == {{sample}}) %>%
    filter(str_detect(parm, "wA")) %>%
    dplyr::select(value) %>%
    as_vector()
  #wB
  wB <- beta_samps %>%
    filter(sample == {{sample}}) %>%
    filter(str_detect(parm, "wB")) %>%
    dplyr::select(value) %>%
    as_vector()
  #wC
  wC <- beta_samps %>%
    filter(sample == {{sample}}) %>%
    filter(str_detect(parm, "wC")) %>%
    dplyr::select(value) %>%
    as_vector()
  
  for(t in 1:n.years){
    for(i in 1:max(n.blobs)){
      for(l in 1:n.lag){
        
        ConeTemp[t,i,l] <- Cone[t,i,l]*wA[l]
        
      }
      
      for(k in 1:n.clag){
        TmaxTemp[t,i,k] <- Temp[t,i,k]*wB[k]
        PPTTemp[t,i,k] <- PPT[t,i,k]*wC[k]
      }
      
      AntCone[t,i] <- sum(ConeTemp[t,i,], na.rm = T)
      AntTmax[t,i] <- sum(TmaxTemp[t,i,], na.rm = T)
      AntPPT[t,i] <- sum(PPTTemp[t,i,], na.rm = T)
    }
  }
  
  #shuffle the covariate of interest
  shuffled_cone <- AntCone[sample(nrow(AntCone)),]
  shuffled_temp <- AntTmax[sample(nrow(AntTmax)),]
  shuffled_ppt <- AntPPT[sample(nrow(AntPPT)),]
  shuffled_monsoon <- Monsoon[sample(nrow(Monsoon)), ]
  shuffled_pba <- PinyonBA[sample(nrow(PinyonBA)), ]
  
  #replace OG covariate with the shuffled one,
  #otherwise, just keep the OG covariate order
  AntCone <- if(covariate == "Cone") shuffled_cone else AntCone
  AntTmax <- if(covariate == "Temp") shuffled_temp else AntTmax
  AntPPT <- if(covariate == "PPT") shuffled_ppt else AntPPT
  Monsoon <- if(covariate == "Monsoon") shuffled_monsoon else cov_list[["Monsoon"]]
  PinyonBA <- if(covariate == "PinyonBA") shuffled_pba else cov_list[["PinyonBA"]]

#posterior means for:
#a0
a0 <- beta_samps %>%
  filter(sample == {{sample}}) %>%
  filter(str_detect(parm, "a0")) %>%
  dplyr::select(value) %>%
  as_vector()
#a
a <- beta_samps %>%
  filter(sample == {{sample}}) %>%
  filter(str_detect(parm, "a")) %>%
  filter(!str_detect(parm, "wA|a0|deviance")) %>%
  dplyr::select(value) %>%
  as_vector()

#c0
c0 <- beta_samps %>%
  filter(sample == {{sample}}) %>%
  filter(str_detect(parm, "c0")) %>%
  dplyr::select(value) %>%
  as_vector()

#c
c <- beta_samps %>%
  filter(sample == {{sample}}) %>%
  filter(str_detect(parm, "c")) %>%
  filter(!str_detect(parm, 'c0|c1|deviance')) %>%
  dplyr::select(value) %>%
  as_vector()

#create empty objects to fill
#blob loops:
log_lambda <- matrix(NA, nrow = n.years, ncol = max(n.blobs))
lambda <- matrix(NA, nrow = n.years, ncol = max(n.blobs))

#checklist loops:
logit_p<- array(NA, dim = c(n.years, max(n.blobs), max(n.ebird.check, na.rm = T)))
p.ebird<- array(NA, dim = c(n.years, max(n.blobs), max(n.ebird.check, na.rm = T))) 
ebird.count.rep <- array(NA, dim = c(n.years, max(n.blobs), max(n.ebird.check, na.rm = T))) 

N<- array(NA, dim = c(n.years, max(n.blobs), max(n.ebird.check, na.rm = T)))

#GOF
resid<- array(NA, dim = c(n.years, max(n.blobs), max(n.ebird.check, na.rm = T)))

sqr<- array(NA, dim = c(n.years, max(n.blobs), max(n.ebird.check, na.rm = T)))

#EVENTUALLY:
#get this to loop through posterior samples for covariate effects
for(t in 1:n.years){
  for(i in 1:n.blobs[t]){
    
    log_lambda[t,i] <- a0 +
      a[1]*AntCone[t,i] +
      a[2]*AntTmax[t,i] +
      a[3]*AntPPT[t,i] +
      a[4]*Monsoon[t,i] +
      a[5]*PinyonBA[t,i] +
      a[6]*AntCone[t,i]*AntTmax[t,i] + 
      a[7]*AntCone[t,i]*AntPPT[t,i] + 
      a[8]*AntCone[t,i]*Monsoon[t,i] + 
      a[9]*AntCone[t,i]*PinyonBA[t,i]
    
    lambda[t,i] <- exp(log_lambda[t,i])
    
    for(r in 1:n.ebird.check[t,i]){
      
      logit_p[t,i,r] <- c0 +          
        #c1[SurveyType[t,i,r]] +
        c[1]*StartTime[t,i,r] +
        c[2]*Duration[t,i,r] +
        c[3]*Speed[t,i,r] +
        c[4]*NumObservers[t,i,r] 
      
      #get p out of logit scale
      p.ebird[t,i,r] <- plogis(logit_p[t,i,r])
      
      #CHECKLIST N CODE
      N[t,i,r] <- rpois(1, lambda[t,i]*listArea[t,i,r])
      
      #to get rmse
      sqr[t,i,r] <- (ebird.count[t,i,r] - p.ebird[t,i,r]*N[t,i,r])^2
      

    }
    
  }
  
}

RMSE <- sqrt(sum(sqr[], na.rm = T)/n.checklists)

importance <- baselineRMSE - RMSE

importance_df <- as.data.frame(cbind(importance = importance)) %>%
  mutate(variable = {{covariate}})

return(importance_df)

}


# Run the function for 10-20 shuffles of each covariate -------------------

t <- predict_fun(sample = 1, covariate = "PinyonBA")

#1000 samples
sample_vec <- 1:1000
sample_vec <- 1:100

cov_vec <- c("Cone", "Temp", "PPT", "Monsoon", 
                "PinyonBA")

apply_matrix <- expand.grid(sampleid = sample_vec, covariateid = cov_vec)

start <- Sys.time()
importance_df_all <- mapply(predict_fun, apply_matrix$sampleid, apply_matrix$covariateid)
end <- Sys.time()
(end-start)

result_df_purrr <- map2_dfr(apply_matrix$sampleid, apply_matrix$covariateid, predict_fun)

importance_df <- result_df_purrr

str(importance_df)

importance_sum <- importance_df %>%
  group_by(variable) %>%
  summarise(mean_importance = mean(importance)) %>%
  ungroup()

total_importance <- importance_sum %>%
  rowwise() %>%
  mutate(mean_importance = mean_importance*-1) %>%
  ungroup() %>%
  summarise(total = sum(mean_importance)) %>%
  dplyr::select(total) %>%
  as_vector()

importance_sum2 <- importance_sum %>%
  rowwise() %>%
  mutate(relative_importance = (mean_importance*-1)/total_importance)

saveRDS(importance_sum2, 
        here('data',
             '06_variable_importance',
             'variable_relative_importance.RDS'))
