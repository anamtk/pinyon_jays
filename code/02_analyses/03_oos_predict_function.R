
# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 'data.table',
                  'corrplot',
                  'sf', 'coda')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data objects -------------------------------------------------------

data <- readRDS(here('data',
                     '03_jags_input_data',
                     'oos',
                     'oos_ebird_data_list_nospuncert.RDS'))


betas <- readRDS(here('monsoon',
                      'outputs',
                      'ebird_abund_model2_summary.RDS'))

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

rmse_df <- bind_rows(as.data.frame(rmse_samples[[1]]),
                     as.data.frame(rmse_samples[[2]]),
                     as.data.frame(rmse_samples[[3]])) %>%
  dplyr::select(RMSE) %>%
  mutate(type = "test")

tidy_oos_df <- readRDS(here('data',
             '01_ebird_data',
             'cleaned_data',
             '05_oos',
             'oos_ebird_check_blob_yr_ids.RDS'))

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

predict_fun <- function(sample){
#posterior samples for:
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
AntCone  <- matrix(NA, nrow = n.years, ncol = max(n.blobs))
AntTmax <- matrix(NA, nrow = n.years, ncol = max(n.blobs))
AntPPT <- matrix(NA, nrow = n.years, ncol = max(n.blobs))
ConeTemp <- array(NA, dim = c(n.years, max(n.blobs), n.lag))
TmaxTemp<- array(NA, dim = c(n.years, max(n.blobs), n.clag))
PPTTemp<- array(NA, dim = c(n.years, max(n.blobs), n.clag))

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
    
    #-------------------------------------## 
    # SAM summing ###
    #-------------------------------------##
    #weight the lags for each covariate based on
    #the lag weight for each lag
    for(l in 1:n.lag){
      ConeTemp[t,i,l] <- Cone[t,i,l]*wA[l]
    }
    
    for(l in 1:n.clag){
      TmaxTemp[t,i,l] <- Temp[t,i,l]*wB[l]
      PPTTemp[t,i,l] <- PPT[t,i,l]*wC[l]
    }
    
    AntCone[t,i] <- sum(ConeTemp[t,i,], na.rm = T)
    AntTmax[t,i] <- sum(TmaxTemp[t,i,], na.rm = T)
    AntPPT[t,i] <- sum(PPTTemp[t,i,], na.rm = T)
    
    #THIS SECTION NEEDS: (Not sure what i was writing to myself here...)
    #A way to convert observed cone/climate lag data
    #to a weighted mean based on posteriors 
    #of weights and the OG data
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
      
      #Get residuals
      
      #CHECKLIST N CODE:
      resid[t,i,r] <- ebird.count[t,i,r] - p.ebird[t,i,r]*N[t,i,r]
      
      #to get rmse
      sqr[t,i,r] <- (ebird.count[t,i,r] - p.ebird[t,i,r]*N[t,i,r])^2
      
      #get replicated data
      ebird.count.rep[t,i,r] <- rbinom(1, N[t,i,r], p.ebird[t,i,r])
      #ebird.count.rep[t,i,r] <- p.ebird[t,i,r]*N[t,i,r]
 
    }
    
  }
  
}

#ebird count and count rep R2
ebird.count.df <- as.data.frame(ebird.count) %>%
  pivot_longer(everything(),
               names_to = "val",
               values_to= "count") 

ebird.count.rep.df <- as.data.frame(ebird.count.rep) %>%
  pivot_longer(everything(),
               names_to = "val",
               values_to= "count.rep") 

# m1 <- cor(ebird.count.rep.df$count.rep, ebird.count.df$count,
#           use = "complete.obs")

model <- lm(ebird.count.rep.df$count.rep ~ 
              ebird.count.df$count)

sum <- summary(model)

R2 <- sum$adj.r.squared
  
#R2 <- m1

RMSE <- sqrt(sum(sqr[], na.rm = T)/n.checklists)

return(data.frame(RMSE, R2))

}


# Run the function for all covariate samples ------------------------------

n.sample <- length(unique(beta_samps$sample))
sample_nums <- 1:n.sample

#get the list of RMSE for each sample
sample_list <- lapply(sample_nums, predict_fun)

#turn into a DF
oos_df <- as.data.frame(do.call(rbind, sample_list)) %>%
  mutate(type = "oos")

#export
saveRDS(oos_RMSE_df, here('data',
                           '04_cross_validation',
                           'oos_RMSE.RDS'))

# Plot both RMSE ----------------------------------------------------------

both_RMSE <- oos_RMSE_df %>%
  bind_rows(rmse_df)
  
theme_set(theme_bw())
ggplot() +
  geom_boxplot(data = both_RMSE, aes(x = type, y = RMSE))

#WHAT IS THIS?? I DUNNO
#maybe looking at observed v predicted for the 
#predicted values??
# df <- bind_rows(y, p, N) %>%
#   rowwise() %>%
#   mutate(y_pred = p*N) 
# 
# mod <- lm(y ~ y_pred,
#           data = df)
# 
# summary(mod)
