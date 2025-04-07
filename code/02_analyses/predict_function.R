#EVENTUALLY:
#get this to loop through posterior samples for covariate effects
for(t in 1:n.years){
  for(i in 1:n.blobs[t]){
    
    #THIS SECTION NEEDS:
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
    
    #transform lambda
    lambda[t,i] <- exp(log_lambda[t,i])
    
    #BLOB N CODE:
    #if blob N:
    N[t,i] <- rpois(lambda[t,i]*blobArea[t,i])
    
    for(r in 1:n.ebird.check[t,i]){
      
      logit_p[t,i,r] <- c0 +          
        c1[SurveyType[t,i,r]] +
        c[2]*StartTime[t,i,r] +
        c[3]*Duration[t,i,r] +
        c[4]*(Distance[t,i,r]/Duration[t,i,r]) +
        c[5]*NumObservers[t,i,r] 
      
      #get p out of logit scale
      p[t,i,r] <- plogis(logit_p[t,i,r])
      
      #CHECKLIST N CODE
      #if checklist N:
      N[t,i,r] <- rpois(lambda[t,i]*listArea[t,i,r])
      
      #BLOB N CODE:
      #for blob N, to get at the total in the 
      #checklist relative to total N (per area of blob/checklist)
      N_check[t,i,r] <- round(N[t,i]*(listArea[t,i,r]/blobArea[t,i]))

      #Get residuals
      
      #CHECKLIST N CODE:
      resid[t,i,r] <- y[t,i,r] - p[t,i,r]*N[t,i,r]
      
      rmse[t,i,r] <- ((p[t,i,r]*N[t,i,r] - y[t,i,r])^2)/n.data
      
      #BLOB N CODE:
      resid[t,i,r] <- y[t,i,r] - p[t,i,r]*N_check[t,i,r]
      
      rmse[t,i,r] <- ((p[t,i,r]*N_check[t,i,r] - y[t,i,r])^2)/n.data
        
    }
  }
}

RMSE <- sqrt(sum(rmse[]))

df <- bind_rows(y, p, N_check, N) %>%
  rowwise() %>%
  mutate(y_pred = p*N_check) #or N, depending on model

mod <- lm(y ~ y_pred,
          data = df)

summary(mod)
