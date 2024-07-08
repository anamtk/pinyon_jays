model{
  
  #loosely based off the joint-likelihood approach in this paper:
  #https://www.nature.com/articles/s41598-022-23603-0
  #and could use this paper as a good reference for ebird cleaning
  #-------------------------------------## 
  # Biological Process Model ###
  #-------------------------------------##
  
  for(i in 1:n.grids){
    for(t in 1:n.years){
    
    #latent "true" abundance is N, with rate parameter lambda
    N[i,t] ~ dpois(lambda[i,t])
    
    #lambda is based on a likelihood with SAM components
    log(lambda[i,t]) <- a0 + #could make this have random effects
      a[1]*AntCone[i,t] +
      a[2]*AntTemp[i,t] +
      a[3]*AntPPT[i,t] #would PPT be highly correlated with cones??
    
    #-------------------------------------## 
    # SAM summing ###
    #-------------------------------------##
    AntCone[i,t] <- sum(ConeTemp[i,t,])
    AntTemp[i,t] <- sum(TempTemp[i,t,])
    AntPPT[i,t] <- sum(PPTTemp[i,t,]) #would PPT be highly correlated with cones??
    
    #weight the lags for each covariate based on
    #the lag weight for each lag
    for(l in 1:n.lag){
      ConeTemp[i,t,l] <- Cone[i,t,l]*wA[l]
      TempTemp[i,t,l] <- Temp[i,t,l]*wB[l]
      PPTTemp[i,t,l] <- PPT[i,t,l]*wC[l]
      
      #any missing data can be imputed
      Cone[i,t,l] ~ dnorm(mu.cone, tau.cone)
      Temp[i,t,l] ~ dnorm(mu.temp, tau.temp)
      PPT[i,t,l] ~ dnorm(mu.ppt, tau.ppt)
    }
    
    #would it be better to do a "time since masting year" covariate? dunno
    
    } #years
  } #grids
  
  #-------------------------------------## 
  # Observation Models ###
  #-------------------------------------##
  
  #-------------------------------------## 
  # BBS Observation Model ###
  #-------------------------------------##
  
  #to deal with the fact taht not all bbs transectts 
  #are surveyed in all years... i think this should work
  #also bbs didn't run in 2020 so...
  for(t in 1:n.bbs.years){
    for(i in 1:n.bbs.trans[t]){
      
      #Detection probability is dependent on covariates
      logit(p.bbs[t,i]) <- b0 + #could make this have random effects
        #observer experience - make this a continous variable 
        #based on observer ID from 1966-onward
        b[1]*ObserverExp[t,i] #+
      #b[2]*AreaSampled[t,i] #look at OG paper and figure out how much of grid was 
      #sampled by each bbs survey - maybe don't include this, might always be the same...
      
      for(r in 1:n.bbs.points[t,i]){
        
        #bbs raw data:
        bbs.count[t,i,r] ~ dbin(p.bbs[t,i], N[bbs.grid[t,i], bbs.year[t]])
        
        #GOODNESS-OF-FIT EVALUATION
        bbs.count.rep[t,i,r] ~ dbin(p.bbs[t,i], N[bbs.grid[t, i], bbs.year[t]])
        
        
      }
    }
  }
  
  #If all years are in dataset for all transects (which it does'nt seem like is true),
  #then this loop should work:
  
  # for(i in 1:n.bbs.trans){
  #   for(t in 1:n.years){ #could make this the number of years for bbs if it's different from total
  #     
  #     #Detection probability is dependent on covariates
  #     logit(p.bbs[i,t]) <- b0 + #could make this have random effects
  #       #observer experience - make this a continous variable 
  #       #based on observer ID from 1966-onward
  #       b[1]*ObserverExp[i,t] #+
  #     #b[2]*AreaSampled[i,t] #look at OG paper and figure out how much of grid was 
  #     #sampled by each bbs survey - maybe don't include this, might always be the same...
  #     
  #     for(r in 1:n.bbs.points[i,t]){
  #   
  #       #bbs raw data:
  #       bbs.count[i,t,r] ~ dbin(p.bbs[i,t], N[bbs.grid[i], t])
  # 
  #       #GOODNESS-OF-FIT EVALUATION
  #       bbs.count.rep[i,t,r] ~ dbin(p.bbs[i,t], N[bbs.grid[i], t])
  #       
  #       } #bbs rep points
  #     } #bbs years
  #   } #bbs trans
  
  #-------------------------------------## 
  # eBIRD Observation Model ###
  #-------------------------------------##
  
  for(t in 1:n.years){ #could make this different for ebird than other loops if the data has more time
    for(i in 1:n.ebird.grids[t]){ #number of grid cells with ebird data i year t
      for(r in 1:n.ebird.check[t,i]){
        
        #ebird raw data
        ebird.count[t,i,r] ~ dbin(p.ebird[t,i,r], N[ebird.grid[t,i], t])
        
        #Detection probability is dependent on covariates
        logit(p.ebird[t,i,r]) <- c0 + #could make this have random effects, maybe observer ID
          #c1[SurveyType[t,i,r]] +
          c[2]*StartTime[t,i,r] +
          c[3]*Duration[t,i,r] +
          c[4]*Distance[t,i,r] +
          c[5]*NumObservers[t,i,r] 
        
        #checklists from June and July only to match up with BBS
        #only complete checklists
        #<5 h duration, <5km distance traveled, <10observers
        #may want to randomly sample one detection and one non-detection checklist
        # of a 5km hexagonal grid 
        #if imbalanced detection-non-detection per grid, subsample
        #only one-non-detection checklist for each cell in each week. 
        #start time, duration, distance, survey type 
        #(only considered stationary and the other one), number of observers,
        #observer ID number 
        #in that paper, they also divided time in a weird way to capture morning
        #and night - but look at htat and decide whether it's worth it, mayb3
        #doesnt matter when only thinking about one species
        
        #GOODNESS-OF-FIT EVALUATION
        ebird.count.rep[t,i,r] ~ dbin(p.ebird[t,i,r], N[ebird.grid[t,i], t])
        
      } 
    } 
  } 
  
  #if all grids have observations in all years, which is unlikely,
  #can do this loop
  
  # for(i in 1:n.ebird.grids){ #number of grid cells???
  #   for(t in 1:n.years){ #could make this different for ebird than other loops if the data has more time
  #     for(r in 1:n.ebird.check[i,t]){
  #   
  #   #ebird raw data
  #   ebird.count[i,t,r] ~ dbin(p.ebird[i,t,r], N[ebird.grid[i], t])
  #   
  #   #Detection probability is dependent on covariates
  #   logit(p.ebird[i,t,r]) <- c0 + #could make this have random effects, maybe observer ID
  #     #c1[SurveyType[i,t,r]] +
  #     c[2]*StartTime[i,t,r] +
  #     c[3]*Duration[i,t,r] +
  #     c[4]*Distance[i,t,r] +
  #     c[5]*NumObservers[i,t,r] 
  #   
  #     #checklists from June and July only to match up with BBS
  #     #only complete checklists
  #     #<5 h duration, <5km distance traveled, <10observers
  #     #may want to randomly sample one detection and one non-detection checklist
  #     # of a 5km hexagonal grid 
  #     #if imbalanced detection-non-detection per grid, subsample
  #     #only one-non-detection checklist for each cell in each week. 
  #     #start time, duration, distance, survey type 
  #   #(only considered stationary and the other one), number of observers,
  #     #observer ID number 
  #     #in that paper, they also divided time in a weird way to capture morning
  #     #and night - but look at htat and decide whether it's worth it, mayb3
  #   #doesnt matter when only thinking about one species
  #   
  #   #GOODNESS-OF-FIT EVALUATION
  #   ebird.count.rep[i,t,r] ~ dbin(p.ebird[i,t,r], N[ebird.grid[i], t])
  #   
  #     } #ebird checklists
  #   } #ebird years
  # } #ebird sites
  #     
  #-------------------------------------## 
  # Priors ###
  #-------------------------------------##
  
  #-------------------------------------## 
  # Priors: Biological process model ###
  #-------------------------------------##
  
  #intercept and slope parameters
  a0 ~ dnorm(0, 1E-2)
  for(i in 1:3){
    a[i] ~ dnorm(0, 1E-2)
  }

  #Antecedent variable priors
  #Employing "delta trick" to give vector of weights dirichlet priors
  #this is doing the dirichlet in two steps 
  #see Ogle et al. 2015 SAM model paper in Ecology Letters
  sumA <- sum(deltaA[])
  sumB <- sum(deltaB[])
  sumC <- sum(deltaC[])
  
  for(l in 1:n.lag){
    #weights for cone lags
    wA[l] <- deltaA[l]/sumA
    #relatively uninformative gamma prior
    deltaA[l] ~ dgamma(1,1)
    #derived quantity of cumulative weight
    cum.cone.wt[l] <- sum(wA[1:l])
    
    #weights for temp lags
    wB[l] <- deltaB[l]/sumB
    #relatively uninformative gamma prior
    deltaB[l] ~ dgamma(1,1)    
    #derived quantity of cumulative weight
    cum.temp.wt[l] <- sum(wB[1:l])
    
    #weights for ppt lags
    wC[l] <- deltaC[l]/sumC
    #relatively uninformative gamma prior
    deltaC[l] ~ dgamma(1,1)
    #derived quantity of cumulative weight
    cum.ppt.wt[l] <- sum(wC[1:l])
    
  }
  
  
  #-------------------------------------## 
  # Priors: BBS Observation model ###
  #-------------------------------------##

  #intercept
  b0 ~ dnorm(0, 1E-2)
  #slopes
  for(i in 1:1){
    b[i] ~ dnorm(0, 1E-2)
  }
  
  #-------------------------------------## 
  # Priors: eBIRD Observation model ###
  #-------------------------------------##
  
  #intercept
  c0 ~ dnorm(0, 1E-2)
  
  #categorical covariate
  # for(i in 2:2){
  #   c1[i] ~ dnorm(0, 1E-2)
  # }
  
  #cell-referenced by setting baseline level to 0
  c1[1] <- 0
  
  #continuous covariate slope paramters
  for(i in 2:5){
    c[i] ~ dnorm(0, 1E-2)
  }
  
  #-------------------------------------## 
  # Priors: Missing data ###
  #-------------------------------------##
  
  mu.cone ~ dunif(-10,10)
  sig.cone ~ dunif(0, 20)
  tau.cone <- pow(sig.cone, -2)
  mu.temp ~ dunif(-10,10)
  sig.temp ~ dunif(0, 20)
  tau.temp <- pow(sig.temp, -2)
  mu.ppt ~ dunif(-10,10)
  sig.ppt ~ dunif(0, 20)
  tau.ppt <- pow(sig.ppt, -2)
  
  #-------------------------------------## 
  # Derived quantities ###
  #-------------------------------------##
  
  for(t in 1:n.years){
   N.tot[t] <- sum(N[,t]) 
  }

  
  
}