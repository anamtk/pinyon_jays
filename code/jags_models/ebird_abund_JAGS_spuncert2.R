model{
  
  #Notes about this model:
  #This is a model evaluating relationships between pinyon jay
  #abundance and cone abundance, mediated by other covariates 
  
  #we are using eBIRD data in this study because of its broad
  #spatial and temporal coverage. We are using cone production 
  #data from Andreas Wion (in review)
  
  #We cleaned ebird data following that paper and guidelines from
  #the Cornell Lab Website
  #https://cornelllabofornithology.github.io/ebird-best-practices/abundance.html
  
  #Components of model:
  # - eBIRD observed counts contribute to estimates of latent abundance
  # - the observation component of the model has covariates to detection
  # - We incorporate spatial uncertainty for eBIRD observations using a 
  ## sampler that selects a "gridID" for each observation in each year based on
  ## the relative coverage of the set of possible 4km grid cells that that observation 
  ## lat-long plus a buffer around it (radius of 1/2 the distance 
  ## variable) could be in. So, for example, if an observation buffer
  ## has a 20% coverage with cell 1 and an 80% coverage with cell 2, the sampler will
  ## have probabilities of selecting grid cell 1 and 2 of .2 and .8 - likely selecting
  ## grid cell 2 most of the time for which to populate covariates for that iteration
  ## of the model
  # - Latent abundance is modeled with environmental covariates
  # - The hypotheses of sam effects for those covariates with this component 
  ## and resulting covariate structures in the SAM model are as follows: 
  
  ## Cone production: pinyon jays may "predict" cone years following good climate years
  ## and/or years with production of other resources (cones come ~2 years after
  ## good climate years; oak and juniper have a ~1 year lag)
  ## pinyon jays may also have an immediate effect with cones if they predict this 
  ## production year or they may have lagged effects if good resource years spur
  ## reproduction
  ## Structure of cone data: cone data have both 'lag' and 'lead' effects, with the 
  ## middle value in the SAM being the current year, sandwiched by years before
  ## and after that year
  
  ## Temp/PPT: pinyon jays may have a lagged effect of climate. a good or bad year
  ## of climate can then trigger a population response
  ## Structure of data: this is more classic "SAM" structure. we have current
  ## year and a set of previous years for this variable

  #-------------------------------------## 
  # Biological Process Model ###
  #-------------------------------------##
  
  for(t in 1:n.years){
    for(i in 1:n.grids[t]){
      
      #latent "true" abundance is N, with rate parameter lambda
      N[t,i] ~ dpois(lambda[t,i])
      
      #lambda is based on a likelihood with SAM components
      log(lambda[t,i]) <- a0 + #could make this have random effects
        a[1]*AntCone[t,i] +
        a[2]*AntTmax[t,i] +
        a[3]*AntPPT[t,i] +
        a[4]*Monsoon[t,i] +
        a[5]*PinyonBA[t,i] +
        a[6]*AntCone[t,i]*AntTmax[t,i] + 
        a[7]*AntCone[t,i]*AntPPT[t,i] + 
        a[8]*AntCone[t,i]*Monsoon[t,i] + 
        a[9]*AntCone[t,i]*PinyonBA[t,i]
      
      #-------------------------------------## 
      # SAM summing ###
      #-------------------------------------##
      AntCone[t,i] <- sum(ConeTemp[t,i,])
      AntTmax[t,i] <- sum(TmaxTemp[t,i,])
      AntPPT[t,i] <- sum(PPTTemp[t,i,])
      
      #weight the lags for each covariate based on
      #the lag weight for each lag
      for(l in 1:n.lag){
        ConeTemp[t,i,l] <- Cone[t,i,l]*wA[l]
        
        #any missing data can be imputed
        Cone[t,i,l] ~ dnorm(mu.cone, tau.cone)
      }
      
      for(l in 1:n.clag){
        TmaxTemp[t,i,l] <- Tmax[t,i,l]*wB[l]
        PPTTemp[t,i,l] <- PPT[t,i,l]*wC[l]
        
        #any missing data can be imputed
        Tmax[t,i,l] ~ dnorm(mu.temp, tau.temp)
        PPT[t,i,l] ~ dnorm(mu.ppt, tau.ppt)
      }
      
      #would it be better to do a "time since masting year" covariate? dunno
      #for any missing basal area data
      PinyonBA[t,i] ~ dnorm(mu.pba, tau.pba)
      
      #for any missing monsoon data
      Monsoon[t,i] ~ dnorm(mu.mon, tau.mon)
      
    } #grids
    
    
    
  } #years
  
  #-------------------------------------## 
  # Observation Model ###
  #-------------------------------------##
  
  #-------------------------------------## 
  # eBIRD Observation Model ###
  #-------------------------------------##
  
  for(t in 1:n.years){ #could make this different for ebird than other loops if the data has more time
    for(i in 1:n.ebird.check[t]){ #number of checklists in year t

      #ebird raw data
      ebird.count[t,i] ~ dbin(p.ebird[t,i], N[t,ebird.grid[t,i]])
        
      #Detection probability is dependent on covariates
      logit(p.ebird[t,i]) <- c0 + #could make this have random effects, maybe observer ID
        c1[SurveyType[t,i]] +
        c[2]*StartTime[t,i] +
        c[3]*Duration[t,i] +
        c[4]*Distance[t,i] +
        c[5]*NumObservers[t,i] 
        
      #checklists from breeding season (late Feb - early May)
      #only complete checklists
      #<5 h duration, <5km distance traveled, <10observers
      #start time, duration, distance, survey type 
      #(only considered stationary and traveling), number of observers,
      #observer ID number 
      #in that paper, they also divided time in a weird way to capture morning
      #and night - but look at htat and decide whether it's worth it, mayb3
      #doesnt matter when only thinking about one species
        
      #GOODNESS-OF-FIT EVALUATION
      ebird.count.rep[t,i] ~ dbin(p.ebird[t,i], N[t,ebird.grid[t,i]])
 
      #Spatial uncertainty in eBIRD locations for populating
      #the covariates
      #pi and grid.array are "data" 
      #pi is the proportion of the likely survey in a given grid
      #cell, grid.array is the gridIDs that correspond to each of those
      #same probabilities and from which a sample grid ID is drawn
      #random index for the grid array - which background point to collect
      ebird.grid.index[t,i] ~ dcat(ebird.pi[t,i,1:n.ebird.cells[t,i]])
      #pulling ut the grid ID for that one
      ebird.grid[t,i] <- ebird.grid.array[t,i,ebird.grid.index[t,i]]  
    } 
    
  } 
  

  #-------------------------------------## 
  # Priors ###
  #-------------------------------------##
  
  #-------------------------------------## 
  # Priors: Biological process model ###
  #-------------------------------------##
  
  #intercept and slope parameters
  a0 ~ dnorm(0, 1E-2)
  
  for(i in 1:9){
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
    
  }
  
  for(l in 1:n.clag){
    
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
  # Priors: eBIRD Observation model ###
  #-------------------------------------##
  
  #intercept
  c0 ~ dnorm(0, 1E-2)
  
  #categorical covariate
  for(i in 2:2){
    c1[i] ~ dnorm(0, 1E-2)
  }
  
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
  mu.mon ~ dunif(-10,10)
  sig.mon ~ dunif(0, 20)
  tau.mon <- pow(sig.mon, -2)
  mu.temp ~ dunif(-10,10)
  sig.temp ~ dunif(0, 20)
  tau.temp <- pow(sig.temp, -2)
  mu.ppt ~ dunif(-10,10)
  sig.ppt ~ dunif(0, 20)
  tau.ppt <- pow(sig.ppt, -2)
  mu.pba ~ dunif(-10,10)
  sig.pba ~ dunif(0, 20)
  tau.pba <- pow(sig.pba, -2)
  
}