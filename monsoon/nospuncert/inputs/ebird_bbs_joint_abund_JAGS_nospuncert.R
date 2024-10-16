model{
  
  #Notes about this model:
  #We are combining data from BBS and eBIRD with a joint likelihood
  
  #loosely based off the joint-likelihood approach in this paper:
  #https://www.nature.com/articles/s41598-022-23603-0
  
  #We cleaned ebird data following that paper and guidelines from
  #the Cornell Lab Website
  #https://cornelllabofornithology.github.io/ebird-best-practices/abundance.html
  
  #We selected BBS data from only the first 10 points so that the 
  #spatial scale would be more in line with the spatial scale of covariates
  #(~5 miles; covariates have a 4km resolution)
  
  #Components of model:
  # - Both eBIRD and BBS contribute equally to predictions of yearly latent population
  ## sizes for pinyon jays
  # - Both of these observation components of the model have covariates to detection
  # - We incorporate spatial uncertainty for both eBIRD and BBS observations using a 
  ## sampler that selects a "gridID" for each observation in each year based on
  ## the relative coverage of the set of possible 4km grid cells that that observation 
  ## lat-long plus a buffer around it (5 mile radius for BBS; radius of 1/2 the distance 
  ## variable for eBIRD) could be in. So, for example, if an observation buffer
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
  
  for(i in 1:n.grids){
    for(t in 1:n.years){
    
    #latent "true" abundance is N, with rate parameter lambda
    N[i,t] ~ dpois(lambda[i,t])
    
    #lambda is based on a likelihood with SAM components
    log(lambda[i,t]) <- a0 + #could make this have random effects
      a[1]*AntCone[i,t] +
      a[2]*AntTmax[i,t] +
      a[3]*AntPPT[i,t] +
      a[4]*Monsoon[i] +
      a[5]*PinyonBA[i,t] +
      a[6]*AntCone[i,t]*AntTmax[i,t] + #would these be three-way with pinyonBA??
      a[7]*AntCone[i,t]*AntPPT[i,t] + #would these be three-way with pinyonBA??
      a[8]*AntCone[i,t]*Monsoon[i] + #would these be three-way with pinyonBA??
      a[9]*AntCone[i,t]*PinyonBA[i,t]
    
    #-------------------------------------## 
    # SAM summing ###
    #-------------------------------------##
    AntCone[i,t] <- sum(ConeTemp[i,t,])
    AntTmax[i,t] <- sum(TmaxTemp[i,t,])
    AntPPT[i,t] <- sum(PPTTemp[i,t,]) #would PPT be highly correlated with cones??

    #weight the lags for each covariate based on
    #the lag weight for each lag
    for(l in 1:n.lag){
      ConeTemp[i,t,l] <- Cone[i,t,l]*wA[l]
      
      #any missing data can be imputed
      Cone[i,t,l] ~ dnorm(mu.cone, tau.cone)
    }
    
    for(l in 1:n.clag){
      TmaxTemp[i,t,l] <- Tmax[i,t,l]*wB[l]
      PPTTemp[i,t,l] <- PPT[i,t,l]*wC[l]
      
      #any missing data can be imputed
      Tmax[i,t,l] ~ dnorm(mu.temp, tau.temp)
      PPT[i,t,l] ~ dnorm(mu.ppt, tau.ppt)
    }
    
    #would it be better to do a "time since masting year" covariate? dunno
    #for any missing basal area data
    PinyonBA[i,t] ~ dnorm(mu.pba, tau.pba)
    
    } #years
    
    #for any missing monsoon data
    Monsoon[i] ~ dnorm(mu.mon, tau.mon)
    
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
      
      #Spatial uncertainty in BBS locations for populating
      #the covariates
      #pi and grid.array are "data" 
      #pi is the proportion of the likely survey in a given grid
      #cell, grid.array is the gridIDs that correspond to each of those
      #same probabilities and from which a sample grid ID is drawn
      #random index for the grid array - which background point to collect
      # bbs.grid.index[t,i] ~ dcat(bbs.pi[t,i,1:n.bbs.cells[t,i]])
      # #pulling ut the grid ID for that one
      # bbs.grid[t,i] <- bbs.grid.array[t,i, bbs.grid.index[t,i]]
    }
  }
  
  #-------------------------------------## 
  # eBIRD Observation Model ###
  #-------------------------------------##
  
  # for(t in 1:n.years){ #could make this different for ebird than other loops if the data has more time
  #   for(i in 1:n.ebird.pairs[t]){ #number of pairs/single ebird observations in year t
  #     for(r in 1:n.ebird.check[t,i]){#number of checklists in each pair - either 1 or 2
  #       
  #       #ebird raw data
  #       ebird.count[t,i,r] ~ dbin(p.ebird[t,i,r], N[ebird.grid[t,i], t])
  #       
  #       #Detection probability is dependent on covariates
  #       logit(p.ebird[t,i,r]) <- c0 + #could make this have random effects, maybe observer ID
  #         c1[SurveyType[t,i,r]] +
  #         c[2]*StartTime[t,i,r] +
  #         c[3]*Duration[t,i,r] +
  #         c[4]*Distance[t,i,r] +
  #         c[5]*NumObservers[t,i,r] 
  #       
  #       #checklists from June and July only to match up with BBS
  #       #only complete checklists
  #       #<5 h duration, <5km distance traveled, <10observers
  #       #may want to randomly sample one detection and one non-detection checklist
  #       # of a 5km hexagonal grid 
  #       #if imbalanced detection-non-detection per grid, subsample
  #       #only one-non-detection checklist for each cell in each week. 
  #       #start time, duration, distance, survey type 
  #       #(only considered stationary and the other one), number of observers,
  #       #observer ID number 
  #       #in that paper, they also divided time in a weird way to capture morning
  #       #and night - but look at htat and decide whether it's worth it, mayb3
  #       #doesnt matter when only thinking about one species
  #       
  #       #GOODNESS-OF-FIT EVALUATION
  #       # I think i need to code ebird grid by replicate as well
  #       ebird.count.rep[t,i,r] ~ dbin(p.ebird[t,i,r], N[ebird.grid[t,i], t])
  #     
  #     
  #     } 
  #     
  #     #Spatial uncertainty in eBIRD locations for populating
  #     #the covariates
  #     #pi and grid.array are "data" 
  #     #pi is the proportion of the likely survey in a given grid
  #     #cell, grid.array is the gridIDs that correspond to each of those
  #     #same probabilities and from which a sample grid ID is drawn
  #     #random index for the grid array - which background point to collect
  #     # ebird.grid.index[t,i] ~ dcat(ebird.pi[t,i,1:n.ebird.cells[t,i]])
  #     # #pulling ut the grid ID for that one
  #     # ebird.grid[t,i] <- ebird.grid.array[t,i,ebird.grid.index[t,i]]  
  #   } 
  #   
  # } 
  # 

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
  
  #-------------------------------------## 
  # Derived quantities ###
  #-------------------------------------##
  
  for(t in 1:n.years){
   N.tot[t] <- sum(N[,t]) 
  }

  
  
}