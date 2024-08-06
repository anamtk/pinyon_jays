
# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Generate data -----------------------------------------------------------


#latent abundance
#indexing
n.grids <- 5 #number of smapling grids
n.years <- 5 #number of dataset years
n.lag <- 5 #number of lags for each covariate

#Evnrionmental covariates
Cone <- array(data = rnorm((n.grids*n.years*n.lag),0,1),
              dim= c(n.grids, n.years, n.lag))
Temp <- array(data = rnorm(n.grids*n.years*n.lag,0,1),
              dim= c(n.grids, n.years, n.lag))
PPT <- array(data = rnorm(n.grids*n.years*n.lag,0,1),
             dim= c(n.grids, n.years, n.lag))

#BBS
#indexing
n.bbs.trans <- 10 #number of bbs transects
#n.years #number of bbs years, hopefully the same as n.years
#number of points on transect i in year t
n.bbs.points <- array(data = 8, dim = c(n.bbs.trans, n.years) )
#vector of the grid from above for each bbs transect
bbs.grid <- c(1,1,2,2,3,3,4,4,5,5)

#response data
#bbs.count[i,t,r] #transect, year, point
bbs.count <- array(data = rpois((n.bbs.trans*n.years*8), lambda = 10),
                   dim = c(n.bbs.trans, n.years, 8))

#bbs covariates
ObserverExp <- array(data = rnorm(n.bbs.trans*n.years),
                     dim = c(n.bbs.trans, n.years))


#eBIRD
n.ebird.grids <- 3 #number of grids with ebird data - could this be the same as the total?
#n.ebird.check[i,t]#number of checklists in that grid
n.ebird.check <- array(data = 2,
                       dim = c(n.ebird.grids, n.years))
#ebird.grid[i] #vector of the ID grid for each ebird grid
ebird.grid <- c(1,2,3)

#response data
#ebird.count[i,t,r] #grid, year, checklist
ebird.count <- array(data = rpois((n.ebird.grids*n.years*n.ebird.check), 10),
                     dim = c(n.ebird.grids, n.years, 2))

#ebird covariates
#SurveyType[i,t,r] #categorical
SurveyType <- array(data = rbinom(n.ebird.grids*n.years*n.ebird.check, size = 1, prob = .6),
                    dim = c(n.ebird.grids, n.years, 2))
#StartTime[i,t,r] #when survey started
StartTime <- array(data = rnorm(n.ebird.grids*n.years*n.ebird.check),
                   dim = c(n.ebird.grids, n.years, 2))
#Duration[i,t,r] 
Duration  <- array(data = rnorm(n.ebird.grids*n.years*n.ebird.check),
                   dim = c(n.ebird.grids, n.years, 2))
#Distance[i,t,r]
Distance  <- array(data = rnorm(n.ebird.grids*n.years*n.ebird.check),
                   dim = c(n.ebird.grids, n.years, 2))
#NumObservers[i,t,r]
NumObservers  <- array(data = rnorm(n.ebird.grids*n.years*n.ebird.check),
                       dim = c(n.ebird.grids, n.years, 2))


# Initial values ----------------------------------------------------------

N<- array(data = 25,
                 dim = c(n.grids, n.years))

# Compile lists -----------------------------------------------------------

data_list <- list(n.grids = n.grids, 
                  n.years = n.years,
                  n.lag = n.lag,
                  Cone = Cone,
                  Temp = Temp,
                  PPT = PPT,
                  n.bbs.trans = n.bbs.trans,
                  n.bbs.points = n.bbs.points,
                  bbs.grid = bbs.grid,
                  bbs.count = bbs.count,
                  ObserverExp = ObserverExp,
                  n.ebird.grids = n.ebird.grids,
                  n.ebird.check = n.ebird.check,
                  ebird.grid = ebird.grid,
                  ebird.count = ebird.count,
                  SurveyType = SurveyType,
                  StartTime = StartTime,
                  Duration = Duration,
                  Distance = Distance,
                  NumObservers = NumObservers
                  )

inits <- list(list(N = N),
              list(N= N),
              list(N= N))

# Save data list ----------------------------------------------------------

saveRDS(data_list, here('data',
                        'simulated_data',
                        'data_list.RDS'))

saveRDS(inits, here('data',
                    'simulated_data',
                    'inits_list.RDS'))

