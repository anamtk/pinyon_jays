#Prep data for model
#Ana Miller-ter Kuile
#July 3, 2024

#prep all datasets for the model


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 'data.table')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Steps -------------------------------------------------------------------

# Load datasets -----------------------------------------------------------

ebird <- read.csv(here('data',
                      'ebird_data',
                      'cleaned_data',
                      'all_ebird_data_cellIDs.csv'))

bbs <- read.csv(here('data',
                     'bbs_data',
                     'cleaned_data',
                     'pjay_data_co_nm_cellIDs.csv'))

#only use this if we use all cells and all years
# cellIDs <- read.csv(here('data',
#                          'spatial_data',
#                          'cleaned_data',
#                          'all_cellIDs.csv'))
cones <- read.csv(here('data',
                       'spatial_data',
                       'cleaned_data',
                       'cone_masting_df.csv'))

# Subset ebird data -------------------------------------------------------

#this paper: 
#https://www.nature.com/articles/s41598-022-23603-0
#subset one 1 and one zero dataset from each cell
#the chapter from cornell lab subsets "10 per site", not specifying whether
#they're 1 or 0 checklists
#https://cornelllabofornithology.github.io/ebird-best-practices/occupancy.html#occupancy-intro

#going the cornell lab route for now - can come back to revisit later
ebird2 <- ebird %>%
  group_by(year, cellID) %>%
  slice_sample(n = 10) %>%
  ungroup()

#the method from the paper
ebird_yes <- ebird %>%
  filter(observation_count > 0) %>%
  group_by(year, cellID) %>%
  slice_sample(n = 1) %>%
  ungroup()

ebird_no <- ebird %>%
  filter(observation_count == 0) %>%
  group_by(year, cellID) %>%
  slice_sample(n = 1) %>%
  ungroup()

ebird2b <- ebird_yes %>%
  rbind(ebird_no)

# Some things about data --------------------------------------------------

#start by subsetting 2010-onward, since these are good ebird years
#but could consider other subsets of years down the road

#regardless, ebird is likely to dominate estimates with the joint likelihood,
#i think. So might be worth considering how to subset it even more and/or
#do another way of data integration with these two datasets

# Subset years for both datasets ------------------------------------------

ebird3 <- ebird2b %>%
  filter(year > 2009 & year < 2023)

bbs2 <- bbs %>%
  filter(Year > 2009)


# Get cell IDs for subset -------------------------------------------------

#Get the unique cell IDs that have data - we will just
#model these cells since it is MUCH smaller than the total 
#dataset

ebird_cells <- ebird3 %>%
  distinct(cellID)

bbs_cells <- bbs2 %>%
  distinct(cellID)

#all the cells in the dataset:
all_cells <- ebird_cells %>%
  full_join(bbs_cells, by = "cellID") %>%
  #get a 1:n() ID for the cells instead of the random
  #ones we extracted, since this will be easier to index
  #in the model than random numbers taht aren't consecutive
  mutate(numID = 1:n())
#2733 - that have ebird and/or bbs AND cone data

# Data objects for jags ---------------------------------------------------


# Latent N loop -----------------------------------------------------------

#Total latent N:
n.grids <- nrow(all_cells)
n.years <- length(2010:2022)

#Evnrionmental covariates
cones2 <- cones %>%
  #get only cell IDs that overlap with bird data
  filter(cell %in% all_cells$cellID) %>%
  dplyr::select(cell, X2005:X2022) %>%
  pivot_longer(X2005:X2022,
               names_to = "year",
               values_to = "cones") %>%
  mutate(year = str_sub(year, 2, nchar(year))) %>%
  mutate(year = as.numeric(year)) %>%
  mutate(yearID = as.numeric(as.factor(year)) - 5) %>%
  left_join(all_cells, by = c("cell" = "cellID")) %>%
  mutate(cones_0 = scale(cones)) %>%
  group_by(cell) %>% 
  arrange(cell, year) %>%
  #this creates a column for every lag 1:4 years ago
  do(data.frame(., setNames(shift(.$cones_0, -6:6),
                            c('cones_l1', 'cones_l2', 'cones_l3', 'cones_l4',
                              "cones_l5", "cones_l6", "cones_l7", 'cones_l8',
                              'cones_l9', 'cones_l10', 'cones_l11', 'cones_l12',
                              'cones_l13')))) %>%
  ungroup() %>%
  filter(yearID >= 1) %>%
  dplyr::select(yearID, numID, cones_l1:cones_l13) %>%
  pivot_longer(cones_l1:cones_l13,
               names_to = 'lagID',
               values_to = "cones") %>%
  mutate(lagID = str_sub(lagID, 8, nchar(lagID))) %>%
  mutate(lagID = as.numeric(lagID))

#Covariates(Eventually)
n.lag <- max(cones2$lagID) #number of lags for each covariate

#now, generate IDs for the for loop where
# we will populate the array
conegrid <- cones2$numID
coneyear <- cones2$yearID
conelag <- cones2$lagID

#make a blank array
Cone <- array(data = NA, dim = c(n.grids, n.years, n.lag))

#fill taht array based on the values in those columns
for(i in 1:dim(cones2)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, grid,
  #and lag ID of row i 
  # populate that space in the array with the column in
  # the dataframe that corresponds to the cone data
  # for that yearxgridxlag combo
  Cone[conegrid[i], coneyear[i], conelag[i]] <- as.numeric(cones2[i,4])
}

# Cone[i,t,l]
# Temp[i,t,l]
# PPT[i,t,l]

# BBS data prep -----------------------------------------------------------

#BBS:
# bbs.trans.id <- bbs2 %>%
#   distinct(StateNum, Route) %>%
#   mutate(TransectID = 1:n()) 
n.bbs.years <- length(c(2010:2019, 2021:2022))
#this is dependent on survey year, since not all transects are 
#surveyed in every year - makes indexing in model easier, hopefully
bbs.trans.id.yr <- bbs2 %>%
  distinct(Year, StateNum, Route) %>%
  group_by(Year) %>%
  mutate(yrtransID = 1:n())
  
#number of transects surveyed in each year
n.bbs.trans <- bbs2 %>%
  distinct(Year, StateNum, Route) %>%
  group_by(Year) %>%
  tally() %>%
  ungroup() %>%
  dplyr::select(n) %>%
  as_vector()

#looks like this variable is skewed toward lower experience
#needs to be year x transect w/in year matrix
ObserverExp <- bbs2 %>%
  mutate(ObserverExp = scale(ObserverExp)) %>%
  left_join(bbs.trans.id.yr, by = c("StateNum", "Route", "Year")) %>%
  dplyr::select(Year, ObserverExp, yrtransID) %>%
  pivot_wider(names_from = yrtransID, 
              values_from = ObserverExp) %>%
  arrange(Year) %>%
  column_to_rownames(var = "Year") %>%
  as.matrix()

#each bbs transect has 50 points - not sure this part is going to work...
n.bbs.points <- matrix(data = 50, nrow = n.bbs.years, ncol = max(n.bbs.trans))


#this is the "Stop" columns in the bbs2 dataset [t,i,r] #year,transect, stop
bbs_count_df <- bbs2 %>%
  dplyr::select(StateNum, Route, Year, Stop1:Stop50) %>%
  pivot_longer(Stop1:Stop50,
               names_to = "Stop",
               values_to = "count") %>%
  #get the transect ID by year
  left_join(bbs.trans.id.yr, by = c("Year", "StateNum", "Route")) %>%
  #set year and stop IDs for looping through for count data below
  mutate(YearID = as.numeric(as.factor(Year)),
         StopID = str_sub(Stop, 5, nchar(Stop))) %>%
  mutate(StopID = as.numeric(StopID))


#now, generate IDs for the for loop where
# we will populate the matrix
yr <- bbs_count_df$YearID #get a yearID for each iteration of the loop
trans <- bbs_count_df$yrtransID #transect ID for each iteration fo the loop
stop <- bbs_count_df$StopID #get a stop ID for each iteration of the loop

#make a blank array
bbs.count <- array(data = NA, dim = c(n.bbs.years, max(n.bbs.trans), 50))

#fill taht array based on the values in those columns
for(i in 1:dim(bbs_count_df)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, transect,
  #and stop ID of row i 
  # populate that space in the array with the column in
  # the dataframe that corresponds to the count data
  # for that yearxtransectxstop combo
  bbs.count[yr[i], trans[i], stop[i]] <- as.numeric(bbs_count_df[i,5])
}

#get the main Grid ID for each BBS count for [t,i]
bbs.grid <- bbs2 %>%
  left_join(bbs.trans.id.yr, by = c("StateNum", "Route", "Year")) %>%
  distinct(StateNum, Route, cellID, Year, yrtransID) %>%
  left_join(all_cells, by = "cellID") %>%
  dplyr::select(Year, yrtransID, numID) %>%
  pivot_wider(names_from = yrtransID,
              values_from = numID) %>%
  arrange(Year) %>%
  column_to_rownames(var = 'Year') %>%
  as.matrix()

#get year infos to get year IDs for bbs 
years <- c(2010:2022)
yearids <- as.data.frame(years) %>%
  mutate(yearnum = 1:n())

#to account for the lack of survey in 2020. but to 
#link back to yearly numbers other than that year
bbs.year <- bbs2 %>%
  left_join(yearids, by = c("Year" = "years")) %>%
  arrange(Year) %>%
  distinct(yearnum) %>%
  as_vector()
  
# eBIRD loop --------------------------------------------------------------

#eBIRD:

#vector of # grids/year
n.ebird.grids <- ebird3 %>%
  distinct(year, cellID) %>%
  group_by(year) %>%
  tally() %>%
  ungroup() %>%
  dplyr::select(n) %>%
  as_vector()

ebird.trans.id.yr <- ebird3 %>%
  distinct(year, cellID) %>%
  group_by(year) %>%
  mutate(yrtransID = 1:n())
  
#number of "replicates" per grid in a given year [t,i]
n.ebird.check <- ebird3 %>%
  group_by(year, cellID) %>%
  tally() %>%
  ungroup() %>%
  left_join(ebird.trans.id.yr, by = c("year", "cellID")) %>%
  dplyr::select(-cellID) %>%
  arrange(yrtransID) %>%
  pivot_wider(names_from = yrtransID,
              values_from = n) %>%
  column_to_rownames(var = 'year') %>%
  as.matrix()

#[t,i] - this is the reference main grid ID 
#to link each ebird observation back to
ebird.grid <- ebird3 %>%
  distinct(year, cellID) %>%
  left_join(all_cells, by = "cellID") %>%
  left_join(ebird.trans.id.yr, by = c("cellID", "year")) %>%
  dplyr::select(-cellID) %>%
  pivot_wider(names_from = yrtransID,
              values_from = numID) %>%
  arrange(year) %>%
  column_to_rownames(var = "year") %>%
  as.matrix()

#create a dataframe that we can referecne for all the loops below
ebird_index_df <- ebird3 %>%
  left_join(all_cells, by = "cellID") %>%
  left_join(ebird.trans.id.yr, by = c("cellID", "year")) %>%
  mutate(yrID = as.numeric(as.factor(year))) %>%
  dplyr::select(observation_count, protocol_type,
                year, time_observations_started,
                duration_minutes, effort_distance_km,
                number_observers, cellID, numID, yrtransID, yrID) %>%
  group_by(cellID, year) %>%
  mutate(checkID = 1:n()) %>%
  ungroup() %>%
  mutate(SurveyType = case_when(protocol_type == "Stationary" ~ 1,
                                   protocol_type == "Traveling" ~ 2)) %>%
  #duration  
  #distance
  #NumObservers
  #time started (hours since midnight)
  mutate(StartTime = scale(time_observations_started),
         Duration = scale(duration_minutes),
         Distance = scale(effort_distance_km),
         NumObservers = scale(number_observers))
  #might need to get a random effect of observer, but let's wait on that 
  #for now and see how it goes with just what we have currently

#now, generate IDs for the for loop where
# we will populate the matrix
yr2 <- ebird_index_df$yrID #get a yearID for each iteration of the loop
grid <- ebird_index_df$yrtransID #transect ID for each iteration fo the loop
check <- ebird_index_df$checkID #get a stop ID for each iteration of the loop

#make a blank array
ebird.count <- array(data = NA, dim = c(n.years, max(n.ebird.grids), 2))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, transect,
  #and stop ID of row i 
  # populate that space in the array with the column in
  # the dataframe that corresponds to the count data
  # for that yearxtransectxstop combo
  ebird.count[yr2[i], grid[i], check[i]] <- as.numeric(ebird_index_df[i,1])
}

#do the same for the covariates below:

SurveyType <- array(data = NA, dim = c(n.years, max(n.ebird.grids), 2))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  SurveyType[yr2[i], grid[i], check[i]] <- as.numeric(ebird_index_df[i,13])
}

StartTime <- array(data = NA, dim = c(n.years, max(n.ebird.grids), 2))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  StartTime[yr2[i], grid[i], check[i]] <- as.numeric(ebird_index_df[i,14])
}

Duration <- array(data = NA, dim = c(n.years, max(n.ebird.grids), 2))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  Duration[yr2[i], grid[i], check[i]] <- as.numeric(ebird_index_df[i,15])
}

Distance <- array(data = NA, dim = c(n.years, max(n.ebird.grids), 2))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  Distance[yr2[i], grid[i], check[i]] <- as.numeric(ebird_index_df[i,16])
}

NumObservers <- array(data = NA, dim = c(n.years, max(n.ebird.grids), 2))
#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  NumObservers[yr2[i], grid[i], check[i]] <- as.numeric(ebird_index_df[i,17])
}

# Values for initials -----------------------------------------------------

#need a starting value for N[i,t]

maxb <- max(bbs.count, na.rm = T)
maxe <- max(ebird.count, na.rm=T)
nmax <- max(c(maxb, maxe))

N <- matrix(nmax, nrow =n.grids, ncol = n.years)

# Compile and export ------------------------------------------------------

data_list <- list(#latent N loop:
                  n.grids = n.grids,
                  n.years = n.years,
                  n.lag = n.lag,
                  Cone = Cone,
                  #Temp = Temp,
                  #PPT = PPT,
                  #BBS loop
                  n.bbs.years = n.bbs.years,
                  n.bbs.trans = n.bbs.trans,
                  ObserverExp = ObserverExp,
                  n.bbs.points = n.bbs.points,
                  bbs.count = bbs.count,
                  bbs.grid = bbs.grid,
                  bbs.year = bbs.year,
                  #ebird loop
                  n.ebird.grids = n.ebird.grids,
                  n.ebird.check = n.ebird.check,
                  ebird.count = ebird.count,
                  ebird.grid = ebird.grid,
                  SurveyType = SurveyType,
                  StartTime = StartTime,
                  Duration = Duration,
                  Distance = Distance,
                  NumObservers = NumObservers)

saveRDS(data_list, here('data',
                        'jags_input_data',
                        'bbs_ebird_joint_data_list.RDS'))

inits_list <- list(list(N = N),
                   list(N = N),
                   list(N = N))

saveRDS(inits_list, here('data',
                        'jags_input_data',
                        'bbs_ebird_joint_init_list.RDS'))


