#Prep data for model
#Ana Miller-ter Kuile
#July 3, 2024

#prep all datasets for the model


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 'data.table',
                  'corrplot')

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

ppt <- read.csv(here('data',
                     'spatial_data',
                     'cleaned_data',
                     'ppt_data_df.csv'))

temp <- read.csv(here('data',
                      'spatial_data',
                      'cleaned_data',
                      'temp_data_df.csv'))

vpd <- read.csv(here('data',
                     'spatial_data',
                     'cleaned_data',
                     'vpd_data_df.csv'))

mons <- read.csv(here('data',
                      'spatial_data',
                      'cleaned_data',
                      'monsoon_data_df.csv'))

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

#subsetting without the last year given lack of overlap and forward
#projectiosn for cones for all data. Insetad of the last lag having
#~16.5% imputing, it now has ~10% and all the others have all
#their data
ebird3 <- ebird2b %>%
  filter(year > 2009 & year < 2022)
#12874 total
bbs2 <- bbs %>% 
  filter(Year > 2009 & Year < 2022)
#760 total

#combine and plot quickly
ebd4 <- ebird3 %>%
  dplyr::select(latitude, longitude, year) %>%
  mutate(datasource = 'ebird') %>%
  rename(Latitude = latitude,
         Longitude = longitude,
         Year = year)

bb3 <- bbs2 %>%
  dplyr::select(Longitude, Latitude, Year) %>%
  mutate(datasource = "bbs")

both <- ebd4 %>%
  rbind(bb3)

ggplot(both) +
  geom_point(aes(x = Longitude, y = Latitude, color = datasource)) +
  facet_wrap(~Year) +
  theme_bw() +
  scale_color_manual(values = c('#d8b365','#5ab4ac'))
  

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
#3870 - that have ebird and/or bbs AND cone data


# Data objects for jags ---------------------------------------------------


# Latent N loop -----------------------------------------------------------

#Total latent N:
n.grids <- nrow(all_cells)

#NOTE: Reducing years to not include 2022 of data so that we 
#can better represent the cone data - when i have the full
#dataset, ~16.5% of the final lag are being imputed for cone data
#and i'm wondering if that's why estimates have been weird.
n.years <- length(2010:2021)

#Evnrionmental covariates
#cones:
cones2 <- cones %>%
  #get only cell IDs that overlap with bird data
  filter(cell %in% all_cells$cellID) %>%
  dplyr::select(cell, X2000:X2023) %>%
  pivot_longer(X2000:X2023,
               names_to = "year",
               values_to = "cones") %>%
  mutate(year = str_sub(year, 2, nchar(year))) %>%
  mutate(year = as.numeric(year)) %>%
  mutate(yearID = as.numeric(as.factor(year)) - 11) %>%
  left_join(all_cells, by = c("cell" = "cellID")) %>%
  mutate(cones_0 = scale(cones)) %>%
  group_by(cell) %>% 
  arrange(cell, year) %>%
  #this creates a column for every lag 3 years previous,
  #this year, and 3 years years into future
  do(data.frame(., setNames(data.table::shift(.$cones_0, 1:3),
                            c('cones_l1', 'cones_l2', 'cones_l3')))) %>%
  do(data.frame(., setNames(data.table::shift(.$cones_0, 1:3, type = "lead"),
                            c('cones_n1', 'cones_n2', 'cones_n3')))) %>%
  ungroup() %>%
  filter(yearID >= 1) %>%
  filter(year < 2022) %>%
  dplyr::select(yearID, numID, cones_l1:cones_l3, cones_0, 
                cones_n1:cones_n3) %>%
  pivot_longer(cones_l1:cones_n3,
               names_to = 'lag',
               values_to = "cones") %>%
  mutate(lagID = case_when(lag == 'cones_l3' ~ 1,
                           lag == 'cones_l2' ~ 2,
                           lag == 'cones_l1' ~ 3,
                           lag == "cones_0" ~ 4,
                           lag == 'cones_n1' ~ 5,
                           lag == 'cones_n2' ~ 6,
                           lag == 'cones_n3' ~ 7))

#number of lags to consider
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

#VPD
#these are now monthly - and we will need to break them up 
#seasonally 

#figure out seasons WRT jays:
## when are they nesting?
## when are they doing other things??

#from Wiggins (2005) report:
#breeding: late Feb-early May
#May-June: dependent young
#August-Feb: fall and early winter: wander to find resources
#what do they do in July??? 
#thus, seasons would be months:
#Season 1: 2-4: breeding
#Season 2: 5-6: feeding dependent young
#Season 3: 7: unknown, just chillin?
#Season 4: 8-1: winter foraging, potentially seeking food elsewhere

# vpd2 <- vpd %>%
#   #get only cell IDs that overlap with bird data
#   filter(cellID %in% all_cells$cellID) %>%
#   dplyr::select(cellID, PRISM_vpdmax_stable_4kmM3_200001_bil:PRISM_vpdmax_stable_4kmM3_202312_bil) %>%
#   pivot_longer(PRISM_vpdmax_stable_4kmM3_200001_bil:PRISM_vpdmax_stable_4kmM3_202312_bil,
#                names_to = "date",
#                values_to = "vpd_l1") %>%
#   mutate(year = str_sub(date, 27, (nchar(date)-6)),
#          month = str_sub(date, 31, (nchar(date)-4))) %>%
#   mutate(year = as.numeric(year),
#          month = as.numeric(month)) %>%
#   dplyr::select(-date) %>%
#   pivot_wider(names_from = month,
#               values_from = vpd_l1) %>%
#   group_by(cellID) %>%
#   arrange(cellID, year) %>%
#   #determine next year january value to compute vpd season 4 below
#   mutate(lead_1 = lead(`1`)) %>%
#   rowwise() %>%
#   #get each seasonal VPD
#   mutate(vpd_1 = max(`2`, `3`, `4`),
#          vpd_2 = max(`5`, `6`),
#          vpd_3 = `7`,
#          vpd_4 = max(`8`, `9`, `10`, `11`, `12`, lead_1, na.rm = T)) %>%
#   dplyr::select(cellID, year, vpd_1:vpd_4) %>%
#   pivot_longer(vpd_1:vpd_4,
#                names_to = "season",
#                values_to = "vpd_l1") %>%
#   mutate(season = str_sub(season, nchar(season))) %>%
#   mutate(season = as.numeric(season)) %>% 
#   mutate(vpd_l1 = scale(vpd_l1)) %>%
#   arrange(cellID, year, season) %>%
#   #this creates a column for every lag 12 seasons previous,
#   do(data.frame(., setNames(data.table::shift(.$vpd_l1, 1:12),
#                             c('vpd_l2', 'vpd_l3', 'vpd_l4',
#                               'vpd_l5', 'vpd_l6', 'vpd_l7',
#                               'vpd_l8', 'vpd_l9', 'vpd_l10',
#                               'vpd_l11', 'vpd_l12', 'vpd_l13')))) %>%
#   ungroup() %>%
#   mutate(yearID = as.numeric(as.factor(year))-10) %>%
#   filter(yearID >= 1) %>%
#   #"current" season is season 2, the season when data started
#   #being collected
#   filter(season == 2) %>%
#   dplyr::select(-season) %>%
#   left_join(all_cells, by = "cellID") %>%
#   dplyr::select(yearID, numID, vpd_l1:vpd_l13) %>%
#   pivot_longer(vpd_l1:vpd_l13,
#                names_to = 'lag',
#                values_to = "vpd") %>%
#   mutate(lagID = str_sub(lag, 6, nchar(lag))) %>%
#   mutate(lagID = as.numeric(lagID)) %>%
#   dplyr::select(yearID, numID, lagID, vpd) %>%
#   filter(yearID < 13)
# 
# #now, generate IDs for the for loop where
# # we will populate the array
# vpdgrid <- vpd2$numID
# vpdyear <- vpd2$yearID
# vpdlag <- vpd2$lagID
# 
# #make a blank array
# VPD <- array(data = NA, dim = c(n.grids, n.years, n.clag))
# 
# #fill taht array based on the values in those columns
# for(i in 1:dim(vpd2)[1]){ #dim[1] = n.rows
#   #using info from the dataframe on the year, grid,
#   #and lag ID of row i 
#   # populate that space in the array with the column in
#   # the dataframe that corresponds to the cone data
#   # for that yearxgridxlag combo
#   VPD[vpdgrid[i], vpdyear[i], vpdlag[i]] <- as.numeric(vpd2[i,4])
# }
# 
# sum(is.na(VPD))

#temperature
temp2 <- temp %>%
  #get only cell IDs that overlap with bird data
  filter(cellID %in% all_cells$cellID) %>%
  dplyr::select(cellID, PRISM_vpdmax_stable_4kmM3_200001_bil:PRISM_vpdmax_stable_4kmM3_202312_bil) %>%
  pivot_longer(PRISM_vpdmax_stable_4kmM3_200001_bil:PRISM_vpdmax_stable_4kmM3_202312_bil,
               names_to = "date",
               values_to = "temp_l1") %>%
  mutate(year = str_sub(date, 27, (nchar(date)-6)),
         month = str_sub(date, 31, (nchar(date)-4))) %>%
  mutate(year = as.numeric(year),
         month = as.numeric(month)) %>%
  dplyr::select(-date) %>%
  pivot_wider(names_from = month,
              values_from = temp_l1) %>%
  group_by(cellID) %>%
  arrange(cellID, year) %>%
  #determine next year january value to compute temp season 4 below
  mutate(lead_1 = lead(`1`)) %>%
  rowwise() %>%
  #get each seasonal VPD
  mutate(temp_1 = max(`2`, `3`, `4`),
         temp_2 = max(`5`, `6`),
         temp_3 = `7`,
         temp_4 = max(`8`, `9`, `10`, `11`, `12`, lead_1, na.rm = T)) %>%
  dplyr::select(cellID, year, temp_1:temp_4) %>%
  mutate(temp_1 = scale(temp_1),
         temp_2 = scale(temp_2),
         temp_3 = scale(temp_3),
         temp_4 = scale(temp_4)) %>%
  pivot_longer(temp_1:temp_4,
               names_to = "season",
               values_to = "temp_l1") %>%
  mutate(season = str_sub(season, nchar(season))) %>%
  mutate(season = as.numeric(season)) %>% 
  #mutate(temp_l1 = scale(temp_l1)) %>%
  arrange(cellID, year, season) %>%
  #this creates a column for every lag 12 seasons previous,
  do(data.frame(., setNames(data.table::shift(.$temp_l1, 1:12),
                            c('temp_l2', 'temp_l3', 'temp_l4',
                              'temp_l5', 'temp_l6', 'temp_l7',
                              'temp_l8', 'temp_l9', 'temp_l10',
                              'temp_l11', 'temp_l12', 'temp_l13')))) %>%
  ungroup() %>%
  mutate(yearID = as.numeric(as.factor(year))-10) %>%
  filter(yearID >= 1) %>%
  #"current" season is season 2, the season when data started
  #being collected
  filter(season == 2) %>%
  dplyr::select(-season) %>%
  left_join(all_cells, by = "cellID") %>%
  dplyr::select(yearID, numID, temp_l1:temp_l13) %>%
  pivot_longer(temp_l1:temp_l13,
               names_to = 'lag',
               values_to = "temp") %>%
  mutate(lagID = str_sub(lag, 6, nchar(lag))) %>%
  mutate(lagID = as.numeric(lagID)) %>%
  dplyr::select(yearID, numID, lagID, temp) %>%
  filter(yearID < 13)

#lag index
n.clag <- max(temp2$lagID)

#now, generate IDs for the for loop where
# we will populate the array
tempgrid <- temp2$numID
tempyear <- temp2$yearID
templag <- temp2$lagID

#make a blank array
Temp <- array(data = NA, dim = c(n.grids, n.years, n.clag))

#fill taht array based on the values in those columns
for(i in 1:dim(temp2)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, grid,
  #and lag ID of row i 
  # populate that space in the array with the column in
  # the dataframe that corresponds to the data
  # for that yearxgridxlag combo
  Temp[tempgrid[i], tempyear[i], templag[i]] <- as.numeric(temp2[i,4])
}

#PPT
ppt2 <- ppt %>%
  #get only cell IDs that overlap with bird data
  filter(cellID %in% all_cells$cellID) %>%
  dplyr::select(cellID, PRISM_vpdmax_stable_4kmM3_200001_bil:PRISM_vpdmax_stable_4kmM3_202312_bil) %>%
  pivot_longer(PRISM_vpdmax_stable_4kmM3_200001_bil:PRISM_vpdmax_stable_4kmM3_202312_bil,
               names_to = "date",
               values_to = "ppt_l1") %>%
  mutate(year = str_sub(date, 27, (nchar(date)-6)),
         month = str_sub(date, 31, (nchar(date)-4))) %>%
  mutate(year = as.numeric(year),
         month = as.numeric(month)) %>%
  dplyr::select(-date) %>%
  pivot_wider(names_from = month,
              values_from = ppt_l1) %>%
  group_by(cellID) %>%
  arrange(cellID, year) %>%
  #determine next year january value to compute ppt season 4 below
  mutate(lead_1 = lead(`1`)) %>%
  rowwise() %>%
  #get each seasonal VPD
  mutate(ppt_1 = sum(`2`, `3`, `4`),
         ppt_2 = sum(`5`, `6`),
         ppt_3 = `7`,
         ppt_4 = sum(`8`, `9`, `10`, `11`, `12`, lead_1, na.rm = T)) %>%
  dplyr::select(cellID, year, ppt_1:ppt_4) %>%
  mutate(ppt_1 = scale(ppt_1),
         ppt_2 = scale(ppt_2),
         ppt_3 = scale(ppt_3),
         ppt_4 = scale(ppt_4)) %>%
  pivot_longer(ppt_1:ppt_4,
               names_to = "season",
               values_to = "ppt_l1") %>%
  mutate(season = str_sub(season, nchar(season))) %>%
  mutate(season = as.numeric(season)) %>% 
  #mutate(ppt_l1 = scale(ppt_l1)) %>%
  arrange(cellID, year, season) %>%
  #this creates a column for every lag 12 seasons previous,
  do(data.frame(., setNames(data.table::shift(.$ppt_l1, 1:12),
                            c('ppt_l2', 'ppt_l3', 'ppt_l4',
                              'ppt_l5', 'ppt_l6', 'ppt_l7',
                              'ppt_l8', 'ppt_l9', 'ppt_l10',
                              'ppt_l11', 'ppt_l12', 'ppt_l13')))) %>%
  ungroup() %>%
  mutate(yearID = as.numeric(as.factor(year))-10) %>%
  filter(yearID >= 1) %>%
  #"current" season is season 2, the season when data started
  #being collected
  filter(season == 2) %>%
  dplyr::select(-season) %>%
  left_join(all_cells, by = "cellID") %>%
  dplyr::select(yearID, numID, ppt_l1:ppt_l13) %>%
  pivot_longer(ppt_l1:ppt_l13,
               names_to = 'lag',
               values_to = "ppt") %>%
  mutate(lagID = str_sub(lag, 6, nchar(lag))) %>%
  mutate(lagID = as.numeric(lagID)) %>%
  dplyr::select(yearID, numID, lagID, ppt) %>%
  filter(yearID < 13)

#now, generate IDs for the for loop where
# we will populate the array
pptgrid <- ppt2$numID
pptyear <- ppt2$yearID
pptlag <- ppt2$lagID

#make a blank array
PPT <- array(data = NA, dim = c(n.grids, n.years, n.clag))

#fill taht array based on the values in those columns
for(i in 1:dim(ppt2)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, grid,
  #and lag ID of row i 
  # populate that space in the array with the column in
  # the dataframe that corresponds to the data
  # for that yearxgridxlag combo
  PPT[pptgrid[i], pptyear[i], pptlag[i]] <- as.numeric(ppt2[i,4])
}




#monsoon
#it seems there is more than one value per "cellID" for monsoonality?
#I'm not really sure why?
Monsoon <- mons %>%
  group_by(cellID) %>%
  summarise(monsoon = mean(PRISM_ppt_30yr_normal_800mM4_07_bil, na.rm = T)) %>%
  ungroup() %>%
  left_join(all_cells, by = "cellID") %>%
  filter(!is.na(numID)) %>%
  mutate(monsoon = scale(monsoon)) %>%
  column_to_rownames(var = "numID") %>%
  dplyr::select(monsoon) %>%
  as_vector()


# BBS data prep -----------------------------------------------------------

#BBS:
# bbs.trans.id <- bbs2 %>%
#   distinct(StateNum, Route) %>%
#   mutate(TransectID = 1:n()) 
n.bbs.years <- length(c(2010:2019, 2021:2021))
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
n.bbs.points <- matrix(data = 10, nrow = n.bbs.years, ncol = max(n.bbs.trans))


#this is the "Stop" columns in the bbs2 dataset [t,i,r] #year,transect, stop
#just going to do the first 10 stops so that we can have a smaller
#spatial scale to match up to the climate and cone data
#scale (~5miles; the grid data are 4km)
bbs_count_df <- bbs2 %>%
  dplyr::select(StateNum, Route, Year, Stop1:Stop10) %>%
  pivot_longer(Stop1:Stop10,
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
bbs.count <- array(data = NA, dim = c(n.bbs.years, max(n.bbs.trans), 10))

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
years <- c(2010:2021)
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
                  n.clag = n.clag,
                  Cone = Cone,
                  Temp = Temp,
                  PPT = PPT,
                  #VPD = VPD,
                  Monsoon = Monsoon,
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


