#Prep data for model
#Ana Miller-ter Kuile
#July 3, 2024

#prep all datasets for the model


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 'data.table',
                  'corrplot',
                  'sf')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Steps -------------------------------------------------------------------

# Load datasets -----------------------------------------------------------

crs_albers <- '+proj=longlat +datum=NAD83 +no_defs +type=crs'

ebird <- read_sf(here('data',
                       'ebird_data',
                       'cleaned_data',
                       'all_ebird_data_conefiltered.shp'))

ebird_gridIDs <- read.csv(here('data',
                               'ebird_data',
                               'cleaned_data',
                               'ebird_cellIDlists.csv'))

bbs <- read_sf(here('data',
                     'bbs_data',
                     'cleaned_data',
                     'all_bbs_data_conefiltered.shp'))

bbs_gridIDs <- read.csv(here('data',
                             'bbs_data',
                             'cleaned_data',
                             'bbs_cellIDlists.csv'))

#Covariates
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

# vpd <- read.csv(here('data',
#                      'spatial_data',
#                      'cleaned_data',
#                      'vpd_data_df.csv'))

mons <- read.csv(here('data',
                      'spatial_data',
                      'cleaned_data',
                      'monsoon_data_df.csv'))

ba <- read.csv(here('data',
                    'spatial_data',
                    'cleaned_data',
                    'pinyonba_data_df.csv'))

#all cells in cell lists to be able to get
#covariate values for all these cells
all_cells <- read.csv(here('data',
                         'spatial_data',
                         'cleaned_data',
                         'cellIDs.csv')) %>%
  rename(numID = X)

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
ebird3 <- ebird %>%
  filter(year > 2009 & year < 2022)

bbs2 <- bbs %>% 
  filter(Year > 2009 & Year < 2022)


#combine and plot quickly
ebd4 <- ebird3 %>%
  dplyr::select(geometry, year) %>%
  mutate(datasource = 'ebird') %>%
  rename(Year = year)

bb3 <- bbs2 %>%
  dplyr::select(geometry, Year) %>%
  mutate(datasource = "bbs")

both <- ebd4 %>%
  rbind(bb3)

ggplot(both) +
  geom_sf(aes(color = datasource)) +
  facet_wrap(~Year) +
  theme_bw() +
  scale_color_manual(values = c('#d8b365','#5ab4ac'))
  

# Get cell IDs for subset -------------------------------------------------

# #Get the unique cell IDs that have data - we will just
# #model these cells since it is MUCH smaller than the total 
# #dataset
# 
# ebird_cells <- ebird3 %>%
#   distinct(cellID)
# 
# bbs_cells <- bbs2 %>%
#   distinct(cellID)
# 
# #all the cells in the dataset:
# all_cells <- ebird_cells %>%
#   full_join(bbs_cells, by = "cellID") %>%
#   #get a 1:n() ID for the cells instead of the random
#   #ones we extracted, since this will be easier to index
#   #in the model than random numbers taht aren't consecutive
#   mutate(numID = 1:n())
# #3870 - that have ebird and/or bbs AND cone data
# 

# Data objects for jags ---------------------------------------------------


# Latent N loop -----------------------------------------------------------

#anywhere BA or Cones == 0, probs 0 because at end of 
#range - set to 0 instead of imputing in the model

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
  filter(cell %in% all_cells$cell) %>%
  dplyr::select(cell, X2000:X2023) %>%
  pivot_longer(X2000:X2023,
               names_to = "year",
               values_to = "cones") %>%
  mutate(year = str_sub(year, 2, nchar(year))) %>%
  mutate(year = as.numeric(year)) %>%
  mutate(yearID = as.numeric(as.factor(year)) - 11) %>%
  left_join(all_cells, by = c("cell")) %>%
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
  #only through 2021 right now
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

sum(is.na(Cone))/(sum(is.na(Cone)) + sum(!is.na(Cone)))

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
  filter(cellID %in% all_cells$cell) %>%
  dplyr::select(cellID, PRISM_tmax_stable_4kmM3_200001_bil:PRISM_tmax_stable_4kmM3_202312_bil) %>%
  pivot_longer(PRISM_tmax_stable_4kmM3_200001_bil:PRISM_tmax_stable_4kmM3_202312_bil,
               names_to = "date",
               values_to = "temp_l1") %>%
  mutate(year = str_sub(date, 25, (nchar(date)-6)),
         month = str_sub(date, 29, (nchar(date)-4))) %>%
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
  ungroup() %>%
  dplyr::select(cellID, year, temp_1:temp_4) %>%
  mutate(temp_1 = base::scale(temp_1),
         temp_2 = base::scale(temp_2),
         temp_3 = base::scale(temp_3),
         temp_4 = base::scale(temp_4)) %>%
  pivot_longer(temp_1:temp_4,
               names_to = "season",
               values_to = "temp_l1") %>%
  mutate(season = str_sub(season, nchar(season))) %>%
  mutate(season = as.numeric(season)) %>% 
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
  left_join(all_cells, by = c("cellID" = 'cell')) %>%
  dplyr::select(yearID, numID, temp_l1:temp_l13) %>%
  pivot_longer(temp_l1:temp_l13,
               names_to = 'lag',
               values_to = "temp") %>%
  mutate(lagID = str_sub(lag, 7, nchar(lag))) %>%
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

sum(is.na(Temp))/(sum(is.na(Temp)) + sum(!is.na(Temp)))
#PPT
ppt2 <- ppt %>%
  #get only cell IDs that overlap with bird data
  filter(cellID %in% all_cells$cell) %>%
  dplyr::select(cellID, PRISM_ppt_stable_4kmM3_200001_bil:PRISM_ppt_stable_4kmM3_202312_bil) %>%
  pivot_longer(PRISM_ppt_stable_4kmM3_200001_bil:PRISM_ppt_stable_4kmM3_202312_bil,
               names_to = "date",
               values_to = "ppt_l1") %>%
  mutate(year = str_sub(date, 24, (nchar(date)-6)),
         month = str_sub(date, 28, (nchar(date)-4))) %>%
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
  ungroup() %>%
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
  left_join(all_cells, by = c("cellID" = "cell")) %>%
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

sum(is.na(PPT))/(sum(is.na(PPT)) + sum(!is.na(PPT)))

#monsoon
Monsoon <- mons %>%
  group_by(cellID) %>%
  summarise(monsoon = mean(PRISM_ppt_30yr_normal_800mM4_07_bil, na.rm = T)) %>%
  ungroup() %>%
  left_join(all_cells, by = c("cellID" = 'cell')) %>%
  filter(!is.na(numID)) %>%
  mutate(monsoon = scale(monsoon)) %>%
  arrange(numID) %>%
  column_to_rownames(var = "numID") %>%
  dplyr::select(monsoon) %>%
  as_vector()

numID <- as.data.frame(1:7104) %>%
  rename(numID = `1:7104`)

#basal area
PinyonBA <- ba %>%
  left_join(all_cells, by = c("cellID" = "cell")) %>%
  filter(!is.na(numID)) %>%
  arrange(numID) %>%
  right_join(numID, by = "numID") %>%
  #filtering out 2022 for now
  dplyr::select(numID, PinyonBA_sqftPerAc_2010:PinyonBA_sqftPerAc_2021) %>%
  pivot_longer(PinyonBA_sqftPerAc_2010:PinyonBA_sqftPerAc_2021,
               names_to = 'year',
               values_to = 'ba') %>%
  mutate(ba = case_when(ba == Inf ~ NA,
                        TRUE ~ ba)) %>%
  mutate(ba = scale(ba)) %>%
  pivot_wider(names_from = year, 
              values_from = ba) %>%
  column_to_rownames(var = 'numID') %>%
  as.matrix()

sum(is.na(PinyonBA))/(sum(is.na(PinyonBA)) + sum(!is.na(PinyonBA)))

ba2 <- ba %>%
  left_join(all_cells, by = c("cellID" = "cell")) %>%
  filter(!is.na(numID)) %>%
  right_join(numID, by = "numID") %>%
  #filtering out 2022 for now
  dplyr::select(numID, PinyonBA_sqftPerAc_2010:PinyonBA_sqftPerAc_2021) %>%
  pivot_longer(PinyonBA_sqftPerAc_2010:PinyonBA_sqftPerAc_2021,
               names_to = 'year',
               values_to = 'ba') %>%
  mutate(ba = case_when(ba == Inf ~ NA,
                        TRUE ~ ba)) %>%
  group_by(year) %>%
  summarise(sum_na = sum(is.na(ba)),
            sum_notna = sum(!is.na(ba)),
            prop = sum_na/(sum_na + sum_notna))

# BBS data prep -----------------------------------------------------------

#BBS:
# bbs.trans.id <- bbs2 %>%
#   distinct(StateNum, Route) %>%
#   mutate(TransectID = 1:n()) 
n.bbs.years <- length(c(2010:2019, 2021:2021))
#this is dependent on survey year, since not all transects are 
#surveyed in every year - makes indexing in model easier, hopefully
bbs.trans.id.yr <- bbs2 %>%
  distinct(Year, StateNm, Route) %>%
  group_by(Year) %>%
  mutate(yrtransID = 1:n())
  
#number of transects surveyed in each year
n.bbs.trans <- bbs2 %>%
  distinct(Year, StateNm, Route) %>%
  group_by(Year) %>%
  tally() %>%
  ungroup() %>%
  dplyr::select(n) %>%
  as_vector()

#looks like this variable is skewed toward lower experience
#needs to be year x transect w/in year matrix
ObserverExp <- as.data.frame(bbs2) %>%
  mutate(ObserverExp = scale(ObsrvrE)) %>%
  left_join(bbs.trans.id.yr, by = c("StateNm", "Route", "Year")) %>%
  dplyr::select(Year, ObserverExp, yrtransID) %>%
  pivot_wider(names_from = 'yrtransID', 
              values_from = 'ObserverExp') %>%
  arrange(Year) %>%
  column_to_rownames(var = "Year") %>%
  as.matrix()

#each bbs transect has 50 points -
#filtering only the first 10 for spatial closeness
n.bbs.points <- matrix(data = 10, nrow = n.bbs.years, ncol = max(n.bbs.trans))


#this is the "Stop" columns in the bbs2 dataset [t,i,r] #year,transect, stop
#just going to do the first 10 stops so that we can have a smaller
#spatial scale to match up to the climate and cone data
#scale (~5miles; the grid data are 4km)
bbs_count_df <- as.data.frame(bbs2) %>%
  dplyr::select(StateNm, Route, Year, Stop1:Stop10) %>%
  pivot_longer(Stop1:Stop10,
               names_to = "Stop",
               values_to = "count") %>%
  #get the transect ID by year
  left_join(bbs.trans.id.yr, by = c("Year", "StateNm", "Route")) %>%
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

#NEED TO UPDATE NEXT SECTION USING
#GRID SAMPLER CODE AS A GUIDELINEEEE
#I think i will need to make the matrix
#to select form and for the proportions
#and make sure that any values past the nubmer of 
#grid cells from that year for pi are == 0

#need
bbs_grid_df <- bbs_gridIDs %>%
  filter(Year >= 2010 &  Year < 2022) %>%
  left_join(all_cells, by = "cell") %>%
  mutate(yearID = as.numeric(as.factor(Year))) %>%
  left_join(bbs.trans.id.yr, by = c("Year", "StateNm",
                                    "Route")) %>%
  group_by(yearID, yrtransID) %>%
  mutate(piID = 1:n()) %>%
  ungroup()
  

bgridyr <- bbs_grid_df$yearID
bgridtra <- bbs_grid_df$yrtransID
bgridid <- bbs_grid_df$piID

bbs.pi <- array(NA, dim = c(n.bbs.years, max(bgridtra), max(bgridid)))

#fill taht array based on the values in those columns
for(i in 1:dim(bbs_grid_df)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, transect,
  #and pi ID of row i
  # populate that space in the array with the column in
  # the dataframe that corresponds to the pi data
  # for that yearxtransectxpi combo
  bbs.pi[bgridyr[i], bgridtra[i], bgridid[i]] <- as.numeric(bbs_grid_df[i,7])
}

#set all values past the possible to 0 
#so the model doesn't break? i think i need this??
bbs.pi[is.na(bbs.pi)] <- 0

bbs.grid.array <- array(NA, dim = c(n.bbs.years, max(bgridtra), max(bgridid)))

#fill taht array based on the values in those columns
for(i in 1:dim(bbs_grid_df)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, transect,
  #and pi ID of row i
  # populate that space in the array with the column in
  # the dataframe that corresponds to the pi data
  # for that yearxtransectxpi combo
  bbs.grid.array[bgridyr[i], bgridtra[i], bgridid[i]] <- as.numeric(bbs_grid_df[i,8])
}


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

#will need to change the pi and indexing array to be 
#one less dimension

#eBIRD:

#vector of # pairs/year (sometimes singleton)
n.ebird.pairs <- ebird3 %>%
  distinct(year, pairID) %>%
  group_by(year) %>%
  tally() %>%
  ungroup() %>%
  dplyr::select(n) %>%
  as_vector()

#geting the IDs for each year-grid combo
#but keeping checklist ID here to blend down below
#more easily
# ebird.trans.id.yr <- ebird3 %>%
#   distinct(year,  pairID) %>%
#   group_by(year) %>%
#   #max for any year is 1831
#   mutate(yrtransID = 1:n()) %>%
#   ungroup()
  
#number of "replicates" per grid in a given year [t,i]
n.ebird.check <- as.data.frame(ebird3) %>%
  group_by(year, pairID) %>%
  tally() %>%
  ungroup() %>%
  left_join(ebird.trans.id.yr, by = c("year", "pairID")) %>%
  dplyr::select(-pairID) %>%
  arrange(yrtransID) %>%
  pivot_wider(names_from = yrtransID,
              values_from = n) %>%
  column_to_rownames(var = 'year') %>%
  as.matrix()

#create a dataframe that we can referecne for all the loops below
ebird_index_df <- as.data.frame(ebird3) %>%
  ungroup() %>%
  #left_join(all_cells, by = c("cellID" = "cell")) %>%
  #left_join(ebird.trans.id.yr, by = c("pairID", "year")) %>%
  mutate(yrID = as.numeric(as.factor(year))) %>%
  dplyr::select(chckls_, #checklist ID
                obsrvtn_c, #observation count
                prtcl_t,#protocol type
                year, 
                tm_bsr_, #time observation started
                drtn_mn, #duration_minutes
                effrt__, #effort distance km
                nmbr_bs, #number observers
               pairID, yrID, geometry) %>%
  group_by(pairID, yrID) %>%
  mutate(checkID = 1:n()) %>%
  ungroup() %>%
  mutate(SurveyType = case_when(prtcl_t == "Stationary" ~ 1,
                                prtcl_t == "Traveling" ~ 2)) %>%
  #duration  
  #distance
  #NumObservers
  #time started (hours since midnight)
  mutate(StartTime = scale(tm_bsr_),
         Duration = scale(drtn_mn),
         Distance = scale(effrt__),
         NumObservers = scale(nmbr_bs))
  #might need to get a random effect of observer, but let's wait on that 
  #for now and see how it goes with just what we have currently

#now, generate IDs for the for loop where
# we will populate the matrix
yr2 <- ebird_index_df$yrID #get a yearID for each iteration of the loop
pair <- ebird_index_df$pairID #pair ID for each iteration fo the loop
check <- ebird_index_df$checkID #get a checklist ID for each iteration of the loop

#make a blank array
ebird.count <- array(data = NA, dim = c(n.years, max(n.ebird.pairs), 2))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, transect,
  #and stop ID of row i 
  # populate that space in the array with the column in
  # the dataframe that corresponds to the count data
  # for that yearxtransectxstop combo
  ebird.count[yr2[i], pair[i], check[i]] <- as.numeric(ebird_index_df[i,2])
}

#do the same for the covariates below:

SurveyType <- array(data = NA, dim = c(n.years, max(n.ebird.pairs), 2))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  SurveyType[yr2[i], pair[i], check[i]] <- as.numeric(ebird_index_df[i,13])
}

StartTime <- array(data = NA, dim = c(n.years, max(n.ebird.pairs), 2))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  StartTime[yr2[i], pair[i], check[i]] <- as.numeric(ebird_index_df[i,14])
}

Duration <- array(data = NA, dim = c(n.years, max(n.ebird.pairs), 2))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  Duration[yr2[i], pair[i], check[i]] <- as.numeric(ebird_index_df[i,15])
}

Distance <- array(data = NA, dim = c(n.years, max(n.ebird.pairs), 2))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  Distance[yr2[i], pair[i], check[i]] <- as.numeric(ebird_index_df[i,16])
}

NumObservers <- array(data = NA, dim = c(n.years, max(n.ebird.pairs), 2))
#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  NumObservers[yr2[i], pair[i], check[i]] <- as.numeric(ebird_index_df[i,17])
}

ebird_index2 <- ebird_index_df %>%
  dplyr::select(chckls_, checkID, pairID)
  
#I need to link up metadata on the 
#"replicate" each checklist represents for each
#cell, and then the indexing will be [i,t,]
ebird_grid_df <- ebird_gridIDs %>%
  filter(year >= 2010 &  year < 2022) %>%
  #get cell numIDs for coding for jags
  left_join(all_cells, by = "cell") %>%
  mutate(yearID = as.numeric(as.factor(year))) %>%
  #get yrtransID index (which transect in year t are we )
  left_join(ebird_index2, by = c("chckls_", 'pairID')) %>%
  filter(!is.na(pairID)) %>%
  group_by(yearID, pairID, checkID) %>%
  mutate(piID = 1:n()) %>%
  ungroup() %>%
  dplyr::select(chckls_, year, cell, prop, numID,
                yearID, pairID, checkID, piID)

egridyr <- ebird_grid_df$yearID
egridpair <- ebird_grid_df$pairID
egridid <- ebird_grid_df$piID

ebird.pi <- array(NA, dim = c(n.years, max(egridpair), max(egridid)))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_grid_df)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, transect,
  #and pi ID of row i
  # populate that space in the array with the column in
  # the dataframe that corresponds to the pi data
  # for that yearxtransectxpi combo
  ebird.pi[egridyr[i], egridpair[i], egridid[i]] <- as.numeric(ebird_grid_df[i,4])
}

#set all values past the possible to 0 
#so the model doesn't break? i think i need this??
ebird.pi[is.na(ebird.pi)] <- 0


ebird.grid.array <- array(NA, dim = c(n.years, max(egridpair), max(egridid)))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_grid_df)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, transect,
  #and pi ID of row i
  # populate that space in the array with the column in
  # the dataframe that corresponds to the pi data
  # for that yearxtransectxpi combo
  ebird.grid.array[egridyr[i], egridpair[i], egridid[i]] <- as.numeric(ebird_grid_df[i,5])
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
                  Monsoon = Monsoon,
                  PinyonBA = PinyonBA,
                  #BBS loop
                  n.bbs.years = n.bbs.years,
                  n.bbs.trans = n.bbs.trans,
                  ObserverExp = ObserverExp,
                  n.bbs.points = n.bbs.points,
                  bbs.count = bbs.count,
                  bbs.pi = bbs.pi,
                  bbs.grid.array = bbs.grid.array,
                  bbs.year = bbs.year,
                  #ebird loop
                  n.ebird.pairs = n.ebird.pairs,
                  n.ebird.check = n.ebird.check,
                  ebird.count = ebird.count,
                  ebird.pi = ebird.pi,
                  ebird.grid.array = ebird.grid.array,
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


