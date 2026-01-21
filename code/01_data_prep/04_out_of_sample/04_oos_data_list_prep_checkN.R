#Prep data for model
#Ana Miller-ter Kuile
#July 3, 2024

#prep all datasets for the model

#NEED TO ADD:
#blobArea[t,i] - multiplies by lambda

#ALSO - NEED TO CHECK DATA DISTRIBUTIONS FOR:
#- are the presence checklists clustered in space/covariate values? or
## are they variable in covariate values? 
#- check priors - do they mesh with data?
#- are presence records geographically distributed narrowly? Could 
## I run a model with only those areas and use that as initials?
#- check distribution of covariates. 

#POTENTIALLY:
#- fix c0 at a reasonable value of ~3% presence and see if it helps


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 'data.table',
                  'corrplot',
                  'sf')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())

# Steps -------------------------------------------------------------------

# Load datasets -----------------------------------------------------------


ebird <- read_sf(here('data',
                      'ebird_data',
                      'cleaned_data',
                      "05_oos",
                      'all_oos_ebird_data_conefiltered.shp')) %>%
  dplyr::select(cellID,
                chckls_, #checklist ID
                obsrvtn_c, #observation count
                prtcl_t,#protocol type
                year, 
                tm_bsr_, #time observation started
                drtn_mn, #duration_minutes
                effrt__, #effort distance km
                nmbr_bs, #number observers
                buffr_m, #size of individual checklist buffer
                geometry) %>%
  rename(blobID = cellID)

ebird_blobIDs <- read.csv(here('data',
                               'ebird_data',
                               'cleaned_data',
                               '05_oos',
                               'ebird_oos_cellIDlists.csv')) %>%
  dplyr::select(-X) %>%
  rename(blobID = cellID) %>%
  rename(cellID = cell)

cones <- read.csv(here('data',
                       'spatial_data',
                       'cleaned_data',
                       '03_oos',
                       'oos_cones_weighted_mean_blob.csv'))

temp <- read.csv(here('data',
                      'spatial_data',
                      'cleaned_data',
                      '03_oos',
                      'oos_temp_weighted_mean_blob.csv'))

tmean <- read.csv(here('data',
                       'spatial_data',
                       'cleaned_data',
                       '03_oos',
                       'oos,tmean_weighted_mean_blob.csv'))

ppt <- read.csv(here('data',
                     'spatial_data',
                     'cleaned_data',
                     '03_oos',
                     'oos_ppt_weighted_mean_blob.csv'))

monsoon <- read.csv(here('data',
                         'spatial_data',
                         'cleaned_data',
                         '03_oos',
                         'oos_monsoon_weighted_mean_blob.csv'))

pinyon <- read.csv(here('data',
                        'spatial_data',
                        'cleaned_data',
                        '03_oos',
                        'oos_pinyonBA_weighted_mean_blob.csv'))


#all blobs in blob lists to be able to get
#indexing
all_blobs <- ebird_blobIDs %>%
  distinct(year, blobnum, area) %>%
  group_by(year) %>%
  mutate(numID = 1:n()) %>%
  ungroup() %>%
  filter(year > 2009 & year < 2023)

# Some things about data --------------------------------------------------

#start by subsetting 2010-onward, since these are good ebird years
#but could consider other subsets of years down the road

# Subset years  ------------------------------------------

#subsetting without the last year given lack of overlap and forward
#projectiosn for cones for all data. Insetad of the last lag having
#~16.5% imputing, it now has ~10% and all the others have all
#their data
ebird3 <- ebird %>%
  filter(year > 2009 & year < 2023)

# Data objects for jags ---------------------------------------------------


# Latent N loop -----------------------------------------------------------

#anywhere BA or Cones == 0, probs 0 because at end of 
#range - set to 0 instead of imputing in the model

#Total latent N:
n.blobs <- all_blobs %>%
  group_by(year) %>%
  summarise(n.blobs = n()) %>%
  dplyr::select(n.blobs) %>%
  as_vector()

#NOTE: Reducing years to not include 2022 of data so that we 
#can better represent the cone data - when i have the full
#dataset, ~16.5% of the final lag are being imputed for cone data
#and i'm wondering if that's why estimates have been weird.
n.years <- length(2010:2022)

#area to multiply lambda by [year, blob]
blobArea <- all_blobs %>%
  dplyr::select(numID, area,year) %>%
  pivot_wider(names_from = numID,
              values_from = area) %>%
  arrange(year) %>%
  column_to_rownames(var = "year") %>%
  as.matrix()
  

#Evnrionmental covariates
#cones:
cones2 <- cones %>%
  filter(blobnum %in% all_blobs$blobnum) %>%
  mutate(cones = scale(cones)) %>%
  mutate(yearID = as.numeric(as.factor(year))) %>%
  left_join(all_blobs, by = c("year", "blobnum"))

cones2 %>%
  summarise(na = sum(is.na(cones)),
            notna = sum(!is.na(cones)),
            prop = na/(na+notna))

#number of lags to consider
n.lag <- max(cones2$lag, na.rm = T) #number of lags for each covariate

#now, generate IDs for the for loop where
# we will populate the array
coneblob <- cones2$numID
coneyear <- cones2$yearID
conelag <- cones2$lag

#make a blank array
Cone <- array(data = NA, dim = c( n.years, max(n.blobs),n.lag))

#fill taht array based on the values in those columns
for(i in 1:dim(cones2)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, grid,
  #and lag ID of row i 
  # populate that space in the array with the column in
  # the dataframe that corresponds to the cone data
  # for that yearxgridxlag combo
  Cone[coneyear[i], coneblob[i], conelag[i]] <- as.numeric(cones2[i,4])
}

#seasons for climate variables::
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

temp2 <- temp %>%
  filter(blobnum %in% all_blobs$blobnum) %>%
  #trying scaling by season
  group_by(season) %>%
  mutate(temp = scale(temp)) %>%
  ungroup() %>%
  mutate(yearID = as.numeric(as.factor(year))) %>%
  left_join(all_blobs, by = c("year", "blobnum"))

temp2 %>%
  summarise(na = sum(is.na(temp)),
            notna = sum(!is.na(temp)),
            prop = na/(na+notna)) 

#lag index
n.clag <- max(temp2$lag)

#now, generate IDs for the for loop where
# we will populate the array
tempblob <- temp2$numID
tempyear <- temp2$yearID
templag <- temp2$lag

#make a blank array
Temp <- array(data = NA, dim = c(n.years,max(n.blobs), n.clag))

#fill taht array based on the values in those columns
for(i in 1:dim(temp2)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, grid,
  #and lag ID of row i 
  # populate that space in the array with the column in
  # the dataframe that corresponds to the data
  # for that yearxgridxlag combo
  Temp[tempyear[i],tempblob[i],  templag[i]] <- as.numeric(temp2[i,5])
}

tmean2 <- tmean %>%
  filter(blobnum %in% all_blobs$blobnum) %>%
  #trying scaling by season
  group_by(season) %>%
  mutate(temp = scale(temp)) %>%
  ungroup() %>%
  mutate(yearID = as.numeric(as.factor(year))) %>%
  left_join(all_blobs, by = c("year", "blobnum"))

tmean2 %>%
  summarise(na = sum(is.na(temp)),
            notna = sum(!is.na(temp)),
            prop = na/(na+notna)) 

#lag index
n.clag <- max(tmean2$lag)

#now, generate IDs for the for loop where
# we will populate the array
tmeanblob <- tmean2$numID
tmeanyear <- tmean2$yearID
tmeanlag <- tmean2$lag

#make a blank array
Tmean <- array(data = NA, dim = c(n.years,max(n.blobs), n.clag))

#fill taht array based on the values in those columns
for(i in 1:dim(tmean2)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, grid,
  #and lag ID of row i 
  # populate that space in the array with the column in
  # the dataframe that corresponds to the data
  # for that yearxgridxlag combo
  Tmean[tmeanyear[i],tmeanblob[i],  tmeanlag[i]] <- as.numeric(tmean2[i,5])
}

ppt2 <- ppt %>%
  filter(blobnum %in% all_blobs$blobnum) %>%
  #trying scaling by season
  group_by(season) %>%
  mutate(ppt = scale(ppt )) %>%
  ungroup() %>%
  mutate(yearID = as.numeric(as.factor(year))) %>%
  left_join(all_blobs, by = c("year", "blobnum"))

ppt2 %>%
  summarise(na = sum(is.na(ppt)),
            notna = sum(!is.na(ppt)),
            prop = na/(na+notna)) 

#now, generate IDs for the for loop where
# we will populate the array
pptblob <- ppt2$numID
pptyear <- ppt2$yearID
pptlag <- ppt2$lag

#make a blank array
PPT <- array(data = NA, dim = c(n.years,max(n.blobs),  n.clag))

#fill taht array based on the values in those columns
for(i in 1:dim(ppt2)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, grid,
  #and lag ID of row i 
  # populate that space in the array with the column in
  # the dataframe that corresponds to the data
  # for that yearxgridxlag combo
  PPT[ pptyear[i], pptblob[i],pptlag[i]] <- as.numeric(ppt2[i,5])
}

#MONSOON
#t,i
Monsoon <- monsoon %>%
  filter(blobnum %in% all_blobs$blobnum) %>%
  mutate(monsoon = scale(wt)) %>%
  mutate(yearID = as.numeric(as.factor(year))) %>%
  left_join(all_blobs, by = c("year", "blobnum")) %>%
  dplyr::select(yearID, numID, monsoon) %>%
  pivot_wider(names_from = 'numID',
              values_from = 'monsoon' ) %>%
  arrange(yearID) %>%
  column_to_rownames(var = 'yearID') %>%
  as.matrix()

#pinyon basal area
PinyonBA <- pinyon %>%
  filter(blobnum %in% all_blobs$blobnum) %>%
  mutate(pinyon = scale(wt)) %>%
  mutate(yearID = as.numeric(as.factor(year))) %>%
  left_join(all_blobs, by = c("year", "blobnum")) %>%
  dplyr::select(yearID, numID, pinyon) %>%
  pivot_wider(names_from = 'numID',
              values_from = 'pinyon' ) %>%
  arrange(yearID) %>%
  column_to_rownames(var = 'yearID') %>%
  as.matrix()

#very small # are NA - we can impute
sum(is.na(PinyonBA))/(sum(is.na(PinyonBA) + sum(!is.na(PinyonBA))))

# eBIRD loop --------------------------------------------------------------

ebird_blobs <- ebird_blobIDs %>%
  dplyr::select(blobID, year, blobnum) %>%
  left_join(all_blobs, by = c("year", "blobnum")) %>%
  filter(!is.na(numID)) %>%
  distinct(blobID, year, blobnum, numID)
#will need to change the pi and indexing array to be 
#one less dimension

#eBIRD:
#matrix of ebird checklists in year t in blob i
n.ebird.check <- as.data.frame(ebird3) %>%
  dplyr::select(blobID, chckls_, year) %>%
  left_join(ebird_blobs, by = c('blobID', "year")) %>%
  group_by(numID, year) %>%
  tally() %>%
  pivot_wider(names_from = numID, 
              values_from = n) %>%
  column_to_rownames(var= 'year') %>%
  as.matrix()
  

#create a dataframe that we can referecne for all the loops below
ebird_index_df <- as.data.frame(ebird3) %>%
  ungroup() %>%
  left_join(ebird_blobs, by =c("blobID", "year")) %>%
  mutate(yrID = as.numeric(as.factor(year))) %>%
  dplyr::select(blobnum,
                numID,
                yrID, 
                year, 
                chckls_, #checklist ID
                obsrvtn_c, #observation count
                prtcl_t,#protocol type
                tm_bsr_, #time observation started
                drtn_mn, #duration_minutes
                effrt__, #effort distance km
                nmbr_bs, #number observers
                buffr_m,
                geometry) %>%
  group_by(yrID, blobnum) %>%
  mutate(checkID = 1:n()) %>%
  ungroup() %>%
  mutate(SurveyType = case_when(prtcl_t == "Stationary" ~ 1,
                                prtcl_t == "Traveling" ~ 2)) %>%
  #duration  
  #distance
  #NumObservers
  #time started (hours since midnight)
  mutate(Speed = scale(effrt__/drtn_mn),
         StartTime = scale(tm_bsr_),
         Duration = scale(drtn_mn),
         Distance = scale(effrt__),
         NumObservers = scale(nmbr_bs)) 
  #might need to get a random effect of observer, but let's wait on that 
  #for now and see how it goes with just what we have currently

n.checklists <- nrow(ebird_index_df)
#now, generate IDs for the for loop where
# we will populate the array
yr2 <- ebird_index_df$yrID #get a yearID for each iteration of the loop
blob <- ebird_index_df$numID
check <- ebird_index_df$checkID #get a checklist ID for each iteration of the loop

#make a blank array t,i,r
ebird.count <- array(data = NA, dim = c(n.years, 
                                      max(n.blobs),
                                      max(n.ebird.check, na.rm = T)))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, transect,
  #and stop ID of row i 
  # populate that space in the array with the column in
  # the dataframe that corresponds to the count data
  # for that yearxtransectxstop combo
  ebird.count[yr2[i], blob[i], check[i]] <- as.numeric(ebird_index_df[i,6])
}

#do the same for the covariates below:

SurveyType  <- array(data = NA, dim = c(n.years, 
                                        max(n.blobs),
                                        max(n.ebird.check, na.rm = T)))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  SurveyType[yr2[i], blob[i], check[i]] <- as.numeric(ebird_index_df[i,15])
}

StartTime   <- array(data = NA, dim = c(n.years, 
                                        max(n.blobs),
                                        max(n.ebird.check, na.rm = T)))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  StartTime[yr2[i], blob[i], check[i]] <- as.numeric(ebird_index_df[i,16])
}

Duration <- array(data = NA, dim = c(n.years, 
                                     max(n.blobs),
                                     max(n.ebird.check, na.rm = T)))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  Duration[yr2[i], blob[i], check[i]] <- as.numeric(ebird_index_df[i,17])
}

Distance <- array(data = NA, dim = c(n.years, 
                                     max(n.blobs),
                                     max(n.ebird.check, na.rm = T)))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  Distance[yr2[i],blob[i], check[i]] <- as.numeric(ebird_index_df[i,18])
}

Speed <- array(data = NA, dim = c(n.years, 
                                     max(n.blobs),
                                     max(n.ebird.check, na.rm = T)))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  Speed[yr2[i],blob[i], check[i]] <- as.numeric(ebird_index_df[i,18])
}

NumObservers <- array(data = NA, dim = c(n.years, 
                                         max(n.blobs),
                                         max(n.ebird.check, na.rm = T)))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  NumObservers[yr2[i],blob[i],check[i]] <- as.numeric(ebird_index_df[i,19])
}

#area for each checklist: listArea[t,i,r]

listArea <- array(data = NA, dim = c(n.years, 
                                     max(n.blobs),
                                     max(n.ebird.check, na.rm = T)))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  listArea[yr2[i],blob[i],check[i]] <- as.numeric(pi*((ebird_index_df[i,12]))^2)
}

# Values for initials -----------------------------------------------------

#need a starting value for N[i,t]

#maxb <- max(bbs.count, na.rm = T)
nmax <- max(ebird.count, na.rm=T)
#nmax <- max(c(maxb, maxe))

ndf <- all_blobs %>%
  mutate(yearID = as.numeric(as.factor(year)))

N <- matrix(NA, nrow = n.years, ncol =max(n.blobs))

nyr <- ndf$yearID
nid <- ndf$numID

for(i in 1:dim(ndf)[1]){ #dim[1] = n.rows
  N[nyr[i], nid[i]] <- nmax
}

# Compile and export ------------------------------------------------------

data_list <- list(#latent N loop:
                  n.blobs = n.blobs,
                  n.years = n.years,
                  blobArea = blobArea,
                  n.lag = n.lag,
                  n.clag = n.clag,
                  Cone = Cone,
                  Temp = Temp,
                  Tmean = Tmean,
                  PPT = PPT,
                  Monsoon = Monsoon,
                  PinyonBA = PinyonBA,
                  #ebird loop
                  n.ebird.check = n.ebird.check,
                  ebird.count = ebird.count,
                  SurveyType = SurveyType,
                  StartTime = StartTime,
                  Duration = Duration,
                  Distance = Distance,
                  Speed = Speed,
                  NumObservers = NumObservers,
                  listArea = listArea,
                  #RMSE
                  n.checklists = n.checklists)

saveRDS(data_list, here('data',
                        'jags_input_data',
                        'oos',
                        'oos_ebird_data_list_nospuncert.RDS'))

# inits_list <- list(list(N = N),
#                    list(N = N),
#                    list(N = N))
# 
# saveRDS(inits_list, here('data',
#                         'jags_input_data',
#                         'ebird_init_list_nospuncert.RDS'))

saveRDS(ebird_index_df, 
        here('data',
             'ebird_data',
             'cleaned_data',
             '05_oos',
             'oos_ebird_check_blob_yr_ids.RDS'))


