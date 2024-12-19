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


ebird <- read_sf(here('data',
                       'ebird_data',
                       'cleaned_data',
                       'all_ebird_data_conefiltered.shp')) %>%
  dplyr::select(cellID,
                chckls_, #checklist ID
                obsrvtn_c, #observation count
                prtcl_t,#protocol type
                year, 
                tm_bsr_, #time observation started
                drtn_mn, #duration_minutes
                effrt__, #effort distance km
                nmbr_bs, #number observers
                geometry) %>%
  rename(blobID = cellID)

ebird_blobIDs <- read.csv(here('data',
                       'ebird_data',
                       'cleaned_data',
                       'ebird_cellIDlists.csv')) %>%
  dplyr::select(-X) %>%
  rename(blobID = cellID) %>%
  rename(cellID = cell)

cones <- read.csv(here('data',
                       'spatial_data',
                       'cleaned_data',
                       'cones_weighted_mean_blob.csv'))

temp <- read.csv(here('data',
                      'spatial_data',
                      'cleaned_data',
                      'temp_weighted_mean_blob.csv'))

ppt <- read.csv(here('data',
                     'spatial_data',
                     'cleaned_data',
                     'ppt_weighted_mean_blob.csv'))

monsoon <- read.csv(here('data',
                         'spatial_data',
                         'cleaned_data',
                         'monsoon_weighted_mean_blob.csv'))

pinyon <- read.csv(here('data',
                        'spatial_data',
                        'cleaned_data',
                        'pinyonBA_weighted_mean_blob.csv'))


#all blobs in blob lists to be able to get
#indexing
all_blobs <- ebird_blobIDs %>%
  distinct(year, blobnum) %>%
  group_by(year) %>%
  mutate(numID = 1:n()) %>%
  ungroup() %>%
  filter(year > 2009 & year < 2022)

# Some things about data --------------------------------------------------

#start by subsetting 2010-onward, since these are good ebird years
#but could consider other subsets of years down the road

# Subset years  ------------------------------------------

#subsetting without the last year given lack of overlap and forward
#projectiosn for cones for all data. Insetad of the last lag having
#~16.5% imputing, it now has ~10% and all the others have all
#their data
ebird3 <- ebird %>%
  filter(year > 2009 & year < 2022)

#combine and plot quickly
ebd4 <- ebird3 %>%
  dplyr::select(geometry, year) %>%
  mutate(datasource = 'ebird') %>%
  rename(Year = year)

ggplot(ebd4) +
  geom_sf(color = '#d8b365') +
  facet_wrap(~Year) +
  theme_bw() 

eb_2020 <- ebd4 %>%
  filter(Year == 2020)

ba_20 <- ba %>%
  filter(PinyonBA_sqftPerAc_2020 != 0)

ggplot() +
  geom_tile(data = ba_20, aes(x = x, y = y, fill = PinyonBA_sqftPerAc_2020))+
  geom_sf(data = eb_2020, color = "white", alpha = 0.6, shape = 2) +
  scale_fill_viridis_c()

ebird_blobIDs %>%
  group_by(blobnum) %>%
  tally() %>%
  ggplot() +
  geom_histogram(aes(x = n)) +
  labs(title = "Number of grids per blob",
       x = "Number of grids",
       y = "Number of blobs")

ebird_blobIDs %>%
  group_by(blobnum) %>%
  tally() %>%
  summarise(mean = mean(n),
            min = min(n),
            max = max(n))

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
n.years <- length(2010:2021)

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
  mutate(StartTime = scale(tm_bsr_),
         Duration = scale(drtn_mn),
         Distance = scale(effrt__),
         NumObservers = scale(nmbr_bs)) 
  #might need to get a random effect of observer, but let's wait on that 
  #for now and see how it goes with just what we have currently

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
  SurveyType[yr2[i], blob[i], check[i]] <- as.numeric(ebird_index_df[i,14])
}

StartTime   <- array(data = NA, dim = c(n.years, 
                                        max(n.blobs),
                                        max(n.ebird.check, na.rm = T)))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  StartTime[yr2[i], blob[i], check[i]] <- as.numeric(ebird_index_df[i,15])
}

Duration <- array(data = NA, dim = c(n.years, 
                                     max(n.blobs),
                                     max(n.ebird.check, na.rm = T)))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  Duration[yr2[i], blob[i], check[i]] <- as.numeric(ebird_index_df[i,16])
}

Distance <- array(data = NA, dim = c(n.years, 
                                     max(n.blobs),
                                     max(n.ebird.check, na.rm = T)))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  Distance[yr2[i],blob[i], check[i]] <- as.numeric(ebird_index_df[i,17])
}

NumObservers <- array(data = NA, dim = c(n.years, 
                                         max(n.blobs),
                                         max(n.ebird.check, na.rm = T)))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  NumObservers[yr2[i],blob[i],check[i]] <- as.numeric(ebird_index_df[i,18])
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
                  n.lag = n.lag,
                  n.clag = n.clag,
                  Cone = Cone,
                  Temp = Temp,
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
                  NumObservers = NumObservers)

saveRDS(data_list, here('data',
                        'jags_input_data',
                        'ebird_data_list_nospuncert.RDS'))

inits_list <- list(list(N = N),
                   list(N = N),
                   list(N = N))

saveRDS(inits_list, here('data',
                        'jags_input_data',
                        'ebird_init_list_nospuncert.RDS'))


