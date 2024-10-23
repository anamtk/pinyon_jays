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
                       'all_ebird_data_conefiltered.shp'))

ebird_gridIDs <- read.csv(here('data',
                               'ebird_data',
                               'cleaned_data',
                               'ebird_cellIDlists.csv'))

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
all_cells <- ebird_gridIDs %>%
  distinct(year, cell) %>%
  group_by(year) %>%
  mutate(numID = 1:n()) %>%
  ungroup() %>%
  filter(year > 2009 & year < 2022)

#get unique cells across all years - we will be
#creating covariate dataframes with data for every
#possible cell in the df for all years
unique_cells <- all_cells %>%
  distinct(cell)
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

# Data objects for jags ---------------------------------------------------


# Latent N loop -----------------------------------------------------------

#anywhere BA or Cones == 0, probs 0 because at end of 
#range - set to 0 instead of imputing in the model

#Total latent N:
n.grids <- all_cells %>%
  group_by(year) %>%
  summarise(n.grids = max(numID)) %>%
  dplyr::select(n.grids) %>%
  as_vector()

#NOTE: Reducing years to not include 2022 of data so that we 
#can better represent the cone data - when i have the full
#dataset, ~16.5% of the final lag are being imputed for cone data
#and i'm wondering if that's why estimates have been weird.
n.years <- length(2010:2021)

#Evnrionmental covariates
#cones:
cones2 <- cones %>%
  #get only cell IDs that overlap with bird data
  filter(cell %in% unique_cells$cell) %>%
  dplyr::select(cell, X2000:X2023) %>%
  pivot_longer(X2000:X2023,
               names_to = "year",
               values_to = "cones") %>%
  mutate(year = str_sub(year, 2, nchar(year))) %>%
  mutate(year = as.numeric(year)) %>%
  #year 1 == 2010. double check this is what we want
  mutate(yearID = as.numeric(as.factor(year)) - 10) %>%
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
  dplyr::select(yearID, year, cell, cones_l1:cones_l3, cones_0, 
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
                           lag == 'cones_n3' ~ 7)) %>%
  #now get the cells and years that match up to the 
  #numIDs
  full_join(all_cells, by = c("cell", "year")) %>%
  filter(!is.na(numID))

cones2 %>%
  summarise(na = sum(is.na(cones)),
            notna = sum(!is.na(cones)),
            prop = na/(na+notna))

#number of lags to consider
n.lag <- max(cones2$lagID, na.rm = T) #number of lags for each covariate

#now, generate IDs for the for loop where
# we will populate the array
conegrid <- cones2$numID
coneyear <- cones2$yearID
conelag <- cones2$lagID

#make a blank array
Cone <- array(data = NA, dim = c( n.years, max(n.grids),n.lag))

#fill taht array based on the values in those columns
for(i in 1:dim(cones2)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, grid,
  #and lag ID of row i 
  # populate that space in the array with the column in
  # the dataframe that corresponds to the cone data
  # for that yearxgridxlag combo
  Cone[coneyear[i], conegrid[i], conelag[i]] <- as.numeric(cones2[i,5])
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

#temperature
temp2 <- temp %>%
  #get only cell IDs that overlap with bird data
  filter(cellID %in% unique_cells$cell) %>%
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
  #"current" season is season 1, the season when data started
  #being collected
  filter(season == 1) %>%
  dplyr::select(-season) %>%
  dplyr::select(yearID, cellID, year, temp_l1:temp_l13) %>%
  pivot_longer(temp_l1:temp_l13,
               names_to = 'lag',
               values_to = "temp") %>%
  mutate(lagID = str_sub(lag, 7, nchar(lag))) %>%
  mutate(lagID = as.numeric(lagID)) %>%
  dplyr::select(yearID,  year, cellID, lagID, temp) %>%
  filter(yearID < 13) %>%
  #now get the cells and years that match up to the 
  #numIDs
  full_join(all_cells, by = c('cellID' = "cell", "year")) %>%
  filter(!is.na(numID))
 
temp2 %>%
  summarise(na = sum(is.na(temp)),
            notna = sum(!is.na(temp)),
            prop = na/(na+notna)) 

#lag index
n.clag <- max(temp2$lagID)

#now, generate IDs for the for loop where
# we will populate the array
tempgrid <- temp2$numID
tempyear <- temp2$yearID
templag <- temp2$lagID

#make a blank array
Temp <- array(data = NA, dim = c(n.years,max(n.grids), n.clag))

#fill taht array based on the values in those columns
for(i in 1:dim(temp2)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, grid,
  #and lag ID of row i 
  # populate that space in the array with the column in
  # the dataframe that corresponds to the data
  # for that yearxgridxlag combo
  Temp[tempyear[i],tempgrid[i],  templag[i]] <- as.numeric(temp2[i,5])
}

#PPT
ppt2 <- ppt %>%
  #get only cell IDs that overlap with bird data
  filter(cellID %in% unique_cells$cell) %>%
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
  #"current" season is season 1, the season when data started
  #being collected
  filter(season == 1) %>%
  dplyr::select(-season) %>%
  dplyr::select(yearID, cellID, year, ppt_l1:ppt_l13) %>%
  pivot_longer(ppt_l1:ppt_l13,
               names_to = 'lag',
               values_to = "ppt") %>%
  mutate(lagID = str_sub(lag, 6, nchar(lag))) %>%
  mutate(lagID = as.numeric(lagID)) %>%
  dplyr::select(yearID, cellID, year, lagID, ppt) %>%
  filter(yearID < 13)%>%
  #now get the cells and years that match up to the 
  #numIDs
  full_join(all_cells, by = c('cellID' = "cell", "year")) %>%
  filter(!is.na(numID))

ppt2 %>%
  summarise(na = sum(is.na(ppt)),
            notna = sum(!is.na(ppt)),
            prop = na/(na+notna)) 

#now, generate IDs for the for loop where
# we will populate the array
pptgrid <- ppt2$numID
pptyear <- ppt2$yearID
pptlag <- ppt2$lagID

#make a blank array
PPT <- array(data = NA, dim = c(n.years,max(n.grids),  n.clag))

#fill taht array based on the values in those columns
for(i in 1:dim(ppt2)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, grid,
  #and lag ID of row i 
  # populate that space in the array with the column in
  # the dataframe that corresponds to the data
  # for that yearxgridxlag combo
  PPT[ pptyear[i], pptgrid[i],pptlag[i]] <- as.numeric(ppt2[i,5])
}


#monsoon [t,i]
Monsoon <- mons %>%
  group_by(cellID) %>%
  summarise(monsoon = mean(PRISM_ppt_30yr_normal_800mM4_07_bil, na.rm = T)) %>%
  ungroup() %>%
  filter(cellID %in% unique_cells$cell) %>%
  #filter(!is.na(numID)) %>%
  mutate(monsoon = scale(monsoon)) %>%
  full_join(all_cells, by = c("cellID" = "cell")) %>%
  dplyr::select(monsoon, year, numID) %>%
  arrange(numID) %>%
  pivot_wider(names_from = numID,
              values_from = monsoon) %>%
  arrange(year) %>%
  column_to_rownames(var = "year") %>%
  as.matrix()

# numID <- as.data.frame(1:6673) %>%
#   rename(numID = `1:6673`)


#basal area
PinyonBA <- ba %>%
  filter(cell %in% unique_cells$cell) %>%
  #filtering out 2022 for now
  dplyr::select(cell, PinyonBA_sqftPerAc_2010:PinyonBA_sqftPerAc_2021) %>%
  pivot_longer(PinyonBA_sqftPerAc_2010:PinyonBA_sqftPerAc_2021,
               names_to = 'year',
               values_to = 'ba') %>%
  mutate(year = str_sub(year, (nchar(year)-3), nchar(year)),
         year = as.numeric(year)) %>%
  full_join(all_cells, by = c("cell" = "cell", "year")) %>%
  filter(!is.na(numID)) %>%
  #19 total NAs
  #replace_na(list(ba = 0)) %>%
  mutate(ba = scale(ba)) %>%
  dplyr::select(-cell) %>%
  arrange(numID) %>%
  pivot_wider(names_from = numID, 
              values_from = ba) %>%
  arrange(year) %>%
  column_to_rownames(var = 'year') %>%
  as.matrix()

#very small # are NA - we can impute
sum(is.na(PinyonBA))/(sum(is.na(PinyonBA) + sum(!is.na(PinyonBA))))

# eBIRD loop --------------------------------------------------------------

#will need to change the pi and indexing array to be 
#one less dimension

#eBIRD:

#vector of # pairs/year (sometimes singleton)
n.ebird.check <- ebird3 %>%
  distinct(year, chckls_) %>%
  group_by(year) %>%
  tally() %>%
  ungroup() %>%
  dplyr::select(n) %>%
  as_vector()

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
                yrID, geometry) %>%
  group_by(yrID) %>%
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
check <- ebird_index_df$checkID #get a checklist ID for each iteration of the loop

#make a blank matrix
ebird.count <- matrix(data = NA, nrow = n.years, ncol = max(n.ebird.check))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, transect,
  #and stop ID of row i 
  # populate that space in the array with the column in
  # the dataframe that corresponds to the count data
  # for that yearxtransectxstop combo
  ebird.count[yr2[i], check[i]] <- as.numeric(ebird_index_df[i,2])
}

#do the same for the covariates below:

SurveyType <- matrix(data = NA, nrow = n.years, ncol = max(n.ebird.check))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  SurveyType[yr2[i], check[i]] <- as.numeric(ebird_index_df[i,12])
}

StartTime <- matrix(data = NA, nrow = n.years, ncol = max(n.ebird.check))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  StartTime[yr2[i], check[i]] <- as.numeric(ebird_index_df[i,13])
}

Duration <- matrix(data = NA, nrow = n.years, ncol = max(n.ebird.check))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  Duration[yr2[i], check[i]] <- as.numeric(ebird_index_df[i,14])
}

Distance <- matrix(data = NA, nrow = n.years, ncol = max(n.ebird.check))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  Distance[yr2[i], check[i]] <- as.numeric(ebird_index_df[i,15])
}

NumObservers <- matrix(data = NA, nrow = n.years, ncol = max(n.ebird.check))
#fill taht array based on the values in those columns
for(i in 1:dim(ebird_index_df)[1]){ #dim[1] = n.rows
  NumObservers[yr2[i],check[i]] <- as.numeric(ebird_index_df[i,16])
}

ebird_index2 <- ebird_index_df %>%
  dplyr::select(chckls_, checkID)
  
#I need to link up metadata on the 
#"replicate" each checklist represents for each
#cell, and then the indexing will be [i,t,]
ebird_grid_df <- ebird_gridIDs %>%
  filter(year >= 2010 &  year < 2022) %>%
  #get cell numIDs for coding for jags
  left_join(all_cells, by = c("cell", 'year')) %>%
  mutate(yearID = as.numeric(as.factor(year))) %>%
  #get yrtransID index (which transect in year t are we )
  left_join(ebird_index2, by = c("chckls_")) %>%
  #filter(!is.na(pairID)) %>%
  group_by(yearID, checkID) %>%
  mutate(piID = 1:n()) %>%
  ungroup() %>%
  dplyr::select(chckls_, year, cell, prop, numID,
                yearID, checkID, piID, checkID)

egridyr <- ebird_grid_df$yearID
egridcheck <- ebird_grid_df$checkID
egridid <- ebird_grid_df$piID

ebird.pi <- array(NA, dim = c(n.years, max(egridcheck), max(egridid)))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_grid_df)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, transect,
  #and pi ID of row i
  # populate that space in the array with the column in
  # the dataframe that corresponds to the pi data
  # for that yearxtransectxpi combo
  ebird.pi[egridyr[i], egridcheck[i], egridid[i]] <- as.numeric(ebird_grid_df[i,4])
}

#set all values past the possible to 0 
#so the model doesn't break? i think i need this??
#ebird.pi[is.na(ebird.pi)] <- 0


ebird.grid.array <- array(NA, dim = c(n.years, max(egridcheck), max(egridid)))

#fill taht array based on the values in those columns
for(i in 1:dim(ebird_grid_df)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the year, transect,
  #and pi ID of row i
  # populate that space in the array with the column in
  # the dataframe that corresponds to the pi data
  # for that yearxtransectxpi combo
  ebird.grid.array[egridyr[i], egridcheck[i], egridid[i]] <- as.numeric(ebird_grid_df[i,5])
}

#matrix t,i
n.ebird.cells <- ebird_grid_df %>%
  distinct(yearID, checkID, piID) %>%
  group_by(yearID, checkID) %>%
  tally() %>%
  pivot_wider(names_from = checkID,
              values_from = n) %>%
  column_to_rownames(var = "yearID") %>%
  as.matrix()
# Values for initials -----------------------------------------------------

#need a starting value for N[i,t]

#maxb <- max(bbs.count, na.rm = T)
nmax <- max(ebird.count, na.rm=T)
#nmax <- max(c(maxb, maxe))

ndf <- all_cells %>%
  mutate(yearID = as.numeric(as.factor(year)))

N <- matrix(NA, nrow = n.years, ncol =max(n.grids))

nyr <- ndf$yearID
nid <- ndf$numID

for(i in 1:dim(ndf)[1]){ #dim[1] = n.rows
  N[nyr[i], nid[i]] <- nmax
}

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
                  #ebird loop
                  n.ebird.check = n.ebird.check,
                  ebird.count = ebird.count,
                  ebird.pi = ebird.pi,
                  ebird.grid.array = ebird.grid.array,
                  n.ebird.cells = n.ebird.cells,
                  SurveyType = SurveyType,
                  StartTime = StartTime,
                  Duration = Duration,
                  Distance = Distance,
                  NumObservers = NumObservers)

saveRDS(data_list, here('data',
                        'jags_input_data',
                        'ebird_data_list.RDS'))

inits_list <- list(list(N = N),
                   list(N = N),
                   list(N = N))

saveRDS(inits_list, here('data',
                        'jags_input_data',
                        'ebird_init_list.RDS'))


