#Prep data for model
#Ana Miller-ter Kuile
#July 3, 2024

#prep all datasets for the model

#NEED TO ADD:
#blobArea[i,t] - multiplies by lambda

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
                  'sf', 'patchwork')

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

tmean <- read.csv(here('data',
                       'spatial_data',
                       'cleaned_data',
                       'tmean_weighted_mean_blob.csv'))

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

# Checks ------------------------------------------------------------------

#get presence-only checklists
ebird_pres <- ebird3 %>%
  filter(obsrvtn_c > 0) 

#get blob IDs to talk to covariate data for 
#presence datasets
blobIDs_pres <- ebird_blobIDs %>%
  distinct(blobID, blobnum, year) %>%
  left_join(ebird_pres, by = c('blobID', 'year')) %>%
  filter(!is.na(chckls_))

blobIDs_pres2 <- blobIDs_pres %>%
  dplyr::select(blobnum, year)

# Spatial clustering of presence? -----------------------------------------

#ANSWER: NO

#- are the presence checklists clustered in space?
states <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE)) %>%
  filter(ID %in% c('arizona', 'colorado', 
                   'utah', 'new mexico'))

#presence checklists are NOT clustered in space
ggplot() +
  geom_sf(data = states, fill = "white") +
  geom_sf(data = ebird_pres) +
  facet_wrap(~year)


# Look at covariate distributions -----------------------------------------
#- are the presence checklists clustered in covariate values? or
## are they variable in covariate values? 
#- check distribution of covariates as well

#cones
pres_cone <- cones %>%
  filter(blobnum %in% (blobIDs_pres2$blobnum)) %>%
  mutate(type = "presence")

all_cone <- cones %>%
  mutate(type = "all") %>%
  bind_rows(pres_cone)
  
#seems ok:
ggplot(all_cone, aes(x = cones, fill = type)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~lag)

#tmax
pres_tmax <- temp %>%
  filter(blobnum %in% (blobIDs_pres2$blobnum)) %>%
  mutate(type = "presence")

all_tmax <- temp %>%
  mutate(type = "all") %>%
  bind_rows(pres_tmax)

#maybe a little "squinched?"
ggplot(all_tmax, aes(x = temp, fill = type)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~lag)

#ppt
pres_ppt <- ppt %>%
  filter(blobnum %in% (blobIDs_pres2$blobnum)) %>%
  mutate(type = "presence")

ppt_all <- ppt %>%
  mutate(type = "all") %>%
  bind_rows(pres_ppt)
 
#seems okay??
ggplot(ppt_all, aes(x = ppt, fill = type)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~lag) 

#tmean
pres_tmean <- tmean %>%
  filter(blobnum %in% (blobIDs_pres2$blobnum)) %>%
  mutate(type = "presence")

all_tmean <- tmean %>%
  mutate(type = "all") %>%
  bind_rows(pres_tmean)
  
#maybe a little "squinched?"
ggplot(all_tmean, aes(x = temp, fill = type)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~lag) 

#monsoon
pres_mon <- monsoon %>%
  filter(blobnum %in% (blobIDs_pres2$blobnum)) %>%
  mutate(type = "presence")

all_mon <- monsoon %>%
  mutate(type = "all") %>%
  bind_rows(pres_mon)

#maybe a little shifted to places
#with more monsoon? 
ggplot(all_mon, aes(x = wt, fill = type))+
  geom_density(alpha = 0.4)

#pinyon ba
pres_ba <- pinyon %>%
  filter(blobnum %in% (blobIDs_pres2$blobnum)) %>%
  mutate(type = "presence")

all_ba <- pinyon %>%
  mutate(type = "all") %>%
  bind_rows(pres_ba)

#seems okay?
ggplot(all_ba, aes(x = wt, fill = type))+
  geom_density(alpha = 0.4)

#- check priors - do they mesh with data?
 
