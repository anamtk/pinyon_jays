#Ana Miller-ter Kuile
#make ebird and bbs raster data
#July 2, 2024

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "sf",  
                  "terra",
                  'readxl',
                  'sf',
                  'exactextractr',
                  'nngeo',
                  'sp')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

set.seed(1)
# Load data ---------------------------------------------------------------

# load ebird --------------------------------------------------------------

ebird <- read.csv(here('data',
                       'ebird_data',
                       'cleaned_data',
                       'all_ebird_data.csv'))

# load bbs ----------------------------------------------------------------

bbs <- read.csv(here('data',
                     'bbs_data',
                     'cleaned_data',
                     'pjay_data_co_nm_az_ut.csv'))

bbs_locations <- read.csv(here('data',
                               'bbs_data',
                               '2023Release_Nor',
                               'routes.csv'))

bbs_locations <- bbs_locations %>%
  dplyr::select(StateNum, Route, Latitude, Longitude)

bbs2 <- bbs %>%
  left_join(bbs_locations, by = c("StateNum", "Route"))

# Load cone dataset -------------------------------------------------------

cones <- terra::rast(here("data",
                          "spatial_data", 
                          "masting_data",
                          "full_pied_masting.tif"))
#plot(cones)

# Set the CRS for both bird datasets -------------------------------------------

#these data are in WGS84, so I'll create this crs to call later
crs_wgs <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

crs_wgs2 <- "+proj=utm +zone=42N +datum=WGS84 +units=km"

crs_albers <- '+proj=longlat +datum=NAD83 +no_defs +type=crs'


# Make spatial ------------------------------------------------------------

#get m radius for buffer.
#some are zero (stationary points)
#so maybe set a value for these, maybe 500m??
#get the others to be half the effort distance - assumption: walk to and
#from a vehicle
ebird2 <- ebird %>%
  rowwise() %>%
  mutate(buffer_m = (effort_distance_km*1000)/2) %>%
  ungroup()%>%
  mutate(buffer_m = case_when(buffer_m == 0 ~ 500,
                              TRUE ~ buffer_m))

#convert all ebird data to an sf object
ebird_spatial <- st_as_sf(ebird2, coords = c("longitude", "latitude"),
                          crs = st_crs(crs_albers))

#convert all bbs data to an sf object
bbs_spatial <- st_as_sf(bbs2, coords = c("Longitude", "Latitude"),
                        crs = st_crs(crs_albers))

# Extract cell IDs from raster ------------------------------------------------------------

##EBIRD
# Use template raster and ebird pts to get cell IDs
ebird_cellIDs <- terra::extract(cones, vect(ebird_spatial), cells = T)

# Add as column to ebird data
ebird_spatial$cellID <- ebird_cellIDs$cell

#some ebird observations don't have masting data cells, so I'll remove those
ebird_spatial <- ebird_spatial %>%
  filter(!is.na(cellID))

##BBS
bbs_cellIDs <- terra::extract(cones, vect(bbs_spatial), cells = T)

bbs_spatial$cellID <- bbs_cellIDs$cell

# Turn cone raster into df with cell IDs ---------------------------------------

cone_df <- terra::as.data.frame(cones, 
                                xy = TRUE,
                                cells = TRUE)

# Filter datasets for cells with cone data --------------------------------

ids <- cone_df %>%
  distinct(cell)

ebird_spatial2 <- ebird_spatial %>%
  filter(cellID %in% ids$cell)

bbs_spatial2 <- bbs_spatial %>%
  filter(cellID %in% ids$cell)

# Subset ebird data -------------------------------------------------------
#this paper: 
#https://www.nature.com/articles/s41598-022-23603-0
#subset one 1 and one zero dataset from each cell
#the chapter from cornell lab subsets "10 per site", not specifying whether
#they're 1 or 0 checklists
#https://cornelllabofornithology.github.io/ebird-best-practices/occupancy.html#occupancy-intro

#get yes-no checklist pairs in a year that are closest together in space

#stratify these pairs AND any solo yes or no points in space with some
#buffer distance to reduce spatial autocorrelation (maybe just some
#scale similar to the 4km grid cells of the covariates)

#types of points we are working with:
# pairs of yes-no observations
# solo no observations
# solo yes observations
# all spatially stratified... 

#find pairs that are closest together then:
#calculate the "middle" point of these as our "x-y" for that "pair" 
# to create a buffer around
#use those middle points to stratify all points (pairs, solo yes and no)???
#what i was doing before to subset them:
# ebird_yes <- ebird_spatial2 %>%
#   filter(observation_count > 0) %>%
#   group_by(year, cellID) %>%
#   slice_sample(n = 1) %>%
#   ungroup()

ebird_yes <- ebird_spatial2 %>%
  filter(observation_count > 0) %>%
  filter(year > 2009) %>%
  ungroup()

ebird_no <- ebird_spatial2 %>%
  filter(observation_count == 0) %>%
  filter(year > 2009)

group_function <- function(year){
  
  set.seed(1)
  yesdf2 <- ebird_yes %>%
    filter(year == {{year}}) 
  
  nodf2 <- ebird_no %>%
    filter(year == {{year}})
  
  #find nearest neighbors for each of the yes-no combos
  #but only those that are within 1km
  nn <- st_nn(x = yesdf2,
              y = nodf2,
              k = 1, 
              maxdist = 1000,
              returnDist = T,
              sparse = T)
  
  #get pairs based on row numbers to keep but then
  #subsample from
  pairIDs <- map_df(nn$nn, ~as.data.frame(.x), .id="id") %>%
    rename(yes = id,
           no = `.x`) %>%
    #group by "no" since multiple match the "no" category
    group_by(no) %>%
    #select only one pair of yes-no for each unique no
    slice_sample(n = 1) %>%
    ungroup() %>%
    mutate(pairID = 1:n())
  
  #get yes dataset row numbers to just remove, since we already
  #know they're too close to the pairs we have
  yes_removes <- map_df(nn$nn, ~as.data.frame(.x), .id="id") %>%
    rename(yes = id,
           no = `.x`) %>%
    anti_join(pairIDs)
  
  #get nos that are not paired
  no_onlys <- nodf2 %>%
    filter(!row_number() %in% pairIDs$no) %>%
    mutate(pairID = NA_real_)
  
  #get yeses that are not paired
  yes_onlys <- yesdf2 %>%
    filter(!row_number() %in% pairIDs$yes) %>%
    filter(!row_number() %in% yes_removes$yes) %>%
    mutate(pairID = NA_real_)
  
  #get those that are in a pair and "yes"
  pair_yes <- yesdf2 %>%
    mutate(yes = as.character(row_number()))%>%
    filter(yes %in% pairIDs$yes) %>%
    #get their pairID column
    left_join(pairIDs, by = 'yes') %>%
    dplyr::select(-yes, -no) 
  
  #get those that are in a pair and "no"
  pair_no <- nodf2 %>%
    mutate(no = row_number()) %>%
    filter(no %in% pairIDs$no) %>%
    #get their pair ID
    left_join(pairIDs, by = "no") %>%
    dplyr::select(-yes, -no) 
  
  #get a dataframe of the pairs together
  pair_df <- pair_yes %>%
    bind_rows(pair_no) %>%
    group_by(pairID) %>%
    #get the buffer to be the maximum for these pairs
    mutate(buffer_m = max(buffer_m) ) %>%
    ungroup()
  
  #set the geometry for each of 
  #the set of a pair as the "centroid" 
  #for those two - so they're sampling from the 
  #same grid cells in the model
  pair_df2 <- pair_df %>%
    group_by(pairID) %>%
    mutate(geometry = st_union(geometry)) %>% 
    mutate(geometry = st_centroid(geometry))
  
  #get only one of the pair set to use in the spatial stratification
  #will have to combine with the rest of the pair ID afterwards
  #get the points for the spatial stratification
  strat <- pair_df2 %>%
    filter(observation_count > 1) %>%
    bind_rows(no_onlys, yes_onlys)
  
  pop <- nrow(strat)
  
  #spatially stratify removing all within 4km of 
  #another one
  strat_df <- subsample.distance(strat, 
                                 d = 4000, 
                                 size = pop-1,
                                 d.max = NULL,
                                 replacement = FALSE)
  
  keep_pairs <- strat_df %>%
    distinct(pairID) %>%
    filter(!is.na(pairID)) %>%
    as_vector()
  
  final_df <- pair_df2 %>%
    filter(observation_count == 0) %>%
    filter(pairID %in% keep_pairs) %>%
    rbind(strat_df)
  
  return(final_df)
  
}

years <- 2010:2024

t <- lapply(years, group_function)

#make into one DF
ebird_spatial3 <- bind_rows(t) %>%
  group_by(year) %>%
  #need to get "group IDs" for 
  #observations that are not in a group
  mutate(val = 1:n(),
         val = val+1000) %>%
  mutate(pairID = case_when(is.na(pairID) ~ val,
                            TRUE ~ pairID)) %>%
  dplyr::select(-val) %>%
  mutate(pairID = as.numeric(as.factor(pairID))) %>%
  ungroup()

# Get polygons around each bird point -------------------------------------

#give all bbs an 8km buffer, 8000m buffer
#give all ebird a buffer equal to sampling distance

#Create buffer zone:
ebird_buffer <- st_buffer(ebird_spatial3, ebird_spatial3$buffer_m)
bbs_buffer <- st_buffer(bbs_spatial2, 8000)#give all bbs an 8km buffer, 8000m buffer
#give all ebird a buffer equal to sampling distance

# Export ------------------------------------------------------------------

st_write(ebird_buffer, here('data',
                             'ebird_data',
                             'cleaned_data',
                             'all_ebird_data_buffercellIDs.shp'))


st_write(bbs_buffer, here('data',
                           'bbs_data',
                           'cleaned_data',
                           'all_bbs_data_buffercellIDs.shp'))

#old code that takes one cell per observation
st_write(ebird_spatial3, here('data',
                      'ebird_data',
                      'cleaned_data',
                      'all_ebird_data_conefiltered.shp'))

st_write(bbs_spatial2, here('data',
                     'bbs_data',
                     'cleaned_data',
                     'all_bbs_data_conefiltered.shp'))


