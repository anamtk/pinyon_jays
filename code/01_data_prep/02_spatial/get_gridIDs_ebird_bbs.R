#Ana Miller-ter Kuile
#make ebird and bbs raster data
#July 2, 2024

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "sf",  
                  "terra",
                  'readxl',
                  'sf',
                  'exactextractr')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

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
ebird2$cellID <- ebird_cellIDs$cell

#some ebird observations don't have masting data cells, so I'll remove those
ebird3 <- ebird2 %>%
  filter(!is.na(cellID))

##BBS
bbs_cellIDs <- terra::extract(cones, vect(bbs_spatial), cells = T)

bbs2$cellID <- bbs_cellIDs$cell

# Turn cone raster into df with cell IDs ---------------------------------------

cone_df <- terra::as.data.frame(cones, 
                                xy = TRUE,
                                cells = TRUE)

# Filter datasets for cells with cone data --------------------------------

ids <- cone_df %>%
  distinct(cell)

ebird4 <- ebird3 %>%
  filter(cellID %in% ids$cell)

bbs3 <- bbs2 %>%
  filter(cellID %in% ids$cell)

# Subset ebird data -------------------------------------------------------
#this paper: 
#https://www.nature.com/articles/s41598-022-23603-0
#subset one 1 and one zero dataset from each cell
#the chapter from cornell lab subsets "10 per site", not specifying whether
#they're 1 or 0 checklists
#https://cornelllabofornithology.github.io/ebird-best-practices/occupancy.html#occupancy-intro

#going the cornell lab route for now - can come back to revisit later
#the method from the paper
ebird_yes <- ebird4 %>%
  filter(observation_count > 0) %>%
  group_by(year, cellID) %>%
  slice_sample(n = 1) %>%
  ungroup()

ebird_no <- ebird4 %>%
  filter(observation_count == 0) %>%
  group_by(year, cellID) %>%
  slice_sample(n = 1) %>%
  ungroup()

ebird5 <- ebird_yes %>%
  rbind(ebird_no)


# Get polygons around each bird point -------------------------------------

#convert all ebird data to an sf object
ebird_spatial2 <- ebird_spatial %>%
  filter(checklist_id %in% ebird5$checklist_id)

#convert all bbs data to an sf object
bbs_spatial2 <- bbs_spatial %>%
  filter(RouteDataID %in% bbs3$RouteDataID)

#give all bbs an 8km buffer, 8000m buffer
#give all ebird a buffer equal to sampling distance

#Create buffer zone:
ebird_buffer <- st_buffer(ebird_spatial2, ebird_spatial2$buffer_m)
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
# write.csv(ebird5, here('data',
#                       'ebird_data',
#                       'cleaned_data',
#                       'all_ebird_data_cellIDs.csv'))
# 
# write.csv(bbs3, here('data',
#                      'bbs_data',
#                      'cleaned_data',
#                      'pjay_data_co_nm_cellIDs.csv'))



# Extract 30m cells under 1km grid and find vlaues and weights


