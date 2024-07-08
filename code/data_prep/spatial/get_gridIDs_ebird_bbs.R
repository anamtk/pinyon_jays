#Ana Miller-ter Kuile
#make ebird and bbs raster data
#July 2, 2024

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "sf",  
                  "terra",
                  'readxl')

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
                     'pjay_data_co_nm.csv'))

bbs_locations <- read.csv(here('data',
                               'bbs_data',
                               '2023Release_Nor',
                               'routes.csv'))

bbs_locations <- bbs_locations %>%
  dplyr::select(StateNum, Route, Latitude, Longitude)

bbs2 <- bbs %>%
  left_join(bbs_locations, by = c("StateNum", "Route"))


# Load cone dataset -------------------------------------------------------

cones <- read_xlsx(here('data',
                        'spatial_data',
                        'cones.xlsx'))

cones <- cones %>%
  mutate(Longitude = -Longitude)
# Read in example grid ----------------------------------------------------

#may need to download different data, we got this here:
#https://prism.oregonstate.edu/recent/
grid <- rast(x = here('data',
                      'spatial_data',
                      'example_grid',
                      'PRISM_ppt_stable_4kmM3_2023_bil.bil'))

#look at it
plot(grid)

# Set the CRS for both bird datasets -------------------------------------------

#these data are in WGS84, so I'll create this crs to call later
crs_wgs <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

# Make spatial ------------------------------------------------------------

#convert all ebird data to an sf object
ebird_spatial <- st_as_sf(ebird, coords = c("longitude", "latitude"),
                          crs = st_crs(crs_wgs))

#convert all bbs data to an sf object
bbs_spatial <- st_as_sf(bbs2, coords = c("Longitude", "Latitude"),
                        crs = st_crs(crs_wgs))

cones_spatial <- st_as_sf(cones, coords = c("Longitude", "Latitude"),
                          crs = st_crs(crs_wgs))

# Extract cell ID from raster ------------------------------------------------------------

# Use template raster and ebird pts to get cell IDs
ebird_cellIDs <- terra::extract(grid, vect(ebird_spatial), cells = T)

# Add as column to ebird data
ebird$cellID <- ebird_cellIDs$cell

#BBS
bbs_cellIDs <- terra::extract(grid, vect(bbs_spatial), cells = T)

bbs2$cellID <- bbs_cellIDs$cell

#cones
cone_cellIDs <- terra::extract(grid, vect(cones_spatial), cells = T)

cones$cellID <- cone_cellIDs$cell
# Export ------------------------------------------------------------------

write.csv(ebird, here('data',
                      'ebird_data',
                      'cleaned_data',
                      'all_ebird_data_cellIDs.csv'))

write.csv(bbs2, here('data',
                     'bbs_data',
                     'cleaned_data',
                     'pjay_data_co_nm_cellIDs.csv'))

write.csv(cones, here('data',
                     'spatial_data',
                     'cleaned_data',
                     'cones_cellIDs.csv'))


#Get the unique cell IDs that have data - we will just
#model these cells since it is MUCH smaller than the total 
#dataset
# 
# ebird_cells <- ebird %>%
#   distinct(cellID)
# 
# bbs_cells <- bbs2 %>%
#   distinct(cellID)
# 
# #all the cells in the dataset:
# all_cells <- ebird_cells %>%
#   full_join(bbs_cells, by = "cellID") 
# #10342 - still a lot but much smaller than 300,000+!
# 
# write.csv(all_cells, here('data',
#                           'spatial_data',
#                           'cleaned_data',
#                           'all_cellIDs.csv'))

# Look at stats for bbs and ebird datasets --------------------------------

#there are 10275 unique 4km cells
ebird %>%
  distinct(cellID) %>%
  tally()

#there are 1-3135 observations in each cell - will probs want
#to subsample in these more populated cells
ebird %>%
  group_by(cellID) %>%
  tally() %>%
  arrange(n)

#1-599 in a single cell in a given year - will want to subset these
ebird %>%
  group_by(year, cellID) %>%
  tally() %>%
  arrange(desc(n))

#216 unique 4km cells
bbs2 %>%
  distinct(cellID) %>%
  tally()

#to get cell number for next model:
#total min: 301078
#total max: 627033
#325955 total cells - EEK, this might take a LOOOONG time to model
#if we use them all...

bbs2 %>%
  summarise(min = min(cellID),
            max = max(cellID))

ebird %>%
  summarise(min = min(cellID),
            max = max(cellID))

