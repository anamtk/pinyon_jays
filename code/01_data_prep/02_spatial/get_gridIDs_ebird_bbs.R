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
plot(cones)

# Load ppt and temp datasets ----------------------------------------------

pptrastlist <- list.files(path = here('data',
                                      'spatial_data',
                                      'prism_monthly_ppt'), pattern='.bil$', 
                          recursive = T, all.files= T, full.names= T)

ppt_rast <- terra::rast(pptrastlist)

ppt_df <- terra::as.data.frame(ppt_rast, 
                               xy = TRUE,
                               cells = TRUE)

#temperature
temprastlist <- list.files(path = here('data',
                                      'spatial_data',
                                      'prism_monthly_tmax'), pattern='.bil$', 
                          recursive = T, all.files= T, full.names= T)

temp_rast <- terra::rast(temprastlist)

temp_df <- terra::as.data.frame(temp_rast, 
                               xy = TRUE,
                               cells = TRUE)

#vpd 
vpdrastlist <- list.files(path = here('data',
                                       'spatial_data',
                                       'prism_monthly_vpd'), pattern='.bil$', 
                           recursive = T, all.files= T, full.names= T)

vpd_rast <- terra::rast(vpdrastlist)

vpd_df <- terra::as.data.frame(vpd_rast, 
                                xy = TRUE,
                                cells = TRUE)

monsoon_rast <- terra::rast(here('data', 'spatial_data',
                                 'monsoon', 'SWMON.tif'))

#terra::plot(monsoon_rast)

monsoon_df <- terra::as.data.frame(monsoon_rast,
                                   xy= TRUE,
                                   cells = TRUE)

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

# Extract cell ID from raster ------------------------------------------------------------

##EBIRD
# Use template raster and ebird pts to get cell IDs
ebird_cellIDs <- terra::extract(cones, vect(ebird_spatial), cells = T)

# Add as column to ebird data
ebird$cellID <- ebird_cellIDs$cell

#some ebird observations don't have masting data cells, so I'll remove those
ebird2 <- ebird %>%
  filter(!is.na(cellID))

##BBS
bbs_cellIDs <- terra::extract(cones, vect(bbs_spatial), cells = T)

bbs2$cellID <- bbs_cellIDs$cell

##PPT and TEMP
ppt_points <- ppt_df %>%
  dplyr::select(x, y)

ppt_spatial <- st_as_sf(ppt_points, coords = c("x", "y"),
                        crs = st_crs(crs_wgs))

ppt_cellIDs <- terra::extract(cones, vect(ppt_spatial), cells = T)

ppt_df$cellID <- ppt_cellIDs$cell

#temp
temp_points <- temp_df %>%
  dplyr::select(x, y)

temp_spatial <- st_as_sf(temp_points, coords = c("x", "y"),
                         crs = st_crs(crs_wgs))

temp_cellIDs <- terra::extract(cones, vect(temp_spatial), cells = T)

temp_df$cellID <- temp_cellIDs$cell

#vpd
vpd_points <- vpd_df %>%
  dplyr::select(x, y)

vpd_spatial <- st_as_sf(vpd_points, coords = c("x", "y"),
                         crs = st_crs(crs_wgs))

vpd_cellIDs <- terra::extract(cones, vect(vpd_spatial), cells = T)

vpd_df$cellID <- vpd_cellIDs$cell

#monsoon
monsoon_points <- monsoon_df %>%
  dplyr::select(x, y)

monsoon_spatial <- st_as_sf(monsoon_points, coords = c("x", "y"),
                        crs = st_crs(crs_wgs))

monsoon_cellIDs <- terra::extract(cones, vect(monsoon_spatial), cells = T)

monsoon_df$cellID <- monsoon_cellIDs$cell

# Turn cone raster into df with cell IDs ---------------------------------------

cone_df <- terra::as.data.frame(cones, 
                                xy = TRUE,
                                cells = TRUE)


# Filter datasets for cells with cone data --------------------------------

ids <- cone_df %>%
  distinct(cell)

ebird3 <- ebird2 %>%
  filter(cellID %in% ids$cell)

bbs3 <- bbs2 %>%
  filter(cellID %in% ids$cell)

ppt_df2 <- ppt_df %>%
  filter(cellID %in% ids$cell) %>%
  dplyr::select(-cell)

temp_df2 <- temp_df %>%
  filter(cellID %in% ids$cell) %>%
  dplyr::select(-cell)

vpd_df2 <- vpd_df %>%
  filter(cellID %in% ids$cell) %>%
  dplyr::select(-cell)

monsoon_df2 <- monsoon_df %>%
  filter(cellID %in% ids$cell) %>%
  dplyr::select(-cell)

# Export ------------------------------------------------------------------

write.csv(ebird3, here('data',
                      'ebird_data',
                      'cleaned_data',
                      'all_ebird_data_cellIDs.csv'))

write.csv(bbs3, here('data',
                     'bbs_data',
                     'cleaned_data',
                     'pjay_data_co_nm_cellIDs.csv'))

write.csv(cone_df, here('data',
                        'spatial_data',
                        'cleaned_data',
                        'cone_masting_df.csv'))

write.csv(temp_df2, here('data',
                         'spatial_data',
                         'cleaned_data',
                         'temp_data_df.csv'))

write.csv(ppt_df2, here('data',
                         'spatial_data',
                         'cleaned_data',
                         'ppt_data_df.csv'))

write.csv(vpd_df2, here('data',
                        'spatial_data',
                        'cleaned_data',
                        'vpd_data_df.csv'))

write.csv(monsoon_df2, here('data',
                        'spatial_data',
                        'cleaned_data',
                        'monsoon_data_df.csv'))

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

#there are 17462 unique 4km cells
ebird2 %>%
  distinct(cellID) %>%
  tally()

#there are 1-7315 observations in each cell - will probs want
#to subsample in these more populated cells
ebird2 %>%
  group_by(cellID) %>%
  tally() %>%
  arrange(desc(n))

#1-945 in a single cell in a given year - will want to subset these
ebird2 %>%
  group_by(year, cellID) %>%
  tally() %>%
  arrange(desc(n))

#388 unique 4km cells
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

ebird2 %>%
  summarise(min = min(cellID),
            max = max(cellID))

