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

# Load cone dataset -------------------------------------------------------

cones <- terra::rast(here("data",
                          "spatial_data", 
                          "masting_data",
                          "full_pied_masting.tif"))
#plot(cones)

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
# vpdrastlist <- list.files(path = here('data',
#                                        'spatial_data',
#                                        'prism_monthly_vpd'), pattern='.bil$', 
#                            recursive = T, all.files= T, full.names= T)
# 
# vpd_rast <- terra::rast(vpdrastlist)
# 
# vpd_df <- terra::as.data.frame(vpd_rast, 
#                                 xy = TRUE,
#                                 cells = TRUE)

monsoon_rast <- terra::rast(here('data', 'spatial_data',
                                 'monsoon', 'SWMON.tif'))

#terra::plot(monsoon_rast)

monsoon_df <- terra::as.data.frame(monsoon_rast,
                                   xy= TRUE,
                                   cells = TRUE)

pinyonba_rast <- terra::rast(here('data',
                                  'spatial_data',
                                  'pinyonBA',
                                  'PinyonBA_4km_sqftPerAc.tif'))

pinyonba_df <- terra::as.data.frame(pinyonba_rast,
                                    xy = TRUE,
                                    cells = TRUE)

# Set the CRS for datasets -------------------------------------------

#these data are in WGS84, so I'll create this crs to call later
crs_wgs <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

crs_wgs2 <- "+proj=utm +zone=42N +datum=WGS84 +units=km"

crs_albers <- '+proj=longlat +datum=NAD83 +no_defs +type=crs'


# Extract cell IDs from raster ------------------------------------------------------------

##PPT and TEMP
ppt_points <- ppt_df %>%
  dplyr::select(x, y)

ppt_spatial <- st_as_sf(ppt_points, coords = c("x", "y"),
                        crs = st_crs(crs_albers))

ppt_cellIDs <- terra::extract(cones, vect(ppt_spatial), cells = T)

ppt_df$cellID <- ppt_cellIDs$cell

#temp
temp_points <- temp_df %>%
  dplyr::select(x, y)

temp_spatial <- st_as_sf(temp_points, coords = c("x", "y"),
                         crs = st_crs(crs_albers))

temp_cellIDs <- terra::extract(cones, vect(temp_spatial), cells = T)

temp_df$cellID <- temp_cellIDs$cell

#vpd
# vpd_points <- vpd_df %>%
#   dplyr::select(x, y)

# vpd_spatial <- st_as_sf(vpd_points, coords = c("x", "y"),
#                          crs = st_crs(crs_albers))
# 
# vpd_cellIDs <- terra::extract(cones, vect(vpd_spatial), cells = T)
# 
# vpd_df$cellID <- vpd_cellIDs$cell

#monsoon
monsoon_points <- monsoon_df %>%
  dplyr::select(x, y)

monsoon_spatial <- st_as_sf(monsoon_points, coords = c("x", "y"),
                            crs = st_crs(crs_albers))

monsoon_cellIDs <- terra::extract(cones, vect(monsoon_spatial), cells = T)

monsoon_df$cellID <- monsoon_cellIDs$cell

#pinyonba
pinyonba_points <- pinyonba_df %>%
  dplyr::select(x, y)

pinyonba_spatial <- st_as_sf(pinyonba_points, coords = c("x", "y"),
                             crs = st_crs(crs_albers))

pinyonba_cellIDs <- terra::extract(cones, vect(pinyonba_spatial), cells = T)

pinyonba_df$cellID <- pinyonba_cellIDs$cell

# Turn cone raster into df with cell IDs ---------------------------------------

cone_df <- terra::as.data.frame(cones, 
                                xy = TRUE,
                                cells = TRUE)


# Filter datasets for cells with cone data --------------------------------

ids <- cone_df %>%
  distinct(cell)

ppt_df2 <- ppt_df %>%
  filter(cellID %in% ids$cell) %>%
  dplyr::select(-cell)

temp_df2 <- temp_df %>%
  filter(cellID %in% ids$cell) %>%
  dplyr::select(-cell)

# vpd_df2 <- vpd_df %>%
#   filter(cellID %in% ids$cell) %>%
#   dplyr::select(-cell)

monsoon_df2 <- monsoon_df %>%
  filter(cellID %in% ids$cell) %>%
  dplyr::select(-cell)

pinyonba_df2 <- pinyonba_df %>%
  filter(cellID %in% ids$cell) %>%
  dplyr::select(-cell)




# Export ------------------------------------------------------------------

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

