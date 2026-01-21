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
                          "02_spatial_data", 
                          "masting_data",
                          "ConePredictions_final.tif"))

cone_df <- terra::as.data.frame(cones, 
                                xy = TRUE,
                                cells = TRUE)

# Load ppt and temp datasets ----------------------------------------------

pptrastlist <- list.files(path = here('data',
                                      '02_spatial_data',
                                      'prism_monthly_ppt'), pattern='.bil$', 
                          recursive = T, all.files= T, full.names= T)

ppt_rast <- terra::rast(pptrastlist)


ppt_rast2 <- project(ppt_rast,
                     crs(cones))

all.equal(crs(cones), crs(ppt_rast2)) ## The same CRS for both layers

ppt_df <- terra::as.data.frame(ppt_rast2, 
                               xy = TRUE,
                               cells = TRUE)

#temperature
temprastlist <- list.files(path = here('data',
                                       '02_spatial_data',
                                       'prism_monthly_tmax'), pattern='.bil$', 
                           recursive = T, all.files= T, full.names= T)

temp_rast <- terra::rast(temprastlist)

temp_rast2 <- project(temp_rast,
                     crs(cones))

temp_df <- terra::as.data.frame(temp_rast2, 
                                xy = TRUE,
                                cells = TRUE)

#tmean
tmeanrastlist <- list.files(path = here('data',
                                       '02_spatial_data',
                                       'prism_monthly_tmean'), pattern='.bil$', 
                           recursive = T, all.files= T, full.names= T)

tmean_rast <- terra::rast(tmeanrastlist)

tmean_rast2 <- project(tmean_rast,
                      crs(cones))

tmean_df <- terra::as.data.frame(tmean_rast2, 
                                xy = TRUE,
                                cells = TRUE)


monsoon_rast <- terra::rast(here('data', '02_spatial_data',
                                 'monsoon', 'SWMON.tif'))

#terra::plot(monsoon_rast)
monsoon_rast2 <- project(monsoon_rast,
                     crs(cones))

monsoon_df <- terra::as.data.frame(monsoon_rast2,
                                   xy= TRUE,
                                   cells = TRUE)

pinyonba_rast <- terra::rast(here('data',
                                  '02_spatial_data',
                                  'pinyonBA',
                                  'PinyonBA_4km_sqmPerHa.tif'))

all.equal(crs(cones), crs(pinyonba_rast)) ## The same CRS for both layers

pinyonba_df <- terra::as.data.frame(pinyonba_rast,
                                    xy = TRUE,
                                    cells = TRUE)

# Extract cell IDs from raster ------------------------------------------------------------

##PPT and TEMP
ppt_points <- ppt_df %>%
  dplyr::select(x, y)

ppt_spatial <- st_as_sf(ppt_points, coords = c("x", "y"),
                        crs = st_crs(cones))

ppt_cellIDs <- terra::extract(pinyonba_rast, vect(ppt_spatial), cells = T)

ppt_df$cellID <- ppt_cellIDs$cell

#temp
temp_points <- temp_df %>%
  dplyr::select(x, y)

temp_spatial <- st_as_sf(temp_points, coords = c("x", "y"),
                         crs = st_crs(cones))

temp_cellIDs <- terra::extract(pinyonba_rast, vect(temp_spatial), cells = T)

temp_df$cellID <- temp_cellIDs$cell

#tmean
tmean_points <- tmean_df %>%
  dplyr::select(x, y)

tmean_spatial <- st_as_sf(tmean_points, coords = c("x", "y"),
                         crs = st_crs(cones))

tmean_cellIDs <- terra::extract(pinyonba_rast, vect(tmean_spatial), cells = T)

tmean_df$cellID <- tmean_cellIDs$cell

#monsoon
monsoon_points <- monsoon_df %>%
  dplyr::select(x, y)

monsoon_spatial <- st_as_sf(monsoon_points, coords = c("x", "y"),
                            crs = st_crs(cones))

monsoon_cellIDs <- terra::extract(pinyonba_rast, vect(monsoon_spatial), cells = T)

monsoon_df$cellID <- monsoon_cellIDs$cell

#pinyonba
pinyonba_points <- pinyonba_df %>%
  dplyr::select(x, y)

pinyonba_spatial <- st_as_sf(pinyonba_points, coords = c("x", "y"),
                             crs = st_crs(cones))

# pinyonba_cellIDs <- terra::extract(cones, vect(pinyonba_spatial), cells = T)
# 
# pinyonba_df$cellID <- pinyonba_cellIDs$cell

# Filter datasets for cells with cone data --------------------------------

ids <- pinyonba_df %>%
  distinct(cell)
# 
ppt_df2 <- ppt_df %>%
  filter(cellID %in% ids$cell) %>%
  dplyr::select(-cell)
# 
temp_df2 <- temp_df %>%
  filter(cellID %in% ids$cell) %>%
  dplyr::select(-cell)

tmean_df2 <- tmean_df %>%
  filter(cellID %in% ids$cell) %>%
  dplyr::select(-cell)

monsoon_df2 <- monsoon_df %>%
  filter(cellID %in% ids$cell) %>%
  dplyr::select(-cell) %>%
  group_by(cellID) %>%
  summarise(PRISM_ppt_30yr_normal_800mM4_07_bil = mean(PRISM_ppt_30yr_normal_800mM4_07_bil, na.rm = T)) %>%
  ungroup()

# Export ------------------------------------------------------------------

write.csv(cone_df, here('data',
                        'spatial_data',
                        'cleaned_data',
                        '01_get_gridIDs',
                        'cone_masting_df.csv'))


write.csv(temp_df2, here('data',
                         'spatial_data',
                         'cleaned_data',
                         '01_get_gridIDs',
                         'temp_data_df.csv'))

write.csv(tmean_df2, here('data',
                         'spatial_data',
                         'cleaned_data',
                         '01_get_gridIDs',
                         'tmean_data_df.csv'))

write.csv(ppt_df2, here('data',
                        'spatial_data',
                        'cleaned_data',
                        '01_get_gridIDs',
                        'ppt_data_df.csv'))


write.csv(monsoon_df2, here('data',
                            'spatial_data',
                            'cleaned_data',
                            '01_get_gridIDs',
                            'monsoon_data_df.csv'))

write.csv(pinyonba_df, here('data',
                             'spatial_data',
                             'cleaned_data',
                            '01_get_gridIDs',
                             'pinyonba_data_df.csv'))

