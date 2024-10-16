
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


# Load CRS ----------------------------------------------------------------

#crs_albers <- '+proj=longlat +datum=NAD83 +no_defs +type=crs'

# Load cone dataset -------------------------------------------------------

cones <- terra::rast(here("data",
                          "spatial_data", 
                          "masting_data",
                          "quantile_pied_predictions_scaled.tif"))

names(cones) <- 1897:2024 ## Actual years of maturation are one year earlier, but adding one year to each to associate it with jay survey period



bbs_buffer <- read_sf(here('data',
                          'bbs_data',
                          'cleaned_data',
                          'all_bbs_data_buffercellIDs.shp'))


ebird_buffer <- read_sf(here('data',
                           'ebird_data',
                           'cleaned_data',
                           'all_ebird_data_buffercellIDs.shp'))



# Turn cones to df --------------------------------------------------------

cone_df <- terra::as.data.frame(cones, 
                                xy = TRUE,
                                cells = TRUE)

# Get grid IDs for BBS ----------------------------------------------------
#the buffer I think gets cells outside what is
#in the cones dataframe (which is what we're masking 
#every dataset to)
#so we can delete those cells outside the range
#for now

cone_cells <- cone_df %>%
  distinct(cell)

bbs_cells <- exact_extract(cones,
                       bbs_buffer,
                      include_cell = T, 
                      include_cols = c("StateNm", "Route", 
                                       "Year", "RPID"),
                      force_df = T)


bbs_df <- bind_rows(bbs_cells)

bbs_df2 <- bbs_df %>%
  dplyr::select(StateNm, Route, Year,
                RPID, cell, coverage_fraction) %>%
  rowwise() %>%
  #get the total area of a 4km x 4km grid cell
  #and the total amount of that that is covered
  #by the polygon
  mutate(area = 1.6e+07*coverage_fraction) %>%
  ungroup() %>%
  group_by(StateNm, Route, Year,
           RPID) %>%
  #get the total area of each sampling 
  #buffer
  mutate(total_area = sum(area)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(prop = area/total_area) %>%
  ungroup() %>%
  dplyr::select(StateNm, Route, Year, 
                RPID, cell, prop) %>%
  filter(cell %in% cone_cells$cell)


# Get eBIRD grid IDs ------------------------------------------------------


ebird_cells <- exact_extract(cones,
                           ebird_buffer,
                           include_cell = T, 
                           include_cols = c("chckls_", "year", 
                                            'pairID',
                                            "buffr_m"),
                           force_df = T)


ebird_df <- bind_rows(ebird_cells)

ebird_df2 <- ebird_df %>%
  dplyr::select(chckls_, year, pairID,
                buffr_m, cell, coverage_fraction) %>%
  rowwise() %>%
  mutate(area = 1.6e+07*coverage_fraction) %>%
  ungroup() %>%
  group_by(chckls_, year) %>%
  mutate(total_area = sum(area)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(prop = area/total_area) %>%
  ungroup() %>%
  dplyr::select(chckls_, pairID, year, cell, prop) %>%
  filter(cell %in% cone_cells$cell)


# Get cell IDs to be able to filter covariates ----------------------------

cell1 <- bbs_df2 %>%
  distinct(cell, Year) %>%
  rename(year = Year)

cell2 <- ebird_df2 %>%
  distinct(cell, year)
  
cells <- cell1 %>%
  bind_rows(cell2) %>%
  distinct(cell, year) 

# Explort -----------------------------------------------------------------

write.csv(bbs_df2,
          here('data',
               'bbs_data',
               'cleaned_data',
               'bbs_cellIDlists.csv'))

write.csv(ebird_df2,
          here('data',
               'ebird_data',
               'cleaned_data',
               'ebird_cellIDlists.csv'))


write.csv(cells, 
          here('data',
               'spatial_data',
               'cleaned_data',
               'cellIDs.csv'))







