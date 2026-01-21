
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


# Load cone dataset -------------------------------------------------------

#for masking to
pinyonba_rast <- terra::rast(here('data',
                                  'spatial_data',
                                  'pinyonBA',
                                  'PinyonBA_4km_sqmPerHa.tif'))

pinyonba_df <- terra::as.data.frame(pinyonba_rast,
                                    xy = TRUE,
                                    cells = TRUE)

ebird_buffer <- read_sf(here('data',
                           'ebird_data',
                           'cleaned_data',
                           "05_oos",
                           'all_oos_ebird_data_buffercellIDs.shp')) %>%
  st_make_valid() %>% 
  rowwise() %>%
  mutate(area = st_area(geometry)) %>%
  ungroup()


# Get eBIRD grid IDs ------------------------------------------------------
ba_cells <- pinyonba_df %>%
  distinct(cell)

ebird_cells <- exact_extract(pinyonba_rast,
                           ebird_buffer,
                           include_cell = T, 
                           include_cols = c("cellID", "year", "area"),
                           force_df = T)

ebird_df <- bind_rows(ebird_cells)

ebird_df2 <- ebird_df %>%
  dplyr::select(cellID, year, cell,area, coverage_fraction) %>%
  rowwise() %>%
  ungroup() %>%
  dplyr::select(cellID, year, cell,area, coverage_fraction) %>%
  filter(cell %in% ba_cells$cell) %>%
  group_by(cellID, year) %>%
  mutate(blobnum = cur_group_id()) %>%
  ungroup() 


# Get cell IDs to be able to filter covariates ----------------------------

cells <- ebird_df2 %>%
  distinct(cell, year)
  
# Explort -----------------------------------------------------------------

write.csv(ebird_df2,
          here('data',
               'ebird_data',
               'cleaned_data',
               '05_oos',
               'ebird_oos_cellIDlists.csv'))


write.csv(cells, 
          here('data',
               'spatial_data',
               'cleaned_data',
               '05_oos',
               'oos_cellIDs.csv'))

