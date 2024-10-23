
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
                           'all_ebird_data_buffercellIDs.shp'))


# Get eBIRD grid IDs ------------------------------------------------------
ba_cells <- pinyonba_df %>%
  distinct(cell)

ebird_cells <- exact_extract(pinyonba_rast,
                           ebird_buffer,
                           include_cell = T, 
                           include_cols = c("chckls_", "year", 
                                            "buffr_m"),
                           force_df = T)


ebird_df <- bind_rows(ebird_cells)

ebird_df2 <- ebird_df %>%
  dplyr::select(chckls_, year,
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
  dplyr::select(chckls_, year, cell, prop) %>%
  filter(cell %in% ba_cells$cell)


# Get cell IDs to be able to filter covariates ----------------------------

cells <- ebird_df2 %>%
  distinct(cell, year)
  
# Explort -----------------------------------------------------------------

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







