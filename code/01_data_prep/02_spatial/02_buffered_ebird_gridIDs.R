
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

sf_use_s2(FALSE)
ebird_buffer$area <- st_area(ebird_buffer)
sf_use_s2(TRUE)

# Get eBIRD grid IDs ------------------------------------------------------
ba_cells <- pinyonba_df %>%
  distinct(cell)

ebird_cells <- exact_extract(pinyonba_rast,
                           ebird_buffer,
                           include_cell = T, 
                           include_cols = c("cellID", "year"),
                           force_df = T)

area <- as.data.frame(ebird_buffer) %>%
  dplyr::select(cellID, year, area)

ebird_df <- bind_rows(ebird_cells) %>%
  left_join(area, by = c("cellID", "year"))

ebird_df2 <- ebird_df %>%
  dplyr::select(cellID, year, cell, coverage_fraction, area) %>%
  rowwise() %>%
  mutate(cell_area = area*coverage_fraction) %>%
  mutate(prop = cell_area/area) %>%
  ungroup() %>%
  dplyr::select(cellID, year, cell, prop) %>%
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







