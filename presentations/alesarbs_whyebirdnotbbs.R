
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

ebird <- read_sf(here('data',
                      'ebird_data',
                      'cleaned_data',
                      'all_ebird_data_conefiltered.shp'))

pinyonba_rast <- terra::rast(here('data',
                                  'spatial_data',
                                  'pinyonBA',
                                  'PinyonBA_4km_sqmPerHa.tif'))

pinyonba_df <- terra::as.data.frame(pinyonba_rast,
                                    xy = TRUE,
                                    cells = TRUE)

cones <- terra::rast(here("data",
                          "spatial_data", 
                          "masting_data",
                          "ConePredictions_final.tif"))

cone_df <- terra::as.data.frame(cones, 
                                xy = TRUE,
                                cells = TRUE)

#this gives an error if default TRUE
sf_use_s2(FALSE)

#get the states we care aboute
states <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE)) %>%
  filter(ID %in% c('arizona', 'colorado', 
                   'utah', 'new mexico')) 

#make them "valid" (e.g. not overlapping)
states <- st_make_valid(states)

#get the total boundary geometry for these states,
#rather than a boundary for each state
sw <- states %>%
  summarise(geometry = sf::st_union(geom))
#switch back - since everything else works
sf_use_s2(TRUE)

sw2 <- sw %>%
  st_transform(st_crs(pinyonba_rast))

arizona <- states %>%
  filter(ID == "arizona")

bbs <- read.csv(here('data',
                         'bbs_data',
                         'cleaned_data',
                         'pjay_data_co_nm_az_ut_wlocations.csv')) %>%
  dplyr::select(RouteDataID:AOU,ObserverExp:Longitude) %>%
  #filter only AZ
  filter(StateNum == 6) 

# Plot ------------------------------------------------------------------

# Plot of basal area in one year -----------------------------------

pinyonba_rast2 <- mask(pinyonba_rast, sw2)

pinyonba_df2 <- terra::as.data.frame(pinyonba_rast2,
                                          xy = TRUE,
                                          cells = TRUE)

ggplot()+
  geom_tile(data = pinyonba_df2, 
            aes(x = x, 
                y = y, 
                fill = `PinyonBA_sqftPerAc_2011`)) +
  geom_sf(data = states, fill = "white", alpha = 0.1, linewidth = 0.5) +
  geom_sf(data = sw2, fill = "white", alpha = 0.2, linewidth = 1,
          color = "black") +
  scale_fill_viridis_c(na.value = 'white',
                       limits = c(0.01, 0.65), 
                       breaks = c(0.01,0.65),
                       labels = c('low', 'high')) +
  labs(fill='Pinyon basal area') 

ggsave(here('pictures',
            'presentations',
            'ales_arb_pinyonba.jpg'),
       height = 5,
       width = 5,
       units = c("in"))
  
# Plot of data availability for ebird bbs ---------------------------------

azebird11 <- ebird %>%
  filter(year == 2011) %>%
  filter(stat_cd == 'US-AZ')

bbs11 <- bbs %>%
  filter(Year == 2011)

bbs_spatial <- st_as_sf(bbs11, coords = c("Longitude", "Latitude"),
                          crs = st_crs(cones))

bbs_cellIDs <- terra::extract(pinyonba_rast2, vect(bbs_spatial), cells = T)

cells_pres <- bbs_cellIDs %>%
  filter_at(vars(-ID, -cell), any_vars(. != 0)) 

cells_abs <- bbs_cellIDs %>%
  anti_join(cells_pres) %>%
  mutate(cell = NA_real_)

cells <- cells_pres %>%
  bind_rows(cells_abs) %>%
  arrange(ID)

# Add as column to ebird data
bbs_spatial$cellID <- cells$cell

#some ebird observations don't have masting data cells, so I'll remove those
bbs_spatial <- bbs_spatial %>%
  filter(!is.na(cellID)) 

arizona2 <- arizona %>%
  st_transform(st_crs(cones))

azba <- mask(pinyonba_rast2, arizona2)

azba_df <- terra::as.data.frame(azba,
                                    xy = TRUE,
                                    cells = TRUE)

theme_set(theme_void())

ggplot()+
  geom_tile(data = azba_df, 
            aes(x = x, 
                y = y, 
                fill = `PinyonBA_sqftPerAc_2011`)) +
  geom_sf(data = arizona2, 
          fill = "grey85", 
          alpha =.2,
          linewidth = 1,
          color= 'black') +
  scale_fill_viridis_c(na.value = 'white',
                       limits = c(0.01, 0.65), 
                       breaks = c(0.01,0.65),
                       labels = c('low', 'high')) +
  labs(fill = "Pinyon basal area")

ggsave(here('pictures',
            'presentations',
            'ales_arb_azpinyonba.jpg'),
       height = 5,
       width = 5,
       units = c("in"))

ggplot()+
  geom_tile(data = azba_df, 
            aes(x = x, 
                y = y, 
                fill = `PinyonBA_sqftPerAc_2011`)) +
  geom_sf(data = arizona2, fill = "grey85", alpha =.2,
          linewidth = 1,
          color= 'black') +
  scale_fill_viridis_c(na.value = 'white',
                       limits = c(0.01, 0.65), 
                       breaks = c(0.01,0.65),
                       labels = c('low', 'high')) +
  geom_sf(data = bbs_spatial, 
          shape = 21, 
          color = "black", 
          fill = "black",
          size = 6) +
  geom_sf(data = bbs_spatial, 
          shape = 21, 
          color = "black", 
          fill = "#737373",
          size = 4) +
  labs(fill = "Pinyon basal area")

ggsave(here('pictures',
            'presentations',
            'ales_arb_azpinyonba_bbs.jpg'),
       height = 5,
       width = 5,
       units = c("in"))

ggplot()+
  geom_tile(data = azba_df, 
            aes(x = x, 
                y = y, 
                fill = `PinyonBA_sqftPerAc_2011`)) +
  geom_sf(data = arizona2, fill = "grey85", alpha =.2,
          linewidth = 1,
          color= 'black') +
  scale_fill_viridis_c(na.value = 'white',
                       limits = c(0.01, 0.65), 
                       breaks = c(0.01,0.65),
                       labels = c('low', 'high')) +
  geom_sf(data = bbs_spatial, 
          shape = 21, 
          color = "black", 
          fill = "black",
          size = 6) +
  geom_sf(data = bbs_spatial, 
          shape = 21, 
          color = "black", 
          fill = "#737373",
          size = 4) +
  geom_sf(data = azebird11,
          shape = 21, 
          color = "black", 
          fill = "black",
          size = 6) +
  geom_sf(data = azebird11,
          shape = 21, 
          color = "black", 
          fill = "#252525",
          size = 4) +
  labs(fill = "Pinyon basal area")

ggsave(here('pictures',
            'presentations',
            'ales_arb_azpinyonba_bbsebird.jpg'),
       height = 5,
       width = 5,
       units = c("in"))
