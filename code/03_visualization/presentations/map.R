

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "sf", 
                  'gpkg',
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
                       '01_ebird_data',
                       'cleaned_data',
                       'all_ebird_data_conefiltered.shp'))

pinyonba_rast <- terra::rast(here('data',
                                  '02_spatial_data',
                                  'pinyonBA',
                                  'PinyonBA_4km_sqmPerHa.tif'))

pinyonba_df <- terra::as.data.frame(pinyonba_rast,
                                    xy = TRUE,
                                    cells = TRUE)

cones <- terra::rast(here("data",
                          "02_spatial_data", 
                          "masting_data",
                          "ConePredictions_final.tif"))

cone_df <- terra::as.data.frame(cones, 
                                xy = TRUE,
                                cells = TRUE)

range <- st_read(here('data',
                      '01_ebird_data',
                      'pinjay_range_2023',
                      'pinjay_range_2023.gpkg'))


range2 <- range %>%
  st_transform(st_crs(pinyonba_rast))


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

sw2 <- sw %>%
  st_transform(st_crs(pinyonba_rast))

sf_use_s2(TRUE)


# Plot?? ------------------------------------------------------------------

ebird11 <- ebird %>%
  filter(year == 2011)

ba10 <- pinyonba_df %>%
  mutate(PinyonBA_sqftPerAc_2010 = case_when(PinyonBA_sqftPerAc_2010 == 0 ~ NA_real_,
                                             TRUE ~ PinyonBA_sqftPerAc_2010)) 

cone11 <- cone_df %>%
  mutate(`2010` = case_when(`2011` == 0 ~ NA_real_,
                            TRUE ~ `2011`))

theme_set(theme_void())

ggplot()+
  geom_sf(data = range2, fill = "purple", alpha = 0.2) +
  geom_tile(data = cone11, 
            aes(x = x, 
                y = y, 
                fill = `2011`)) +
    geom_sf(data = states, fill = "grey", alpha =.2) +
  scale_fill_viridis_c(na.value = 'transparent',
                       limits = c(0.01, 0.65), 
                       breaks = c(0.01,0.65),
                       labels = c('low', 'high')) +
  geom_sf(data = ebird11, shape = 1, alpha = 0.4) +
    labs(fill = "Cone availability")

ggplot()+
  geom_sf(data = range2, fill = "purple", alpha = 0.2) +
  # geom_tile(data = pinyonba_df, 
  #           aes(x = x, 
  #               y = y, 
  #               fill = `PinyonBA_sqftPerAc_2010`)) +
  geom_sf(data = states, fill = "grey", alpha =.2) +
  scale_fill_viridis_c(na.value = 'transparent',
                       limits = c(0.01, 0.65), 
                       breaks = c(0.01,0.65),
                       labels = c('low', 'high')) +
  geom_sf(data = ebird11, shape = 1, alpha = 0.4) +
  labs(fill = "Pinyon basal area")
  
ggplot() +
    geom_sf(data = range2, fill = "purple", alpha = 0.2) +
    geom_sf(data = states2, alpha = 0.2)
  
intersect_pct <- st_intersection(range2, sw2) 

ggplot() +
  geom_sf(data = intersect_pct)


st_area(intersect_pct)/st_area(range2)
