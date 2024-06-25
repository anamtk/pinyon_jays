
# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "auk", 'sf', 'mapdata',
                  'usmap')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

ebirdco <- read_delim(here('data',
                         'ebird_data',
                         'ebd_US-CO_pinjay_smp_relMay-2024',
                         'ebd_US-CO_pinjay_smp_relMay-2024.txt'))

checkco <- read_delim(here('data',
                           'ebird_data',
                           'ebd_US-CO_pinjay_smp_relMay-2024',
                           'ebd_US-CO_pinjay_smp_relMay-2024_sampling.txt'))

ebirdnm <- read_delim(here('data',
                           'ebird_data',
                           'ebd_US-NM_pinjay_smp_relMay-2024',
                           'ebd_US-NM_pinjay_smp_relMay-2024.txt'))

checknm <- read_delim(here('data',
                           'ebird_data',
                           'ebd_US-NM_pinjay_smp_relMay-2024',
                           'ebd_US-NM_pinjay_smp_relMay-2024_sampling.txt'))


points <- read_sf(here('data',
                       'spatial_data',
                       'pts4ana',
                       'pts4ana.shp'))

state <- map_data('state')

sw <- state %>%
  filter(region %in% c('colorado',
                       'new mexico'))

usa <- map_data('usa')
canada <- map_data('worldHires', 'Canada')
mexico <- map_data('worldHires', 'Mexico')


# Where are plots ---------------------------------------------------------

ggplot() +
  geom_polygon(data = sw,
               aes(x = long, y = lat, group = group), 
               fill = "white", color = "black") +
  geom_sf(data = points) 

# Filter four corners area ------------------------------------------------

ebird_small <- ebird %>%
  filter(STATE %in% c("Colorado", "New Mexico")) %>%
  mutate(year = str_sub(`LAST EDITED DATE`, 1, 4))

obs <- ebird_small %>%
  group_by(year) %>%
  tally()

ebird_small <- ebird_small %>%
  left_join(obs, by = "year") %>%
  filter(year > 2005)

# Plot --------------------------------------------------------------------

ebird_sf <- st_as_sf(ebird_small, coords = c("LONGITUDE", "LATITUDE"))

ggplot() +
  geom_polygon(data = sw,
               aes(x = long, y = lat, group = group), 
               fill = "white", color = "black") +
  #geom_sf(data = ebird_sf) +
  stat_density2d(data = ebird_small, aes(x = LONGITUDE,
                                         y = LATITUDE,
                                         fill = after_stat(density)), geom = 'tile',
                 contour = F, alpha = 0.5, bins = 100) +
  scale_fill_gradient2(low = "white", high = 'red') +
  geom_sf(data = points, alpha = 0.6) 
  
ggplot() +
  geom_polygon(data = sw,
               aes(x = long, y = lat, group = group), 
               fill = "white", color = "black") +
  geom_point(data = ebird_small, aes(x = LONGITUDE,
                                  y = LATITUDE), 
             color = '#1c9099',
             alpha = 0.6) +
  facet_wrap(~ year, nrow = 3) +
  geom_sf(data = points, shape = 1) 



