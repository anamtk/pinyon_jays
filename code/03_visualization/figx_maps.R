

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "sf", 
                  'gpkg',
                  "terra",
                  'readxl',
                  'sf',
                  'exactextractr',
                  'hexbin',
                  'ggtext',
                  'patchwork')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_void())

# Load data ---------------------------------------------------------------

ebird <- read_sf(here('data',
                       '01_ebird_data',
                       'cleaned_data',
                      '03_subsampled',
                       'all_ebird_data_conefiltered.shp'))

pinyonba_rast <- terra::rast(here('data',
                                  '02_spatial_data',
                                  'pinyonBA',
                                  'PinyonBA_4km_sqmPerHa.tif'))

pinyonba_df <- terra::as.data.frame(pinyonba_rast,
                                    xy = TRUE,
                                    cells = TRUE)

# cones <- terra::rast(here("data",
#                           "02_spatial_data", 
#                           "masting_data",
#                           "ConePredictions_final.tif"))
# 
# cone_df <- terra::as.data.frame(cones, 
#                                 xy = TRUE,
#                                 cells = TRUE)

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

#get all the states
states2 <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE)) 

#get the total boundary geometry for these states,
#rather than a boundary for each state
sw <- states %>%
  summarise(geometry = sf::st_union(geom))

sw2 <- sw %>%
  st_transform(st_crs(pinyonba_rast))

canada <- st_as_sf(maps::map('world', 
                             region = c("Canada"),
                             fill=TRUE, 
                             plot =FALSE),
                   coords = c("long", "lat"))

mexico <- st_as_sf(maps::map('world', 
                            region = c("Mexico"),
                            fill=TRUE, 
                            plot =FALSE),
                   coords = c("long", "lat"))

mexico2 <- mexico %>%
  st_transform(st_crs(pinyonba_rast))

canada2 <- canada %>%
  st_transform(st_crs(pinyonba_rast))

pinyonba_rast2 <- terra::mask(pinyonba_rast, sw2)

pinyonba_df <- terra::as.data.frame(pinyonba_rast2,
                                    xy = TRUE,
                                    cells = TRUE)
#switch back - since everything else works
sf_use_s2(TRUE)

#get one year of basal area data
ba15 <- pinyonba_df %>%
  mutate(PinyonBA_sqftPerAc_2015 = case_when(PinyonBA_sqftPerAc_2015 == 0 ~ NA_real_,
                                             TRUE ~ PinyonBA_sqftPerAc_2015))

# Plot full range ---------------------------------------------------------

# ggplot()+
#   # geom_tile(data = cone15, 
#   #           aes(x = x, 
#   #               y = y, 
#   #               fill = `2011`)) +
#   geom_tile(data = ba15, 
#             aes(x = x, 
#                 y = y, 
#                 fill = PinyonBA_sqftPerAc_2015)) +
#     geom_sf(data = states, fill = "white", alpha =.2) +
#   # scale_fill_viridis_c(na.value = 'transparent',
#   #                      limits = c(5e-9, 22),
#   #                      breaks = c(5e-9,22),
#   #                      labels = c('low', 'high')) +
#   scale_fill_continuous(na.value = 'transparent',
#                         low = 'grey',
#                         high = 'grey') +
#   #scale_fill_viridis_c(na.value = 'transparent') +
#   geom_sf(data = ebird, shape = 1, alpha = 0.5) +
#     labs(fill = "Pinyon basal area") +
#   theme(legend.position = "none")

(total_range_map <- ggplot() +
    geom_sf(data = range2, 
            fill = "#df65b0", 
            color = 'transparent') +
    geom_sf(data = states2, 
            fill = 'transparent',
            color = 'grey58') +
  geom_sf(data = mexico, 
          fill = 'transparent',
          color= 'grey69') +
  geom_sf(data = canada, 
          fill = 'transparent',
          color= 'grey69') +
  geom_sf(data = sw2, alpha = 0.01,
          color = 'grey35',
          linewidth = 0.4) +
  #coord_sf(xlim = c(-125,-60), ylim = c(25, 55)) +
  coord_sf(crs = 'EPSG:5070',
           xlim = c(-2271217.626143, 2201865.603405),
           ylim = c(300000,  3245749.91))) 

ggsave(filename = here('pictures',
                       'R',
                       'range_map.pdf'),
       height = 2, 
       width = 2,
       units = "in")

# Plot range in SW --------------------------------------------------------

intersect_pct <- st_intersection(range2, sw2) 

sw_range_map <- ggplot()+
  geom_sf(data = intersect_pct,
          fill = "#df65b0", 
          color = 'transparent') +
  geom_sf(data = states, 
          fill = "transparent") +
  coord_sf() 

#how much of their range is covered by this region??
st_area(intersect_pct)/st_area(range2)

# Basal area map ----------------------------------------------------------

ba_map <- ggplot() +
  geom_tile(data = ba15, 
            aes(x = x, 
                y = y, 
                fill = PinyonBA_sqftPerAc_2015)) +
  scale_fill_continuous(na.value = 'transparent',
                        low = 'grey',
                        high = 'grey') +
  geom_sf(data = states, 
          fill = 'transparent') +
  coord_sf() +
  theme_void() +
  theme(legend.position = "none")

# Sampling density hex map ------------------------------------------------

ebird_coords <- st_coordinates(ebird) %>%
  as.data.frame() %>%
  rename(lon = X, lat = Y)

(hex_map <- ggplot() +
  geom_sf(data = states, fill = "white", alpha =.2) +
  geom_hex(data = ebird_coords, 
           binwidth = 0.3,
           aes(x = lon, 
               y = lat, 
               fill = stat(count))) +
  scale_fill_viridis_c(name = "Sampling density<br> (checklists km<sup>-2</sup>)",
                       breaks = c(77.9, 233.8, 389.7),
                       labels = c(0.1, 0.3, 0.5)) +
  coord_sf() +
  theme_void() +
  theme(legend.title = element_markdown(size = 8),
        #legend.position = "bottom",
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"),
        legend.key.height = unit(0.2, 'cm')) +
    guides(fill = guide_colorbar(title.position = "top")))


# Plot together -----------------------------------------------------------

total_range_map/ba_map/hex_map+
  plot_annotation(tag_levels = "a",
                  tag_suffix = ")")

ggsave(filename = here('pictures',
                       'final',
                       'maps.jpg'),
       height = 5, 
       width = 3,
       units = "in")

total_range_map+ba_map+hex_map+
  plot_annotation(tag_levels = "a",
                  tag_suffix = ")")

ggsave(filename = here('pictures',
                       'final',
                       'maps_wide.jpg'),
       height = 2, 
       width = 5,
       units = "in")

