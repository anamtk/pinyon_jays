
# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "sf",  
                  "terra",
                  'readxl')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

az <- read.csv(here('data',
                 'bbs_data',
                 '2023Release_Nor',
                 'States',
                 'Arizona.csv')) 

co <- read.csv(here('data',
                    'bbs_data',
                    '2023Release_Nor',
                    'States',
                    'Colorad.csv')) 

nm <- read.csv(here('data',
                    'bbs_data',
                    '2023Release_Nor',
                    'States',
                    'NMexico.csv')) 

ut <- read.csv(here('data',
                    'bbs_data',
                    '2023Release_Nor',
                    'States',
                    'Utah.csv')) 

routes <- read.csv(here('data',
                        'bbs_data',
                        '2023Release_Nor',
                        'routes.csv'))


# Get presence-absence for jays for 10-stops ------------------------------

#get all the routes and years for colorado
CO_routes <- co %>%
  distinct(RouteDataID, CountryNum, StateNum, Route, RPID, Year)

#get all the routes/years where jays were observed
CO_jay <- co %>%
  filter(AOU == 4920)

#blend together and fill all route/years without Jay observations with 0
#across all stops on that route
CO_all <- CO_jay %>%
  full_join(CO_routes, by = c("RouteDataID", "CountryNum", "StateNum",
                              "Route", "RPID", "Year")) %>%
  replace_na(list(AOU = 4940)) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  dplyr::select(RouteDataID, CountryNum, StateNum,
                Route, RPID, Year, Count10) %>%
  mutate(state = "CO")


#get all the routes and years for az
AZ_routes <- az %>%
  distinct(RouteDataID, CountryNum, StateNum, Route, RPID, Year)

#get all the routes/years where jays were observed
AZ_jay <- az %>%
  filter(AOU == 4920)

#blend together and fill all route/years without Jay observations with 0
#across all stops on that route
AZ_all <- AZ_jay %>%
  full_join(AZ_routes, by = c("RouteDataID", "CountryNum", "StateNum",
                              "Route", "RPID", "Year")) %>%
  replace_na(list(AOU = 4940)) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  dplyr::select(RouteDataID, CountryNum, StateNum,
                Route, RPID, Year, Count10)%>%
  mutate(state = "AZ")


#get all the routes and years for ut
UT_routes <- ut %>%
  distinct(RouteDataID, CountryNum, StateNum, Route, RPID, Year)

#get all the routes/years where jays were observed
UT_jay <- ut %>%
  filter(AOU == 4920)

#blend together and fill all route/years without Jay observations with 0
#across all stops on that route
UT_all <- UT_jay %>%
  full_join(UT_routes, by = c("RouteDataID", "CountryNum", "StateNum",
                              "Route", "RPID", "Year")) %>%
  replace_na(list(AOU = 4940)) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  dplyr::select(RouteDataID, CountryNum, StateNum,
                Route, RPID, Year, Count10)%>%
  mutate(state = "UT")

#get all the routes and years for nm
NM_routes <- nm %>%
  distinct(RouteDataID, CountryNum, StateNum, Route, RPID, Year)

#get all the routes/years where jays were observed
NM_jay <- nm %>%
  filter(AOU == 4920)

#blend together and fill all route/years without Jay observations with 0
#across all stops on that route
NM_all <- NM_jay %>%
  full_join(NM_routes, by = c("RouteDataID", "CountryNum", "StateNum",
                              "Route", "RPID", "Year")) %>%
  replace_na(list(AOU = 4940)) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  dplyr::select(RouteDataID, CountryNum, StateNum,
                Route, RPID, Year, Count10)%>%
  mutate(state = "NM")


# Combine -----------------------------------------------------------------

all_bbs <- CO_all %>%
  bind_rows(NM_all, UT_all, AZ_all) %>%
  left_join(routes, by = c("CountryNum", "StateNum",
                           "Route")) %>%
  dplyr::select(RouteDataID, CountryNum,
                StateNum, Route, RPID, Year,
                Count10, state, Latitude,
                Longitude)


# Import cone dataset for filtering bbs -----------------------------------

cones <- terra::rast(here("data",
                          "spatial_data", 
                          "masting_data",
                          "quantile_pied_predictions_scaled.tif"))

# Set the CRS for birds -------------------------------------------

#these data are in WGS84, so I'll create this crs to call later
crs_wgs <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

# Make spatial ------------------------------------------------------------

#convert all bbs data to an sf object
bbs_spatial <- st_as_sf(all_bbs, coords = c("Longitude", "Latitude"),
                        crs = st_crs(crs_wgs))

# Extract cell ID from raster ------------------------------------------------------------

##BBS
bbs_cellIDs <- terra::extract(cones, vect(bbs_spatial), cells = T)

all_bbs$cellID <- bbs_cellIDs$cell

# Turn cone raster into df with cell IDs ---------------------------------------

cone_df <- terra::as.data.frame(cones, 
                                xy = TRUE,
                                cells = TRUE)


# Filter datasets for cells with cone data --------------------------------

ids <- cone_df %>%
  distinct(cell)

all_bbs2 <- all_bbs %>%
  filter(cellID %in% ids$cell) %>%
  rename(Count = Count10)

write.csv(all_bbs2, here('data',
                         'spatial_data',
                         'cleaned_data',
                         'allbbsyears_10count.csv'))


