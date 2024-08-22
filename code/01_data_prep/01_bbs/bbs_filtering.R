# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Notes on data filtering -------------------------------------------------

#Filter to just include pinyon jay observations
#pinyon jay AOU: 04920
#filter to only include colorado and new mexico
#STate Nums: 17 (Colorado) and 60 (New Mexico)
#Arizona: 06, Utah: 85

# Load datasets -----------------------------------------------------------

#Dataset structure:
#50 stops on each BBS route (0.5 miles apart, ~24.5 miles total on a route?)
#start with the 50-stop data (only available for 1997-onward)

#colorado in fifty2.csv
#new mexico in fifty6.csv
#arizona in fifty1.csv
#utah in fifty9.csv

#get colorado data read in
CO <- read.csv(here('data',
                     'bbs_data',
                     '2023Release_Nor',
                     '50-StopData',
                     'fifty2.csv'))

#new mexico data
NM <- read.csv(here('data',
                    'bbs_data',
                    '2023Release_Nor',
                    '50-StopData',
                    'fifty6.csv'))

#AZ data
AZ <- read.csv(here('data',
                    'bbs_data',
                    '2023Release_Nor',
                    '50-StopData',
                    'fifty1.csv'))

#UT data
UT <- read.csv(here('data',
                    'bbs_data',
                    '2023Release_Nor',
                    '50-StopData',
                    'fifty9.csv'))

#observer IDs for each dataset - for experience covariate
obs <- read.csv(here('data',
                     'bbs_data',
                     '2023Release_Nor',
                     'weather.csv'))


# Filter colorado data ----------------------------------------------------

#get all the routes and years for colorado
CO_routes <- CO %>%
  filter(StateNum == 17) %>%
  distinct(RouteDataID, CountryNum, StateNum, Route, RPID, Year)

#get all the routes/years where jays were observed
CO_jay <- CO %>%
  filter(StateNum == 17) %>%
  filter(AOU == 4920)

#jays seen on 432/2737 route/years (16%)

#blend together and fill all route/years without Jay observations with 0
#across all stops on that route
CO_all <- CO_jay %>%
  full_join(CO_routes, by = c("RouteDataID", "CountryNum", "StateNum",
                              "Route", "RPID", "Year")) %>%
  replace_na(list(AOU = 4940)) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  mutate(State = "CO")


# Filter new mexico data --------------------------------------------------

#get all the routes and years for colorado
NM_routes <- NM %>%
  filter(StateNum == 60) %>%
  distinct(RouteDataID, CountryNum, StateNum, Route, RPID, Year)

#get all the routes/years where jays were observed
NM_jay <- NM %>%
  filter(StateNum == 60) %>%
  filter(AOU == 4920)

#jays observed on 466/1478 route/years (32%)

#blend together and fill all route/years without Jay observations with 0
#across all stops on that route
NM_all <- NM_jay %>%
  full_join(NM_routes, by = c("RouteDataID", "CountryNum", "StateNum",
                              "Route", "RPID", "Year")) %>%
  replace_na(list(AOU = 4940)) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  mutate(State = "NM")


# Filter AZ data ----------------------------------------------------------

#get all the routes and years for colorado
AZ_routes <- AZ %>%
  filter(StateNum == 6) %>%
  distinct(RouteDataID, CountryNum, StateNum, Route, RPID, Year)

#get all the routes/years where jays were observed
AZ_jay <- AZ %>%
  filter(StateNum == 6) %>%
  filter(AOU == 4920)

#jays observed on 145/1129 route/years (13%)

#blend together and fill all route/years without Jay observations with 0
#across all stops on that route
AZ_all <- AZ_jay %>%
  full_join(AZ_routes, by = c("RouteDataID", "CountryNum", "StateNum",
                              "Route", "RPID", "Year")) %>%
  replace_na(list(AOU = 4940)) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  mutate(State = "AZ")


# Filter UT data ----------------------------------------------------------

#get all the routes and years for colorado
UT_routes <- UT %>%
  filter(StateNum == 85) %>%
  distinct(RouteDataID, CountryNum, StateNum, Route, RPID, Year)

#get all the routes/years where jays were observed
UT_jay <- UT %>%
  filter(StateNum == 85) %>%
  filter(AOU == 4920)

#jays observed on 623/1901 route/years (33%)

#blend together and fill all route/years without Jay observations with 0
#across all stops on that route
UT_all <- UT_jay %>%
  full_join(UT_routes, by = c("RouteDataID", "CountryNum", "StateNum",
                              "Route", "RPID", "Year")) %>%
  replace_na(list(AOU = 4940)) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  mutate(State = "UT")

# Combine and export ------------------------------------------------------

BBS_all <- CO_all %>%
  rbind(NM_all, AZ_all, UT_all)


# Get observer IDs and experience for each of these -----------------------

obs_exp <- obs %>%
  distinct(Year, ObsN) %>%
  group_by(ObsN) %>%
  arrange(Year) %>%
  mutate(ObserverExp = 1:n()) %>%
  ungroup() 

obs2 <- obs %>%
  left_join(obs_exp, by = c("Year", "ObsN")) %>%
  dplyr::select(RouteDataID, StateNum, Route,
                Year, ObsN, ObserverExp)

# combine with bbs data ---------------------------------------------------

BBS_all2 <- BBS_all %>%
  left_join(obs2, by = c("Year", "RouteDataID", "StateNum", "Route"))


write.csv(BBS_all2, here('data',
                        'bbs_data',
                        'cleaned_data',
                        'pjay_data_co_nm_az_ut.csv'))


# Explore -----------------------------------------------------------------

BBS_all %>%
  pivot_longer(Stop1:Stop50,
               names_to = "Stop", 
               values_to = "count") %>%
  ggplot() +
  geom_histogram(aes(x = count)) +
  scale_y_log10() +
  theme_bw() +
  labs(x = "Number of birds observed", 
       y = "Number of points")

BBS_all %>%
  pivot_longer(Stop1:Stop50,
               names_to = "Stop", 
               values_to = "count") %>%
  group_by(Year, Route) %>%
  summarise(total = sum(count)) %>%
  ggplot() +
  geom_histogram(aes(x = total)) +
  scale_y_log10() +
  theme_bw() +
  labs(x = "Number of birds observed", 
       y = "Number of route-years")

  