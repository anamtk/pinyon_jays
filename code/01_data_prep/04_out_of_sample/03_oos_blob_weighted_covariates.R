#Get weighted average covariates per ebird blob
#Ana Miller-ter Kuile
#November 13, 2024

#this script calculates weighted average covariate effects for 
#each "blob" of points in the ebird dataset based on the coverage 
#of the blob on different grid cells


# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse", 
                  "sf",  
                  "terra",
                  'readxl',
                  'sf')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

ebird <- read.csv(here('data',
                       'ebird_data',
                       'cleaned_data',
                       '05_oos',
                       'ebird_oos_cellIDlists.csv')) %>%
  dplyr::select(-X) %>%
  rename(blobID = cellID) %>%
  rename(cellID = cell) 

cells <- ebird %>%
  distinct(cellID)

cones <- read.csv(here('data',
                        'spatial_data',
                        'cleaned_data',
                       '01_get_gridIDs',
                        'cone_masting_df.csv'))%>%
  dplyr::select(cell, 
                X2000:X2024) %>%
  rename(cellID = cell) %>%
  filter(cellID %in% cells$cellID)


temp <- read.csv(here('data',
                         'spatial_data',
                         'cleaned_data',
                      '01_get_gridIDs',
                         'temp_data_df.csv')) %>%
  dplyr::select(cellID, 
                PRISM_tmax_stable_4kmM3_200001_bil:PRISM_tmax_stable_4kmM3_202312_bil) %>%
  filter(cellID %in% cells$cellID)

tmean <- read.csv(here('data',
                       'spatial_data',
                       'cleaned_data',
                       '01_get_gridIDs',
                       'tmean_data_df.csv')) %>%
  dplyr::select(cellID, 
                PRISM_tmean_stable_4kmM3_200001_bil:PRISM_tmean_stable_4kmM3_202312_bil) %>%
  filter(cellID %in% cells$cellID)


ppt <- read.csv(here('data',
                        'spatial_data',
                        'cleaned_data',
                     '01_get_gridIDs',
                        'ppt_data_df.csv')) %>%
  dplyr::select(cellID, 
                PRISM_ppt_stable_4kmM3_200001_bil:PRISM_ppt_stable_4kmM3_202312_bil)%>%
  filter(cellID %in% cells$cellID)

monsoon <- read.csv(here('data',
                            'spatial_data',
                            'cleaned_data',
                         '01_get_gridIDs',
                            'monsoon_data_df.csv')) %>%
  dplyr::select(cellID, PRISM_ppt_30yr_normal_800mM4_07_bil) %>%
  rename(monsoon = PRISM_ppt_30yr_normal_800mM4_07_bil) %>%
  filter(cellID %in% cells$cellID)

pinyon <- read.csv(here('data',
                            'spatial_data',
                            'cleaned_data',
                        '01_get_gridIDs',
                            'pinyonba_data_df.csv'))%>%
  dplyr::select(cell, 
                PinyonBA_sqftPerAc_2010:PinyonBA_sqftPerAc_2022) %>%
  rename(cellID = cell) %>%
  filter(cellID %in% cells$cellID)

# Join climate data with ebird data ---------------------------------------


# Cones -------------------------------------------------------------------


cones2 <- ebird %>%
  left_join(cones, by = "cellID") %>%
  pivot_longer(X2000:X2024,
               names_to = 'cone_yr',
               values_to = "cones") %>%
  group_by(blobnum, year, cone_yr) %>%
  summarise(wt = stats::weighted.mean(cones, coverage_fraction, na.rm = T)) %>%
  ungroup() %>%
  mutate(cone_yr = str_sub(cone_yr, 2, nchar(cone_yr)),
         cone_yr = as.numeric(cone_yr)) %>%
  group_by(blobnum, year)%>%
  rowwise() %>%
  mutate(lagnum = year - cone_yr) %>%
  ungroup() %>%
  filter(lagnum %in% c(-1:3)) %>%
  #fix lag indexing - 4 is current time
  mutate(lag = case_when(lagnum == -1 ~ 1,
                         lagnum == -0 ~ 2,
                         lagnum == 1 ~ 3,
                         lagnum == 2 ~ 4,
                         lagnum == 3 ~ 5)) %>%
  dplyr::select(blobnum, year, wt, lag) %>%
  rename(cones = wt)

write.csv(cones2, here('data',
              'spatial_data',
              'cleaned_data',
              '03_oos',
              'oos_cones_weighted_mean_blob.csv'))

# temp --------------------------------------------------------------------

#seasons for climate variables::
#from Wiggins (2005) report:
#breeding: late Feb-early May
#May-June: dependent young
#August-Feb: fall and early winter: wander to find resources
#what do they do in July??? 
#thus, seasons would be months:
#Season 1: 2-4: breeding
#Season 2: 5-6: feeding dependent young
#Season 3: 7: unknown, just chillin?
#Season 4: 8-1: winter foraging, potentially seeking food elsewhere

#months 2-5 is when we're sampling ebird, so maybe months
#starting in the 
temp2 <- ebird %>%
  left_join(temp, by = c("cellID")) %>%
  pivot_longer(PRISM_tmax_stable_4kmM3_200001_bil:PRISM_tmax_stable_4kmM3_202312_bil,
               names_to = 'temp_yrmonth',
               values_to = "temp") %>%
  mutate(temp_yr = str_sub(temp_yrmonth, 25, (nchar(temp_yrmonth))-6),
         temp_month = str_sub(temp_yrmonth, 29, (nchar(temp_yrmonth))-4),
         temp_yr = as.numeric(temp_yr),
         temp_month = as.numeric(temp_month)) %>%
  group_by(blobnum) %>%
  filter(temp_yr %in% c((year-3), (year-2), (year-1), year, (year+1))) %>%
  dplyr::select(-temp_yrmonth) %>%
  pivot_wider(names_from = temp_month,
              values_from = temp) %>%
  arrange(blobnum, temp_yr) %>%
  #determine next year january value to compute temp season 4 below
  mutate(lead_1 = lead(`1`)) %>%
  filter(temp_yr %in% c((year-3), (year-2), (year-1), year)) %>%
  rowwise() %>%
  #get each seasonal VPD
  mutate(temp_1 = max(`2`, `3`, `4`),
         temp_2 = max(`5`, `6`),
         temp_3 = `7`,
         temp_4 = max(`8`, `9`, `10`, `11`, `12`, lead_1, na.rm = T)) %>%
  dplyr::select(-c(`1`:`12`), -lead_1) %>%
  ungroup() %>%
  pivot_longer(temp_1:temp_4,
               names_to = "season",
               values_to = "temp") %>%
  mutate(season =str_sub(season, 6, nchar(season)),
         season = as.numeric(season)) %>%
  group_by(blobnum,year, temp_yr,season) %>%
  summarise(wt = stats::weighted.mean(temp, coverage_fraction, na.rm = T)) %>%
  ungroup() %>%
  group_by(blobnum) %>%
  filter((!(season %in% c(2,3,4) & temp_yr == year))) %>%
  arrange(blobnum, -temp_yr, -season) %>%
  mutate(lag = 1:n()) %>%
  ungroup() %>%
  dplyr::select(blobnum, year, season, wt, lag) %>%
  rename(temp = wt)
  
write.csv(temp2, here('data',
                       'spatial_data',
                       'cleaned_data',
                      '03_oos',
                       'oos_temp_weighted_mean_blob.csv'))

# tmean --------------------------------------------------------------------

#seasons for climate variables::
#from Wiggins (2005) report:
#breeding: late Feb-early May
#May-June: dependent young
#August-Feb: fall and early winter: wander to find resources
#what do they do in July??? 
#thus, seasons would be months:
#Season 1: 2-4: breeding
#Season 2: 5-6: feeding dependent young
#Season 3: 7: unknown, just chillin?
#Season 4: 8-1: winter foraging, potentially seeking food elsewhere

#months 2-5 is when we're sampling ebird, so maybe months
#starting in the 
tmean2 <- ebird %>%
  left_join(tmean, by = c("cellID")) %>%
  pivot_longer(PRISM_tmean_stable_4kmM3_200001_bil:PRISM_tmean_stable_4kmM3_202312_bil,
               names_to = 'temp_yrmonth',
               values_to = "temp") %>%
  mutate(temp_yr = str_sub(temp_yrmonth, 26, (nchar(temp_yrmonth))-6),
         temp_month = str_sub(temp_yrmonth, 30, (nchar(temp_yrmonth))-4),
         temp_yr = as.numeric(temp_yr),
         temp_month = as.numeric(temp_month)) %>%
  group_by(blobnum) %>%
  filter(temp_yr %in% c((year-3), (year-2), (year-1), year, (year+1))) %>%
  dplyr::select(-temp_yrmonth) %>%
  pivot_wider(names_from = temp_month,
              values_from = temp) %>%
  arrange(blobnum, temp_yr) %>%
  #determine next year january value to compute temp season 4 below
  mutate(lead_1 = lead(`1`)) %>%
  filter(temp_yr %in% c((year-3), (year-2), (year-1), year)) %>%
  rowwise() %>%
  #get each seasonal VPD
  mutate(temp_1 = mean(`2`, `3`, `4`),
         temp_2 = mean(`5`, `6`),
         temp_3 = `7`,
         temp_4 = mean(`8`, `9`, `10`, `11`, `12`, lead_1, na.rm = T)) %>%
  dplyr::select(-c(`1`:`12`), -lead_1) %>%
  ungroup() %>%
  pivot_longer(temp_1:temp_4,
               names_to = "season",
               values_to = "temp") %>%
  mutate(season =str_sub(season, 6, nchar(season)),
         season = as.numeric(season)) %>%
  group_by(blobnum,year, temp_yr,season) %>%
  summarise(wt = stats::weighted.mean(temp, coverage_fraction, na.rm = T)) %>%
  ungroup() %>%
  group_by(blobnum) %>%
  filter((!(season %in% c(2,3,4) & temp_yr == year))) %>%
  arrange(blobnum, -temp_yr, -season) %>%
  mutate(lag = 1:n()) %>%
  ungroup() %>%
  dplyr::select(blobnum, year, season, wt, lag) %>%
  rename(temp = wt)

write.csv(tmean2, here('data',
                      'spatial_data',
                      'cleaned_data',
                      '03_oos',
                      'oos,tmean_weighted_mean_blob.csv'))

# PPT ---------------------------------------------------------------------

ppt2 <- ebird %>%
  left_join(ppt, by = c("cellID")) %>%
  pivot_longer(PRISM_ppt_stable_4kmM3_200001_bil:PRISM_ppt_stable_4kmM3_202312_bil,
               names_to = 'ppt_yrmonth',
               values_to = "ppt") %>%
  mutate(ppt_yr = str_sub(ppt_yrmonth, 24, (nchar(ppt_yrmonth))-6),
         ppt_month = str_sub(ppt_yrmonth, 28, (nchar(ppt_yrmonth))-4),
         ppt_yr = as.numeric(ppt_yr),
         ppt_month = as.numeric(ppt_month)) %>%
  group_by(blobnum) %>%
  filter(ppt_yr %in% c((year-3), (year-2), (year-1), year, (year+1))) %>%
  dplyr::select(-ppt_yrmonth) %>%
  pivot_wider(names_from = ppt_month,
              values_from = ppt) %>%
  arrange(blobnum, ppt_yr) %>%
  #determine next year january value to compute ppt season 4 below
  mutate(lead_1 = lead(`1`)) %>%
  filter(ppt_yr %in% c((year-3), (year-2), (year-1), year)) %>%
  rowwise() %>%
  #get each seasonal VPD
  mutate(ppt_1 = max(`2`, `3`, `4`),
         ppt_2 = max(`5`, `6`),
         ppt_3 = `7`,
         ppt_4 = max(`8`, `9`, `10`, `11`, `12`, lead_1, na.rm = T)) %>%
  dplyr::select(-c(`1`:`12`), -lead_1) %>%
  ungroup() %>%
  pivot_longer(ppt_1:ppt_4,
               names_to = "season",
               values_to = "ppt") %>%
  mutate(season =str_sub(season, 5, nchar(season)),
         season = as.numeric(season)) %>%
  group_by(blobnum,year, ppt_yr,season) %>%
  summarise(wt = stats::weighted.mean(ppt, coverage_fraction, na.rm = T)) %>%
  ungroup() %>%
  group_by(blobnum) %>%
  filter((!(season %in% c(2,3,4) & ppt_yr == year))) %>%
  arrange(blobnum, -ppt_yr, -season) %>%
  mutate(lag = 1:n()) %>%
  ungroup() %>%
  dplyr::select(blobnum, year, season, wt, lag) %>%
  rename(ppt = wt)

write.csv(ppt2, here('data',
                      'spatial_data',
                      'cleaned_data',
                     '03_oos',
                      'oos_ppt_weighted_mean_blob.csv'))


# Monsoon -----------------------------------------------------------------

monsoon2 <- ebird %>%
  left_join(monsoon, by = "cellID") %>%
  group_by(blobnum, year) %>%
  summarise(wt = stats::weighted.mean(monsoon, coverage_fraction, na.rm = T)) %>%
  ungroup()
  
write.csv(monsoon2, here('data',
                     'spatial_data',
                     'cleaned_data',
                     '03_oos',
                     'oos_monsoon_weighted_mean_blob.csv'))


# Pinyon BA ---------------------------------------------------------------

pinyon2 <- ebird %>%
  left_join(pinyon, by = "cellID") %>%
  pivot_longer(PinyonBA_sqftPerAc_2010:PinyonBA_sqftPerAc_2022,
               names_to = 'pinyon_yr',
               values_to = "pinyon") %>%
  mutate(pinyon_yr = str_sub(pinyon_yr, 20, (nchar(pinyon_yr))),
         pinyon_yr = as.numeric(pinyon_yr)) %>%
  filter(pinyon_yr == year) %>%
  group_by(blobnum, year) %>%
  summarise(wt = stats::weighted.mean(pinyon, coverage_fraction, na.rm = T)) %>%
  ungroup()
  
write.csv(pinyon2, here('data',
                         'spatial_data',
                         'cleaned_data',
                        '03_oos',
                         'oos_pinyonBA_weighted_mean_blob.csv'))
  
  
  
