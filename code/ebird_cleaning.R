
# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "auk")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

ebdco <- auk_ebd(here('data', 'ebird_data',
                    'ebd_US-CO_pinjay_smp_relMay-2024',
                    "ebd_US-CO_pinjay_smp_relMay-2024.txt"), 
               file_sampling = here('data', 'ebird_data',
                                    'ebd_US-CO_pinjay_smp_relMay-2024',
                                    "ebd_US-CO_pinjay_smp_relMay-2024_sampling.txt"))

ebdco2 <- ebdco %>% 
  # restrict to the standard traveling and stationary count protocols
  auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
  auk_complete()

data_dir <- "data"
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}
f_ebd <- file.path(data_dir, "ebd_woothr_june_bcr27.txt")
f_sampling <- file.path(data_dir, "ebd_checklists_june_bcr27.txt")

# only run if the files don't already exist
if (!file.exists(f_ebd)) {
  auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling)
}

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


# Select columns of interest from all dfs ---------------------------------

# Observation datasets ----------------------------------------------------

#observation datasets:
ebirdco1 <- ebirdco %>%  
  filter(`ALL SPECIES REPORTED` == 1) %>% # (1 = yes, 0 = no)
  dplyr::select(`SAMPLING EVENT IDENTIFIER`) %>% 
  mutate(presence = 1)

ebirdnm1 <- ebirdnm %>%  
  filter(`ALL SPECIES REPORTED` == 1) %>% # (1 = yes, 0 = no)
  dplyr::select(`SAMPLING EVENT IDENTIFIER`) %>% 
  mutate(presence = 1)



# Checklist datasets ------------------------------------------------------

checkco1 <- checkco %>%
  dplyr::select(`SAMPLING EVENT IDENTIFIER`,
                `STATE CODE`, 
                COUNTY,
                `BCR CODE`,
                `LOCALITY`,
                `LOCALITY ID`,
                LATITUDE,
                LONGITUDE,
                `TIME OBSERVATIONS STARTED`,
                `OBSERVER ID`,
                `OBSERVATION DATE`,
                `PROTOCOL TYPE`,
                `DURATION MINUTES`,
                `EFFORT DISTANCE KM`,
                `NUMBER OBSERVERS`, 
                `ALL SPECIES REPORTED`) %>% # (1 = yes, 0 = no)
  filter(`ALL SPECIES REPORTED` == 1) 

all_co <- ebirdco1 %>%
  full_join(checkco1, by = 'SAMPLING EVENT IDENTIFIER') %>%
  replace_na(list(presence = 0))
#1922076 checklists

all_co2 <- all_co %>%
  filter(!is.na(`DURATION MINUTES`)) 
#1848096 checklists

all_co3 <- all_co2 %>%
  filter(!is.na(`TIME OBSERVATIONS STARTED`)) 
#1846945 checklists

all_co4 <- all_co3 %>%
  filter(`PROTOCOL TYPE` == "Traveling") 
#1221976 checklists

all_co5 <- all_co4 %>%
  mutate(year = str_sub(`OBSERVATION DATE`, 1, 4))

yearly_co <- all_co5 %>%
  mutate(presence = as.factor(presence)) %>%
  group_by(year, presence) %>%
  tally() %>%
  pivot_wider(names_from = presence,
              values_from = n)

#NM
checknm1 <- checknm %>%
  dplyr::select(`SAMPLING EVENT IDENTIFIER`,
                `STATE CODE`, 
                COUNTY,
                `BCR CODE`,
                `LOCALITY`,
                `LOCALITY ID`,
                LATITUDE,
                LONGITUDE,
                `TIME OBSERVATIONS STARTED`,
                `OBSERVER ID`,
                `OBSERVATION DATE`,
                `PROTOCOL TYPE`,
                `DURATION MINUTES`,
                `EFFORT DISTANCE KM`,
                `NUMBER OBSERVERS`, 
                `ALL SPECIES REPORTED`) %>% # (1 = yes, 0 = no)
  filter(`ALL SPECIES REPORTED` == 1) 

all_nm <- ebirdnm1 %>%
  full_join(checknm1, by = 'SAMPLING EVENT IDENTIFIER') %>%
  replace_na(list(presence = 0))
#743823 checklists

all_nm2 <- all_nm %>%
  filter(!is.na(`DURATION MINUTES`)) 
#709697 checklists

all_nm3 <- all_nm2 %>%
  filter(!is.na(`TIME OBSERVATIONS STARTED`)) 
#709353 checklists

all_nm4 <- all_nm3 %>%
  filter(`PROTOCOL TYPE` == "Traveling") 
#453221 checklists

all_nm5 <- all_nm4 %>%
  mutate(year = str_sub(`OBSERVATION DATE`, 1, 4))

yearly_nm <- all_nm5 %>%
  mutate(presence = as.factor(presence)) %>%
  group_by(year, presence) %>%
  tally() %>%
  pivot_wider(names_from = presence,
              values_from = n)


# Subsample for occupancy modeling ----------------------------------------

#Following tutorial here:
#https://cornelllabofornithology.github.io/ebird-best-practices/occupancy.html

#for this, I'm keeping sites where they've been sampled at least 2x/year in
#one or more years (but not necessarily all years)
#Subsettign sites with >10 observations in a year so that we only have 
#10 of those sites (so they're not over-represented)

#CO 
#figure out how often a site has been sampled within a year
all_co6 <- all_co5 %>%
  group_by(year, `LOCALITY ID`) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  group_by(`LOCALITY ID`) %>%
  filter(any(count > 1)) %>%
  ungroup()
  
#subsample when there's lots of sampling in a site and year
largeco <- all_co6 %>%
  group_by(`LOCALITY ID`, year) %>%
  filter(count > 9) %>%
  sample_n(10) %>%
  ungroup()

smallco <- all_co6 %>%
  group_by(`LOCALITY ID`, year) %>%
  filter(count < 10) %>%
  ungroup()

all_co7 <- largeco %>%
  bind_rows(smallco)
#428670 observations

#NM

#figure out how often a site has been sampled within a year
all_nm6 <- all_nm5 %>%
  group_by(year, `LOCALITY ID`) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  group_by(`LOCALITY ID`) %>%
  filter(any(count > 1)) %>%
  ungroup()

#subsample when there's lots of sampling in a site and year
largenm <- all_nm6 %>%
  group_by(`LOCALITY ID`, year) %>%
  filter(count > 9) %>%
  sample_n(10) %>%
  ungroup()

smallnm <- all_nm6 %>%
  group_by(`LOCALITY ID`, year) %>%
  filter(count < 10) %>%
  ungroup()

all_nm7 <- largenm %>%
  bind_rows(smallnm)
#167722

