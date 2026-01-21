# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "auk", 'lubridate')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Extract using auk -------------------------------------------------------

#FOllowing protocol here:
#https://cornelllabofornithology.github.io/ebird-best-practices/ebird.html

# Colorado ----------------------------------------------------------------


ebd_co <- auk_ebd(file = here('data',
                              '01_ebird_data',
                              'ebd_US-CO_pinjay_smp_relMay-2024',
                              'ebd_US-CO_pinjay_smp_relMay-2024.txt'),
                  file_sampling = here('data',
                                       'ebird_data',
                                       'ebd_US-CO_pinjay_smp_relMay-2024',
                                       'ebd_US-CO_pinjay_smp_relMay-2024_sampling.txt'))

#define the filters I will want
ebd_filters <- ebd_co %>%
  #filter just the months of interest, I think
  #we will want june and july, but need to verify
  #auk_date(date = c("*-06-01", "*-07-31")) %>%
  #filtering without bbs will be the nesting season
  #"late Feb to early May"
  auk_date(date = c("*-02-15", "*-05-15")) %>%
#restirct to the stationary and traveling protocols
  auk_protocol(protocol = c("Stationary", "Traveling")) %>%
  #less than 5 h
  auk_duration(duration = c(0, 300)) %>%
  # less than 5km distance traveled
  auk_distance(distance = c(0, 5)) %>%
  #less than 10 observers
  #not sure how to do this one...
  #only use complete checklists
  auk_complete()

f_ebd_co <- file.path(here('data',
                           '01_ebird_data',
                           'cleaned_data',
                           '01_auk_filtering',
                           'ebd_co_filtered.txt'))
f_sampling_co <- file.path(here('data',
                            '01_ebird_data',
                            'cleaned_data',
                            '01_auk_filtering',
                            'ebd_co_checklists.txt'))

# only run if the files don't already exist
if (!file.exists(f_ebd_co)) {
  auk_filter(ebd_filters, file = f_ebd_co, file_sampling = f_sampling_co)
}

# New Mexico --------------------------------------------------------------

ebd_nm <- auk_ebd(file = here('data',
                              '01_ebird_data',
                              'ebd_US-NM_pinjay_smp_relMay-2024',
                              'ebd_US-NM_pinjay_smp_relMay-2024.txt'),
                  file_sampling = here('data',
                                       '01_ebird_data',
                                       'ebd_US-NM_pinjay_smp_relMay-2024',
                                       'ebd_US-NM_pinjay_smp_relMay-2024_sampling.txt'))

#define the filters I will want
ebd_filters <- ebd_nm %>%
  #filter just the months of interest, I think
  #we will want june and july, but need to verify
  #auk_date(date = c("*-06-01", "*-07-31")) %>%
  #filtering without bbs will be the nesting season
  #"late Feb to early May"
  auk_date(date = c("*-02-15", "*-05-15")) %>%
  #restirct to the stationary and traveling protocols
  auk_protocol(protocol = c("Stationary", "Traveling")) %>%
  #less than 5 h
  auk_duration(duration = c(0, 300)) %>%
  # less than 5km distance traveled
  auk_distance(distance = c(0, 5)) %>%
  #less than 10 observers
  #not sure how to do this one...
  #only use complete checklists
  auk_complete()

f_ebd_nm <- file.path(here('data',
                           '01_ebird_data',
                           'cleaned_data',
                           '01_auk_filtering',
                           'ebd_nm_filtered.txt'))
f_sampling_nm <- file.path(here('data',
                                '01_ebird_data',
                                'cleaned_data',
                                '01_auk_filtering',
                                'ebd_nm_checklists.txt'))

# only run if the files don't already exist
if (!file.exists(f_ebd_nm)) {
  auk_filter(ebd_filters, file = f_ebd_nm, file_sampling = f_sampling_nm)
}


# Arizona -----------------------------------------------------------------

ebd_az <- auk_ebd(file = here('data',
                              '01_ebird_data',
                              'ebd_US-AZ_pinjay_smp_relJun-2024',
                              'ebd_US-AZ_pinjay_smp_relJun-2024.txt'),
                  file_sampling = here('data',
                                       '01_ebird_data',
                                       'ebd_US-AZ_pinjay_smp_relJun-2024',
                                       'ebd_US-AZ_pinjay_smp_relJun-2024_sampling.txt'))

#define the filters I will want
ebd_filters <- ebd_az %>%
  #filter just the months of interest, I think
  #we will want june and july, but need to verify
  #auk_date(date = c("*-06-01", "*-07-31")) %>%
  #filtering without bbs will be the nesting season
  #"late Feb to early May"
  auk_date(date = c("*-02-15", "*-05-15")) %>%
  #restirct to the stationary and traveling protocols
  auk_protocol(protocol = c("Stationary", "Traveling")) %>%
  #less than 5 h
  auk_duration(duration = c(0, 300)) %>%
  # less than 5km distance traveled
  auk_distance(distance = c(0, 5)) %>%
  #less than 10 observers
  #not sure how to do this one...
  #only use complete checklists
  auk_complete()

f_ebd_az <- file.path(here('data',
                           '01_ebird_data',
                           'cleaned_data',
                           '01_auk_filtering',
                           'ebd_az_filtered.txt'))
f_sampling_az <- file.path(here('data',
                                '01_ebird_data',
                                'cleaned_data',
                                '01_auk_filtering',
                                'ebd_az_checklists.txt'))

# only run if the files don't already exist
if (!file.exists(f_ebd_az)) {
  auk_filter(ebd_filters, file = f_ebd_az, file_sampling = f_sampling_az)
}


# UTah --------------------------------------------------------------------

ebd_ut <- auk_ebd(file = here('data',
                              '01_ebird_data',
                              'ebd_US-UT_pinjay_smp_relJun-2024',
                              'ebd_US-UT_pinjay_smp_relJun-2024.txt'),
                  file_sampling = here('data',
                                       '01_ebird_data',
                                       'ebd_US-UT_pinjay_smp_relJun-2024',
                                       'ebd_US-UT_pinjay_smp_relJun-2024_sampling.txt'))

#define the filters I will want
ebd_filters <- ebd_ut %>%
  #filter just the months of interest, I think
  #we will want june and july, but need to verify
  #auk_date(date = c("*-06-01", "*-07-31")) %>%
  #filtering without bbs will be the nesting season
  #"late Feb to early May"
  auk_date(date = c("*-02-15", "*-05-15")) %>%
  #restirct to the stationary and traveling protocols
  auk_protocol(protocol = c("Stationary", "Traveling")) %>%
  #less than 5 h
  auk_duration(duration = c(0, 300)) %>%
  # less than 5km distance traveled
  auk_distance(distance = c(0, 5)) %>%
  #less than 10 observers
  #not sure how to do this one...
  #only use complete checklists
  auk_complete()

f_ebd_ut <- file.path(here('data',
                           '01_ebird_data',
                           'cleaned_data',
                           '01_auk_filtering',
                           'ebd_ut_filtered.txt'))
f_sampling_ut <- file.path(here('data',
                                '01_ebird_data',
                                'cleaned_data',
                                '01_auk_filtering',
                                'ebd_ut_checklists.txt'))

# only run if the files don't already exist
if (!file.exists(f_ebd_ut)) {
  auk_filter(ebd_filters, file = f_ebd_ut, file_sampling = f_sampling_ut)
}

# Combine and zero-fill ---------------------------------------------------

# f_ebd_co <- read.delim(file = here('data',
#                                    'ebird_data',
#                                    'cleaned_data',
#                                    'ebd_co_filtered.txt'))
# f_ebd_co_check <- read.delim(file = here('data',
#                                          'ebird_data',
#                                          'cleaned_data',
#                                          'ebd_co_checklists.txt'))

ebd_co_zf <- auk_zerofill(f_ebd_co, f_sampling_co, collapse = TRUE)
ebd_nm_zf <- auk_zerofill(f_ebd_nm, f_sampling_nm, collapse = TRUE)
ebd_az_zf <- auk_zerofill(f_ebd_az, f_sampling_az, collapse = TRUE)
ebd_ut_zf <- auk_zerofill(f_ebd_ut, f_sampling_ut, collapse = TRUE)

ebd_all <- ebd_co_zf %>%
  rbind(ebd_nm_zf,
        ebd_az_zf,
        ebd_ut_zf)


# Extra filtering/cleaining -----------------------------------------------

# function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

# clean up variables
ebd_zf <- ebd_all %>% 
  mutate(
    # convert X to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )

# additional filtering
ebd_zf_filtered <- ebd_zf %>% 
  filter(
    # effort filters
    duration_minutes <= 5 * 60,
    effort_distance_km <= 5,
    # 10 or fewer observers
    number_observers <= 10)


# Select only variables of interest and expoort ---------------------------

ebird <- ebd_zf_filtered %>% 
  select(checklist_id, observer_id, sampling_event_identifier,
         scientific_name,
         observation_count, species_observed, 
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, year, day_of_year,
         time_observations_started, 
         duration_minutes, effort_distance_km,
         number_observers)

write.csv(ebird, here('data',
               '01_ebird_data',
               'cleaned_data',
               '02_all_auk_filtered',
               'all_ebird_data.csv'))


# Explore for ebird eyars -------------------------------------------------

ebdyrs <- ebird %>%
  filter(year > 2001) %>%
  group_by(year, species_observed) %>%
  tally() %>%
  pivot_wider(names_from = species_observed,
              values_from = n) %>%
  rowwise() %>%
  mutate(total = `TRUE` + `FALSE`) %>%
  ungroup() 

#looks like 2009 onward have more checklists
#with pinyon jays observed - so this could be a
#a good first year to filter based on??
