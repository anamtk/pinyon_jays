#Look at distribution of in vs. all of all covariates
#Ana Miller-ter Kuile
#June 26, 2025

#this script looks at the distribution of covariate
#data both with ebird checklists and without ebird checklists to
#determine whether ebird sampled in a relatively unbiased way
#along the range of the covariates

#what we'd want to see is that there are checklists across the range
#of covariate conditions, not whether those distributions are ~the same
#since we've already subsampled checklists in some ways that bias them
#towards being in locations where we have cone/pinyon data

#this is a requirement for being able to use this kind of semi-structured
#data in occupancy-type models with a "space-for-time" substitution, 
#treating all repeat observations in a given space as "replicates" or 
#"repeat surveys" like in a normal occupancy model.

#e.g., this paper: https://academic.oup.com/auk/article/140/4/ukad035/7229681

#summary:
# - good job spanning climate variables
# - missing "really high" values for both cones and pinyon basal area
##   - these "really high" values are minority of the total available covariate values
##     (most checklist-nochecklist values are in the similar part of the distribution
##      of both covariates) 

# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse", 
                  "sf",  
                  "terra",
                  'readxl',
                  'sf',
                  'patchwork')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())
# Load data ---------------------------------------------------------------

ebird <- read.csv(here('data',
                       '01_ebird_data',
                       'cleaned_data',
                       '03_subsampled',
                       'ebird_cellIDlists.csv')) %>%
  dplyr::select(-X) %>%
  rename(blobID = cellID) %>%
  rename(cellID = cell) 



# Join climate data with ebird data ---------------------------------------

# Cones -------------------------------------------------------------------

cones <- read.csv(here('data',
                       '02_spatial_data',
                       'cleaned_data',
                       '01_get_gridIDs',
                       'cone_masting_df.csv'))%>%
  rename(cellID = cell)

cones2 <- cones %>%
  left_join(ebird, by = "cellID") %>%
  dplyr::select(-(X1897:X2007)) %>%
  pivot_longer(X2008:X2024,
               names_to = 'cone_yr',
               values_to = "cones") %>%
  mutate(type = case_when(!is.na(blobnum) ~ "checklist",
                          TRUE ~ "no_checklist")) %>%
  mutate(logcones = log(cones + 0.000001)) 

cone_plot <- ggplot(cones2, aes(x = type, y = logcones)) +
  geom_violin() +
  labs(x = "Checklist present in grid",
       y = "Cone abundance") +
  scale_x_discrete(labels = c("Yes", "No"))

#maybe don't have good coverage in the realllly high range - which 
#i think is generally super rare anyway? so probably hard for anyone to
#find IRL honestly
max_checklist_cones <- cones2 %>%
  filter(type == "checklist") %>%
  filter(cones == max(cones, na.rm = T)) %>%
  distinct(cones) %>%
  as_vector()

nocheck_num <- cones2 %>%
  filter(type == "no_checklist")

nocheck_cones <- cones2 %>%
  filter(cones > max_checklist_cones)

nrow(nocheck_cones)/nrow(nocheck_num)

#0.0000235
#or, 0.00235% of cells - SO FEW

ggsave(here("pictures",
            "final",
            'si',
            "space_for_time_cones.jpg"),
       width = 3,
       height = 2,
       units = 'in',
       dpi = 300
)

rm(cones, cones2, cone_plot)
# temp --------------------------------------------------------------------

temp <- read.csv(here('data',
                      '02_spatial_data',
                      'cleaned_data',
                      '01_get_gridIDs',
                      'temp_data_df.csv')) %>%
  dplyr::select(cellID, 
                PRISM_tmax_stable_4kmM3_200001_bil:PRISM_tmax_stable_4kmM3_202312_bil) 


temp2 <- temp %>%
  left_join(ebird, by = c("cellID")) %>%
  pivot_longer(PRISM_tmax_stable_4kmM3_200701_bil:PRISM_tmax_stable_4kmM3_202212_bil,
               names_to = 'temp_yrmonth',
               values_to = "temp") %>%
mutate(type = case_when(!is.na(blobnum) ~ "checklist",
                        TRUE ~ "no_checklist")) 

temp_plot <- ggplot(temp2, aes(x = type, y = temp)) +
  geom_violin()  +
  labs(x = "Checklist present in grid",
       y = "Maximum temperature") +
  scale_x_discrete(labels = c("Yes", "No"))

#this one looks good!
ggsave(here("pictures",
            "final",
            'si',
            "space_for_time_temp.jpg"),
       width = 3,
       height = 2,
       units = 'in',
       dpi = 300
)

rm(temp2, temp, temp_plot)
# PPT ---------------------------------------------------------------------

ppt <- read.csv(here('data',
                     '02_spatial_data',
                     'cleaned_data',
                     '01_get_gridIDs',
                     'ppt_data_df.csv')) %>%
  dplyr::select(cellID, 
                PRISM_ppt_stable_4kmM3_200001_bil:PRISM_ppt_stable_4kmM3_202312_bil)

ppt2 <- ppt %>%
  left_join(ebird, by = c("cellID")) %>%
  pivot_longer(PRISM_ppt_stable_4kmM3_200701_bil:PRISM_ppt_stable_4kmM3_202212_bil,
               names_to = 'ppt_yrmonth',
               values_to = "ppt") %>%
  mutate(type = case_when(!is.na(blobnum) ~ "checklist",
                          TRUE ~ "no_checklist"))

ppt2 <- ppt2 %>% 
  group_by(type) %>%
  slice_sample(n = 100000)

(ppt_plot <- ggplot(ppt2, aes(x = type, y = ppt)) +
  geom_violin() +
  labs(x = "Checklist present in grid",
       y = "Precipitation") +
  scale_x_discrete(labels = c("Yes", "No")))

#maybe not so good at getting surveys during the really wet periods - but also
#seems fine honestly

ggsave(here("pictures",
            "final",
            'si',
            "space_for_time_ppt.jpg"),
       width = 3,
       height = 2,
       units = 'in',
       dpi = 300
)


rm(ppt, ppt2, ppt_plot)

# Monsoon -----------------------------------------------------------------
monsoon <- read.csv(here('data',
                         '02_spatial_data',
                         'cleaned_data',
                         '01_get_gridIDs',
                         'monsoon_data_df.csv')) %>%
  dplyr::select(cellID, PRISM_ppt_30yr_normal_800mM4_07_bil) %>%
  rename(monsoon = PRISM_ppt_30yr_normal_800mM4_07_bil)

monsoon2 <- monsoon %>%
  left_join(ebird, by = "cellID") %>%
  mutate(type = case_when(!is.na(blobnum) ~ "checklist",
                          TRUE ~ "no_checklist")) 
  
monsoon_plot <- ggplot(monsoon2, aes(x = type, y = monsoon)) +
  geom_violin() +
  labs(x = "Checklist present in grid",
       y = "Monsoonality") +
  scale_x_discrete(labels = c("Yes", "No"))

#looks ok to me
ggsave(here("pictures",
            "final",
            'si',
            "space_for_time_monsoon.jpg"),
       width = 3,
       height = 2,
       units = 'in',
       dpi = 300
)
# Pinyon BA ---------------------------------------------------------------

pinyon <- read.csv(here('data',
                        '02_spatial_data',
                        'cleaned_data',
                        '01_get_gridIDs',
                        'pinyonba_data_df.csv'))%>%
  dplyr::select(cell, 
                PinyonBA_sqftPerAc_2010:PinyonBA_sqftPerAc_2022) %>%
  rename(cellID = cell) 

pinyon2 <- pinyon %>%
  left_join(ebird, by = "cellID") %>%
  pivot_longer(PinyonBA_sqftPerAc_2010:PinyonBA_sqftPerAc_2022,
               names_to = 'pinyon_yr',
               values_to = "pinyon") %>%
  mutate(type = case_when(!is.na(blobnum) ~ "checklist",
                          TRUE ~ "no_checklist")) 

pba_plot <- ggplot(pinyon2, aes(x = type, y = pinyon)) +
  geom_violin()+
  labs(x = "Checklist present in grid",
       y = "Pinyon basal area") +
  scale_x_discrete(labels = c("Yes", "No"))

ggsave(here("pictures",
            "final",
            'si',
            "space_for_time_pba.jpg"),
       width = 3,
       height = 2,
       units = 'in',
       dpi = 300
)

#i think this high BA missing value is correlated with the missing 
#values of cone production, i bet. it seems like it is *so little* 
#of the distribution - it doesn't seem like birders are going out 
#and finding "optimal" pinyon jay habitat to sample
  
max_checklist_ba <- pinyon2 %>%
  filter(type == "checklist") %>%
  filter(pinyon == max(pinyon, na.rm = T)) %>%
  distinct(pinyon) %>%
  as_vector()

nocheck_num <- pinyon2 %>%
  filter(type == "no_checklist")

nocheck_pinyon <- pinyon2 %>%
  filter(pinyon > max_checklist_ba)

nrow(nocheck_pinyon)/nrow(nocheck_num)
#222 grid cells in this category versus ~800,000 available grid cells
#0.026% of total possible pinyon basal area values in this edge of the range
  
# cone_plot + ppt_plot + temp_plot + monsoon_plot + pba_plot +
#   plot_annotation(tag_levels = "a")
# 
# ggsave(here("pictures",
#             "final",
#             'si',
#             "space_for_time.jpg"),
#        width = 7,
#        height = 4,
#        units = 'in',
#        dpi = 300
#        )
