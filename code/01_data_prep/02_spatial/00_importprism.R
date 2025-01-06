library(prism)
library(terra)
library(here)
library(tidyverse)

prism_set_dl_dir(path =
                   here('data',
                        'spatial_data',
                        'prism_monthly_tmax'), create = TRUE)

get_prism_monthlys(type = "tmax", year = 2000:2023, keepZip = FALSE)

prism_set_dl_dir(path =
                   here('data',
                        'spatial_data',
                        'prism_monthly_ppt'), create = TRUE)

get_prism_monthlys(type = "ppt", year = 2000:2023, keepZip = FALSE)

# prism_set_dl_dir(path =
#                    here('data',
#                         'spatial_data',
#                         'prism_vpd'), create = TRUE)
# 
# get_prism_annual(type = "vpdmax", year = 2000:2023, keepZip = FALSE)
# 
# prism_set_dl_dir(path =
#                    here('data',
#                         'spatial_data',
#                         'prism_monthly_vpd'), create = TRUE)
# 
# get_prism_monthlys(type = "vpdmax", years = 2000:2023, keepZip = FALSE)

prism_set_dl_dir(path =
                   here('data',
                        'spatial_data',
                        'prism_monthly_tmean'), create = TRUE)

get_prism_monthlys(type = "tmean", year = 2000:2023, keepZip = FALSE)
