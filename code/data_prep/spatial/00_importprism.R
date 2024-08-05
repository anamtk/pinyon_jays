library(prism)
library(terra)
library(here)
library(tidyverse)

prism_set_dl_dir(path =
                   here('data',
                        'spatial_data',
                        'prism_temp'), create = TRUE)

get_prism_annual(type = "tmean", year = 2000:2023, keepZip = FALSE)

prism_set_dl_dir(path =
                   here('data',
                        'spatial_data',
                        'prism_ppt'), create = TRUE)

get_prism_annual(type = "ppt", year = 2000:2023, keepZip = FALSE)

prism_set_dl_dir(path =
                   here('data',
                        'spatial_data',
                        'prism_vpd'), create = TRUE)

get_prism_annual(type = "vpdmax", year = 2000:2023, keepZip = FALSE)

