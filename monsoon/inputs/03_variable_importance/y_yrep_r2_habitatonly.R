
# Load packages -----------------------------------------------------------

package.list <- c("jagsUI", "coda",
                  'dplyr', 'stringr',
                  'magrittr', 'tidyr',
                  'mcmcplots','ggplot2')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# -------------------------------------------------------------------------

data <- readRDS("/scratch/atm234/pinyon_jays/ebird/nospuncert/inputs/ebird_check_blob_yr_ids.RDS")

data <- data %>%
  dplyr::select(numID, yrID, obsrvtn_c, checkID)

yrep_samps <- readRDS("/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model_yrepsamples_habitatonly.RDS")

yrep_samps_df <- bind_rows(as.data.frame(yrep_samps[[1]]),
                           as.data.frame(yrep_samps[[2]]),
                           as.data.frame(yrep_samps[[3]])) %>%
  mutate(sample = 1:n()) %>%
  pivot_longer(-sample,
               names_to = "parm",
               values_to = "count") %>%
  filter(parm != "deviance") 

yrep_samps_df2 <- yrep_samps_df %>%
  separate(parm, 
           into = c("yrID", "numID", "checkID"),
           sep = ",") %>%
  mutate(yrID = str_sub(yrID, 17, nchar(yrID)),
         checkID =str_sub(checkID, 1, (nchar(checkID)-1)),
         yrID =as.numeric(yrID),
         numID = as.numeric(numID),
         checkID = as.numeric(checkID))

# R2 from samples ---------------------------------------------------------

y_yrep_samps_df <- yrep_samps_df2 %>%
  left_join(data, by = c("yrID", "numID", "checkID"))

y_yrep_r2 <- function(sample){
  
  yrep <- y_yrep_samps_df %>%
    filter(sample == {{sample}}) 
  
  #combine that yrep with observed y
  
  R2 <- cor(log(yrep$obsrvtn_c+.01), 
            log(yrep$count+.01),
      use = 'complete.obs')^2
  
  R2_nolog <- cor(yrep$obsrvtn_c, yrep$count,
                  use = 'complete.obs')^2
  
  R2_df <- as.data.frame(cbind(R2 = R2,
                               R2_nolog = R2_nolog))
  
  return(list(R2_df))
  
}

samps <- unique(y_yrep_samps_df$sample)

r2_list <- lapply(samps, y_yrep_r2)

r2_df <- as.data.frame(do.call(rbind, r2_list))

saveRDS(r2_df, "/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model_yyrepr2_habitatonly.RDS")

