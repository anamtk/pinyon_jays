#Model summary figures
#March 3, 2025
#Ana Miller-ter Kuile

#exploreing results of the Ebird model without lags

# Load packages -----------------------------------------------------------

package.list <- c('tidyverse',
                  'here', 'patchwork', 'coda', 
                  "FNN", 'ggtext',
                  'grid') 

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())
# Load model summary ------------------------------------------------------
# 
betas <- readRDS(here('monsoon',
                      'outputs',
                      'ebird_abund_model2_summary.RDS'))
# 


beta_samples <- readRDS(here('monsoon',
                           'outputs',
                           'ebird_abund_model_covariate_effect_samples.RDS'))

beta_samps <- bind_rows(as.data.frame(beta_samples[[1]]),
                          as.data.frame(beta_samples[[2]]),
                          as.data.frame(beta_samples[[3]])) %>%
  mutate(sample = 1:n()) %>%
  pivot_longer(-sample,
               names_to = "parm",
               values_to = "value")

data_list <- readRDS(here('data',
                          '03_jags_input_data',
                          'ebird_data_list_nospuncert.RDS'))

ebird_blobIDs <- read.csv(here('data',
                               '01_ebird_data',
                               'cleaned_data',
                               '03_subsampled',
                               'ebird_cellIDlists.csv')) %>%
  dplyr::select(-X) %>%
  rename(blobID = cellID) %>%
  rename(cellID = cell)
#all blobs in blob lists to be able to get
#indexing
all_blobs <- ebird_blobIDs %>%
  distinct(year, blobnum, area) %>%
  group_by(year) %>%
  mutate(numID = 1:n()) %>%
  ungroup() %>%
  filter(year > 2009 & year < 2022)

cones <- read.csv(here('data',
                       '02_spatial_data',
                       'cleaned_data',
                       '02_weighted_blob_covariates',
                       'cones_weighted_mean_blob.csv')) 

temp <- read.csv(here('data',
                      '02_spatial_data',
                      'cleaned_data',
                      '02_weighted_blob_covariates',
                      'temp_weighted_mean_blob.csv'))

ppt <- read.csv(here('data',
                     '02_spatial_data',
                     'cleaned_data',
                     '02_weighted_blob_covariates',
                     'ppt_weighted_mean_blob.csv'))

monsoon <- read.csv(here('data',
                         '02_spatial_data',
                         'cleaned_data',
                         '02_weighted_blob_covariates',
                         'monsoon_weighted_mean_blob.csv')) %>%
  filter(blobnum %in% all_blobs$blobnum) 

pinyon <- read.csv(here('data',
                        '02_spatial_data',
                        'cleaned_data',
                        '02_weighted_blob_covariates',
                        'pinyonBA_weighted_mean_blob.csv')) %>%
  filter(blobnum %in% all_blobs$blobnum)

# Filter and scale all covariate datasets ---------------------------------

cones2 <- cones %>%
  filter(blobnum %in% all_blobs$blobnum) %>%
  mutate(cones = scale(cones)) %>%
  mutate(yearID = as.numeric(as.factor(year))) %>%
  left_join(all_blobs, by = c("year", "blobnum")) %>%
  dplyr::select(blobnum, cones, lag) %>%
  arrange(lag) %>%
  pivot_wider(names_from = "lag",
              values_from = "cones") %>%
  column_to_rownames(var= "blobnum") %>%
  as.matrix()

#tmax
temp2 <- temp %>%
  filter(blobnum %in% all_blobs$blobnum) %>%
  #trying scaling by season
  group_by(season) %>%
  mutate(temp = scale(temp)) %>%
  ungroup() %>%
  mutate(yearID = as.numeric(as.factor(year))) %>%
  left_join(all_blobs, by = c("year", "blobnum"))%>%
  dplyr::select(blobnum, temp, lag) %>%
  arrange(lag) %>%
  pivot_wider(names_from = "lag",
              values_from = "temp") %>%
  column_to_rownames(var= "blobnum") %>%
  as.matrix()

#ppt
ppt2 <- ppt %>%
  filter(blobnum %in% all_blobs$blobnum) %>%
  #trying scaling by season
  group_by(season) %>%
  mutate(ppt = scale(ppt)) %>%
  ungroup() %>%
  mutate(yearID = as.numeric(as.factor(year))) %>%
  left_join(all_blobs, by = c("year", "blobnum"))%>%
  dplyr::select(blobnum, ppt, lag) %>%
  arrange(lag) %>%
  pivot_wider(names_from = "lag",
              values_from = "ppt") %>%
  column_to_rownames(var= "blobnum") %>%
  as.matrix()


# Correction function per sample ------------------------------------------

correct_wt_fun <- function(sample){
  df <- beta_samps %>%
    filter(sample == {{sample}})
  
  #cone median weights
  conewt <- df %>%
    filter(str_detect(parm, "wA")) %>%
    dplyr::select(value) %>%
    as_vector()
  
  coneVals <- apply(cones2, MARGIN = 1, FUN = function(x){sum(x*conewt, na.rm = T)})
  
  #tmax
  tmxwt <- df %>%
    filter(str_detect(parm, "wB")) %>%
    dplyr::select(value) %>%
    as_vector()
  
  tempVals <- apply(temp2, MARGIN = 1, FUN = function(x){sum(x*tmxwt, na.rm = T)})
  
  #ppt
  pptwt <- df %>%
    filter(str_detect(parm, "wC")) %>%
    dplyr::select(value) %>%
    as_vector()
  
  pptVals <- apply(ppt2, MARGIN = 1, FUN = function(x){sum(x*pptwt, na.rm = T)})
  
  covariate <- c('a0', "Cones", "Tmax", "PPT", "Monsoon",
                 "PinyonBA", "ConexTmax", "ConexPPT",
                 'ConexMonsoon', "ConexBA")
  
  beta_df <- df %>%
    filter(str_detect(parm, "a")) %>%
    filter(!parm %in% c("deviance")) %>%
    filter(!parm %in% c("wA")) %>%
    bind_cols(covariate = covariate) %>%
    mutate(value_corrected = case_when(covariate == "Cones" ~ value*sd(coneVals),
                                       covariate == "Tmax" ~ value*sd(tempVals),
                                       covariate == "	PPT" ~ value*sd(pptVals),
                                       covariate == "ConexTmax" ~ value*sd(coneVals)*sd(tempVals),
                                       covariate == "	ConexPPT" ~ value*sd(coneVals)*sd(pptVals),
                                       TRUE ~ value)) %>%
    mutate(sample = {{sample}}) 
  
  return(beta_df)
  
}


# Run on all samples ------------------------------------------------------
sample_vect <- 1:max(beta_samps$sample)

beta_list <- lapply(sample_vect, correct_wt_fun)

corr_beta_df <- as.data.frame(do.call(rbind, beta_list))

corr_beta_sum <- corr_beta_df %>%
  group_by(covariate) %>%
  summarise(median = median(value_corrected),
         lci = quantile(value_corrected, probs = c(0.025),
                        type = 8),
         uci = quantile(value_corrected, probs = c(0.975),
                        type = 8)) %>%
  mutate(type = case_when(covariate %in% c("Cones", "Monsoon",
                                           "PPT", "PinyonBA",
                                           "Tmax") ~ "Main effects",
                          covariate %in% c("ConexBA",
                                           "ConexMonsoon",
                                           "ConexPPT",
                                           "ConexTmax") ~ "Interactions",
                          TRUE ~ NA_character_)) %>%
  mutate(covariate = case_when(covariate == "Monsoon" ~ "Monsoonality",
                               covariate == "Cones" ~ "Pinyon (seed) cone availability",
                               covariate == "PinyonBA" ~ "Pinyon basal area (BA)",
                               covariate == "PPT" ~ "Precipitation (PPT)",
                               covariate == "Tmax" ~ "Maximum temperature (Tmax)",
                               covariate == "ConexBA" ~ "Cones x BA",
                               covariate == "ConexMonsoon" ~ "Cones x monsoonality",
                               covariate == "ConexPPT" ~ "Cones x PPT",
                               covariate == "ConexTmax" ~ "Cones x Tmax",
                               TRUE ~ covariate))
  

saveRDS(corr_beta_sum, here('data',
                            '04_jags_output_data',
                            'corrected_covariate_effects.RDS'))
# Figure 2: covariate effects ---------------------------------------------


(covariate_betas <- corr_beta_sum %>%
  filter(covariate != 'a0') %>%
   mutate(type = factor(type, levels = c("Main effects", "Interactions"))) %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_pointrange(aes(x = median, 
                      y = reorder(covariate, median),
                      xmin = lci,
                      xmax = uci),
                  size = 0.25) + 
  labs(x = "Covariate effect\n(median and 95% CI)",
       y = "Covariate") +
   facet_grid(type~., scales ="free") +
   theme(strip.background = element_rect(fill = "white",
                                         color = "white"),
         strip.text = element_text(size = 12),
         axis.title = element_text(size = 15),
         axis.text = element_text(size = 10))
)

ggsave(plot = covariate_betas,
       here('pictures',
                   'final',
                   'covariate_effects.jpg'),
       width = 6.25,
       height = 4,
       dpi = 300,
       units = "in")

