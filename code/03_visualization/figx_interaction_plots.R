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

monsoon2 <- monsoon %>%
  mutate(wt= scale(wt))

pinyon2 <- pinyon %>%
  mutate(wt = scale(wt))

# Interaction plots -------------------------------------------------------


# WEighted covariate function ---------------------------------------------

weighted_cov_fun <- function(sample){
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
  
  all_cov_df <- as.data.frame(bind_cols(cones = coneVals,
                                        ppt = pptVals,
                                        tmax = tempVals,
                                        monsoon = monsoon2$wt,
                                        pinyonba = pinyon2$wt)) %>%
    mutate(sample = {{sample}}) %>%
    bind_cols(all_blobs)
  
  return(all_cov_df)
  
}

sample_vect <- 1:max(beta_samps$sample)

weighted_cov_list <- lapply(sample_vect, weighted_cov_fun)

weighted_cov_df <- as.data.frame(do.call(rbind, weighted_cov_list))
#weighted climate actual data (at blobs)
#coneVals
#pptVals
#tempVals
#others:
#monsoon
#pinyon

weighted_cov_sum_df <- weighted_cov_df %>%
  dplyr::select(-area) %>%
  group_by(year, blobnum, numID) %>%
  summarise(cones = mean(cones, na.rm = T),
            ppt = mean(ppt, na.rm = T),
            tmax = mean(tmax,na.rm = T),
            monsoon = mean(monsoon, na.rm = T),
            pinyonba = mean(pinyonba, na.rm = T)) %>%
  ungroup() 

# Make simulated range of predicted data ----------------------------------

#need a continuous cones, ppt, and temp based on the range of the data
#and all combinations in a set of 2 DFs
cone_sim <- seq(from = min(weighted_cov_df$cones), to = max(weighted_cov_df$cones), length.out = 100)
ppt_sim <- seq(from = min(weighted_cov_df$ppt), to = max(weighted_cov_df$ppt), length.out = 100)
temp_sim <- seq(from = min(weighted_cov_df$tmax), to = max(weighted_cov_df$tmax), length.out = 100)
mons_sim <- seq(from = min(monsoon2$wt, na.rm = T), 
                to = max(monsoon2$wt, na.rm = T), 
                length.out = 100)
ba_sim <- seq(from = min(pinyon2$wt, na.rm = T), 
              to = max(pinyon2$wt),
              length.out= 100)


#IDs for covariates:  

covariate <- c('a0', "Cones", "Tmax", "PPT", "Monsoon",
               "PinyonBA", "ConexTmax", "ConexPPT",
               'ConexMonsoon', "ConexBA")

beta_median_df <- beta_samps %>%
  filter(str_detect(parm, "a")) %>%
  filter(!parm %in% c("deviance")) %>%
  filter(!parm %in% c("wA")) %>%
  mutate(covariate = case_when(parm == "a0" ~ "a0",
                               parm == "a[1]" ~ 'Cones',
                               parm == "a[2]" ~ "Tmax",
                               parm == "a[3]" ~ 'PPT',
                               parm == "a[4]" ~ 'Monsoon',
                               parm == "a[5]" ~ 'PinyonBA',
                               parm == "a[6]" ~ "ConesxTmax",
                               parm == 'a[7]' ~ "ConesxPPT",
                               parm == "a[8]" ~ "ConesxMonsoon",
                               parm == "a[9]" ~ "ConesxBA",
                               TRUE ~ NA_character_)) %>%
  group_by(covariate) %>%
  summarise(median = median(value)) %>%
  ungroup()

cov_function <- function(cov){
  vec <- beta_median_df %>%
    filter(covariate == {{cov}}) %>%
    dplyr::select(median) %>%
    as_vector() 
  
  return(vec)
  
}


#a0 
a0 <- cov_function(cov = "a0")
#a[1]:cones
aCones <- cov_function(cov = "Cones")
#a[2]:tmax
aTmax <- cov_function(cov =  "Tmax")
#a[3]:ppt
aPPT <- cov_function(cov = "PPT")
#a[4]:monsoon
aMons <- cov_function(cov = "Monsoon")
#a[5]:ba
aBA <- cov_function(cov = "PinyonBA")
#a[6]: conextmax
aConeTmax <- cov_function(cov = "ConesxTmax")
#a[7]: conexppt
aConePPT <- cov_function(cov = "ConesxPPT")
#a[8]: conexmonsoon
aConeMons <- cov_function(cov = "ConesxMonsoon")
#a[9]:conexba
aConeBA <- cov_function(cov = "ConesxBA")


# Get scaled-unscaled values for plotting axes ----------------------------

cones_mean <- mean(cones$cones, na.rm =T)
cones_sd <- sd(cones$cones, na.rm = T)
mons_mean <- mean(monsoon$wt, na.rm = T)
mons_sd <- sd(monsoon$wt, na.rm = T)
tmax_mean <- mean(temp$temp, na.rm = T)
tmax_sd <- sd(temp$temp, na.rm = T)
ppt_mean <- mean(ppt$ppt, na.rm = T)
ppt_sd <- sd(ppt$ppt, na.rm = T)
ba_mean <- mean(pinyon$wt, na.rm = T)
ba_sd <- sd(pinyon$wt, na.rm = T)

# Figure 4: interaction plots ---------------------------------------------

weighted_cov_sum_df2 <- weighted_cov_sum_df %>%
  sample_n(1500)

cone_interaction_function <- function(int_cov,
                                      cov_name,
                                      beta,
                                      int_beta,
                                      cov_mean,
                                      cov_sd){
  
  int_df <- expand.grid(cone_sim, int_cov) %>%
    rename(cones = Var1,
           cov = Var2) %>%
    mutate(loglambda = a0 + aCones*cones +
             beta*cov + int_beta*cones*cov,
           lambda = exp(loglambda))
  
  cov_vector <- weighted_cov_sum_df %>%
    dplyr::select({{cov_name}})
  
  int_df$nnDists <- knnx.dist(data = cbind(weighted_cov_sum_df$cones, 
                                           cov_vector),
                              query = cbind(int_df$cones, int_df$cov),
                              k = 1)    
  
  # And set values far from observed data (> 0.1 distance in z-scores across both axes) to NA
  int_df$loglambda[int_df$nnDists > 0.75] <- NA  
  
 # lab <- expression(paste(italic("log"), "(", lambda, ")"))
 
  #lab <- expression(paste("birds" %.% m^{-2}))
  
  lab <- expression(paste(birds %.% ha^{-1}))
  
  unscaled_x_labels <- function(breaks) {
    paste(round((breaks*cones_sd) + cones_mean, digits = 1))
  }
  
  unscaled_y_labels <- function(breaks) {
    paste(round((breaks*cov_sd) + cov_mean, digits = 1))
  }
  
  plot <- ggplot(int_df) +
    geom_tile(aes(x = cones, y = cov, fill = loglambda)) +
    scale_fill_distiller(type = "seq",
                         palette = "PuRd",
                         direction = 1,
                         breaks = c(-9.2, -13.8, -18.4),
                         labels = c(1, 0.01, 0.0001), 
                         na.value = "transparent") +
    geom_contour(aes(x = cones, y = cov, z = loglambda), color = "lightgrey", alpha = 0.5) +
    geom_point(data = weighted_cov_sum_df2, 
               aes(x = cones, y = {{cov_name}}), 
               color = "black", 
               alpha = 0.2, shape = 1) +
    scale_x_continuous(labels = unscaled_x_labels) +
    scale_y_continuous(labels = unscaled_y_labels) +
    labs(fill = lab) +
    theme(axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.text = element_text(size = 10))
  
  return(plot)
  
}

# scale_fill_distiller(type = "seq",
#                      palette = "PuRd",
#                      direction = 1,
#                      breaks = c(-9, -14, -20),
#                      labels = c(1e-04, 1e-06, 1e-09)) 
#meters

#conextmax:
(a <- cone_interaction_function(int_cov = temp_sim,
                               beta = aTmax,
                               int_beta = aConeTmax,
                               cov_name = tmax,
                               cov_mean = tmax_mean,
                               cov_sd = tmax_sd) +
  labs(x = "Cones", y = "Maximum temperature (\u00B0C)") +
    theme(axis.title.x = element_blank()))

#conexppt:
(b <- cone_interaction_function(int_cov = ppt_sim,
                               beta = aPPT,
                               int_beta = aConePPT,
                               cov_name = ppt,
                               cov_mean = ppt_mean,
                               cov_sd = ppt_sd)+
  labs(x = "Cones", y = "Precipitation (mm)")+
  theme(legend.position = "none",
        axis.title.x = element_blank()))

#conexmonsoon
(c <- cone_interaction_function(int_cov = mons_sim,
                               beta = aMons,
                               int_beta = aConeMons,
                               cov_name = monsoon,
                               cov_mean = mons_mean,
                               cov_sd = mons_sd)+
  labs(x = "Cones", y = "Monsoonality (%)")+
  theme(legend.position = "none"))

#conexba
balab <- expression(paste('Pinyon basal area (', m^2 %.% ha^-1,")" ))

(d <- cone_interaction_function(int_cov = ba_sim,
                               beta = aBA,
                               int_beta = aConeBA,
                               cov_name = pinyonba,
                               cov_mean = ba_mean,
                               cov_sd = ba_sd)+
  labs(x = "Cones", y = balab) +
  theme(legend.position = "none"))

b + a+ c+ d +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a",
                  tag_suffix = ")")

ggsave(here('pictures',
            'final',
            'interaction_plots.jpg'),
       width = 6.5,
       height = 5,
       units = 'in',
       dpi = 300)



