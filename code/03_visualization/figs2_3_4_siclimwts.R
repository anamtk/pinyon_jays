#Model summary figures
#March 3, 2025
#Ana Miller-ter Kuile

#exploreing results of the Ebird model without lags

# Load packages -----------------------------------------------------------

package.list <- c('tidyverse',
                  'here', 'patchwork', 'coda', "FNN") 

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())
# Load model summary ------------------------------------------------------
# 
betas <- readRDS(here('monsoon',
                      'ebird',
                      'nospuncert',
                      'outputs',
                      'ebird_abund_model2_summary.RDS'))
# 


beta_samples <- readRDS(here('monsoon',
                           'ebird',
                           'nospuncert',
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
                       'cones_weighted_mean_blob.csv')) 

temp <- read.csv(here('data',
                      '02_spatial_data',
                      'cleaned_data',
                      'temp_weighted_mean_blob.csv'))

ppt <- read.csv(here('data',
                     '02_spatial_data',
                     'cleaned_data',
                     'ppt_weighted_mean_blob.csv'))

monsoon <- read.csv(here('data',
                         '02_spatial_data',
                         'cleaned_data',
                         'monsoon_weighted_mean_blob.csv')) %>%
  filter(blobnum %in% all_blobs$blobnum) 

pinyon <- read.csv(here('data',
                        '02_spatial_data',
                        'cleaned_data',
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
  mutate(covariate = case_when(covariate == "Monsoon" ~ "Monsoonality \n (Monsoon)",
                               covariate == "Cones" ~ "Pinyon pine seed \n availability (Cones)",
                               covariate == "PinyonBA" ~ "Pinyon basal area \n (PBA)",
                               covariate == "PPT" ~ "Precipitation \n (PPT)",
                               covariate == "Tmax" ~ "Maximum temperature \n (Tmax)",
                               covariate == "ConexBA" ~ "Cones by PBA",
                               covariate == "ConexMonsoon" ~ "Cones by Monsoon",
                               covariate == "ConexPPT" ~ "Cones by PPT",
                               covariate == "ConexTmax" ~ "Cones by Tmax",
                               TRUE ~ covariate))
  


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
                  size = 0.05) + 
  labs(x = "Covariate effect\n(median and 95% BCI)",
       y = "Covariate") +
   facet_grid(type~., scales ="free") +
   theme(strip.background = element_rect(fill = "white",
                                         color = "white"))
)

ggsave(plot = covariate_betas,
       here('pictures',
                   'final',
                   'covariate_effects.jpg'),
       width = 5,
       height = 3.5,
       dpi = 300,
       units = "in")



# Weights -----------------------------------------------------------------

weights_df <- as.data.frame(betas$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "w")) %>%
  mutate(covariate = case_when(str_detect(parm, "wA") ~ "Cones",
                               str_detect(parm, "wB") ~ "Temperature",
                               str_detect(parm, 'wC') ~ "Precipitation")) %>%
  mutate(lag = str_sub(parm, 4, (nchar(parm)-1)),
         lag = as.numeric(lag)) 



# Figure 3: cone wts ------------------------------------------------------


start <- c(0.5,2.55, 3.55)
end <- c(2.45, 3.45, 5.5)
times <- c("1. Quicker resources ('predict')", 
           "2. Cone caching ('lag')", 
           "3. More breeders in population ('lag')")

box_df <- as.data.frame(start) %>%
  bind_cols(end = end, times = times) %>%
  mutate(times = factor(times, levels = c("1. Quicker resources ('predict')", 
                                          "2. Cone caching ('lag')", 
                                          "3. More breeders in population ('lag')")))

(cone_weight_plot <- weights_df %>%
  filter(covariate == "Cones") %>%
ggplot(aes(x = lag, y = `50%`)) +
  geom_rect(data = box_df, aes(xmin = start,
                               xmax = end,
                               group = times,
                               ymin = 0,
                               ymax = 1,
                               fill = times),
            inherit.aes = F, alpha = 0.6) +
  scale_fill_manual(values = c('#02818A',
                               '#67A9CF',
                               "#A6BDDB")) +
  #geom_hline(yintercept = 1/5, linetype = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`,
                    ymax = `97.5%`),
                width = 0.2) +
  scale_x_continuous(breaks = c(1:5),
                     labels = c( "2", "1", "1", 
                                "2", "3")) +
  annotate(geom = "text", x = 1.5, y = 0.9, label = "Years before cones") +
  annotate(geom = "text", x = 4.5, y = 0.9, label = "Years after cones") +
    #geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 1)) +
  labs(x = "",
       y = "Importance weight \n (median and 95% BCI)") )

ggsave(here('pictures',
            'R',
            'cone_weights.pdf'),
       width = 6,
       height = 3,
       units = 'in',
       dpi = 300)
#green - lag
#ccebc5
#blue - predictive
#b3cde3
#yellow - immediate
#ffffcc

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
                                        monsoon = monsoon$wt,
                                        pinyonba = pinyon$wt)) %>%
    mutate(sample = {{sample}}) %>%
    bind_cols(all_blobs)
  
  return(all_cov_df)
  
}

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
mons_sim <- seq(from = min(monsoon$wt, na.rm = T), 
                to = max(monsoon$wt, na.rm = T), length.out = 100)
ba_sim <- seq(from = min(pinyon$wt, na.rm = T), to = max(pinyon$wt),
              length.out= 100)


#IDs for covariates:  

covariate <- c('a0', "Cones", "Tmax", "PPT", "Monsoon",
               "PinyonBA", "ConexTmax", "ConexPPT",
               'ConexMonsoon', "ConexBA")

# a[1]*AntCone[t,i] +
#   a[2]*AntTmax[t,i] +
#   a[3]*AntPPT[t,i] +
#   a[4]*Monsoon[t,i] +
#   a[5]*PinyonBA[t,i] +
#   a[6]*AntCone[t,i]*AntTmax[t,i] + 
#   a[7]*AntCone[t,i]*AntPPT[t,i] + 
#   a[8]*AntCone[t,i]*Monsoon[t,i] + 
#   a[9]*AntCone[t,i]*PinyonBA[t,i]

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


# Figure 4: interaction plots ---------------------------------------------

cone_interaction_function <- function(int_cov,
                                      cov_name,
                                      beta,
                                      int_beta){
  
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
  
  #return(int_df)
  #c <- noquote('")')
  lab <- expression(paste(italic("log"), "(", lambda, ")"))
  #lab <- expression(paste("ln(E(x))"))
  #lab <- (expression(paste(lambda), "\n (birds \cdot m^{-2})"))
  
  #ylab(expression(Anthropogenic~SO[4]^{"2-"}~(ngm^-3))) 
  lab <- expression(paste("birds" %.% m^{-2}))
  
  plot <- ggplot(int_df) +
    geom_tile(aes(x = cones, y = cov, fill = loglambda)) +
    scale_fill_distiller(type = "seq",
                         palette = "PuRd",
                         direction = 1,
                         breaks = c(-9.2, -13.8, -18.4),
                         labels = c(1e-04, 1e-06, 1e-08), 
                         na.value = "transparent") +
    geom_contour(aes(x = cones, y = cov, z = loglambda), color = "lightgrey", alpha = 0.5) +
    geom_point(data = weighted_cov_sum_df, 
               aes(x = cones, y = {{cov_name}}), 
               color = "black", 
               alpha = 0.2, shape = 1) +
    labs(fill = lab) 
  
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
                               cov_name = tmax) +
  labs(x = "Cones", y = "Maximum temperature"))

#conexppt:
b <- cone_interaction_function(int_cov = ppt_sim,
                               beta = aPPT,
                               int_beta = aConePPT,
                               cov_name = ppt)+
  labs(x = "Cones", y = "Precipitation")+
  theme(legend.position = "none")

#conexmonsoon
(c <- cone_interaction_function(int_cov = mons_sim,
                               beta = aMons,
                               int_beta = aConeMons,
                               cov_name = monsoon)+
  labs(x = "Cones", y = "Monsoonality")+
  theme(legend.position = "none"))

#conexba
d <- cone_interaction_function(int_cov = ba_sim,
                               beta = aBA,
                               int_beta = aConeBA,
                               cov_name = pinyonba)+
  labs(x = "Cones", y = "Pinyon basal area") +
  theme(legend.position = "none")

b + a+ c+ d +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a",
                  tag_suffix = ")")

ggsave(here('pictures',
            'final',
            'interaction_plots.jpg'),
       width = 6,
       height = 5,
       units = 'in',
       dpi = 300)

# SI fig: climate weights -------------------------------------------------


#Season 1: 2-4: breeding
#Season 2: 5-6: feeding dependent young
#Season 3: 7: unknown, just chillin?
#Season 4: 8-1: winter foraging, potentially seeking food elsewhere

weights_df %>%
  filter(!covariate %in%  c("Cones", "Temperature")) %>%
  arrange(desc(`50%`))

weights_df %>%
  filter(!covariate %in%  c("Cones", "Precipitation")) %>%
  arrange(desc(`50%`))

(clim_weights_plot <- weights_df %>%
  filter(covariate != "Cones") %>%
ggplot(aes(x = lag, y = `50%`)) +
  geom_hline(yintercept = 1/13, linetype = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`,
                    ymax = `97.5%`),
                width = 0.2) +
  facet_wrap(~covariate, scales = "free_x")+
  scale_x_continuous(breaks = c(1:13),
                     labels = c("0(Br)", "1(SW)", "2(Su)", "3(Fl)", "4(Br)", 
                                "5(SW)", "6(Su)", "7(Fl)", '8(Br)',
                                "9(SW)", "10(Su)", "11(Fl)", "12(Br)")) +
  labs(x = "Season into the past",
       y = "Importance weight \n (median and 95% BCI)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)))

ggsave(here('pictures',
            'final',
            'climate_weight_plots.jpg'),
       width = 5,
       height = 3,
       units = 'in',
       dpi = 300)

