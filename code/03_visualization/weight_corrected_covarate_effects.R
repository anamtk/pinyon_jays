#Model summary figures
#March 3, 2025
#Ana Miller-ter Kuile

#exploreing results of the Ebird model without lags

# Load packages -----------------------------------------------------------

package.list <- c('tidyverse',
                  'here', 'patchwork') 

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_classic())
# Load model summary ------------------------------------------------------

betas <- readRDS(here('monsoon',
                      'ebird',
                      'nospuncert',
                      'outputs',
                      'ebird_abund_model2_summary.RDS'))

data_list <- readRDS(here('data',
                          'jags_input_data',
                          'ebird_data_list_nospuncert.RDS'))

ebird_blobIDs <- read.csv(here('data',
                               'ebird_data',
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
                       'spatial_data',
                       'cleaned_data',
                       'cones_weighted_mean_blob.csv')) 

temp <- read.csv(here('data',
                      'spatial_data',
                      'cleaned_data',
                      'temp_weighted_mean_blob.csv'))

ppt <- read.csv(here('data',
                     'spatial_data',
                     'cleaned_data',
                     'ppt_weighted_mean_blob.csv'))

# Get weighted values for lagged terms ------------------------------------

#cone median weights
conewt <- as.data.frame(betas$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "wA")) %>%
  dplyr::select(`50%`) %>%
  as_vector()

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

coneVals <- apply(cones2, MARGIN = 1, FUN = function(x){sum(x*conewt, na.rm = T)})

#tmax
tmxwt <- as.data.frame(betas$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "wB")) %>%
  dplyr::select(`50%`) %>%
  as_vector()

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

tempVals <- apply(temp2, MARGIN = 1, FUN = function(x){sum(x*tmxwt, na.rm = T)})


#ppt
pptwt <- as.data.frame(betas$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "wC")) %>%
  dplyr::select(`50%`) %>%
  as_vector()

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

pptVals <- apply(ppt2, MARGIN = 1, FUN = function(x){sum(x*pptwt, na.rm = T)})

# Pull out beta DF --------------------------------------------------------

covariates <- c('a0', "Cones", "Tmax", "PPT", "Monsoon",
                "PinyonBA", "ConexTmax", "ConexPPT",
                'ConexMonsoon', "ConexBA")
beta_df <- as.data.frame(betas$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "a")) %>%
  filter(!parm %in% c("deviance")) %>%
  bind_cols(covariates) %>%
  rename(covariate = `...7`) %>%
  dplyr::select(covariate, `2.5%`, `50%`, `97.5%`) %>%
  mutate(med = case_when(covariate == "Cones" ~ `50%`*sd(coneVals),
                           covariate == "Tmax" ~ `50%`*sd(tempVals),
                           covariate == "	PPT" ~ `50%`*sd(pptVals),
                           covariate == "ConexTmax" ~ `50%`*sd(coneVals)*sd(tempVals),
                           covariate == "	ConexPPT" ~ `50%`*sd(coneVals)*sd(pptVals),
                           TRUE ~ `50%`))

covariate_betas <- beta_df %>%
  filter(covariate != 'a0') %>%
  ggplot(aes(x = `50%`,
                    y = reorder(covariate, `50%`))) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point() +
  labs(x = "Covariate effect\n(median and 95% BCI)",
       y = "Covariate")

# Weights -----------------------------------------------------------------

weights_df <- as.data.frame(betas$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "w")) %>%
  mutate(covariate = case_when(str_detect(parm, "wA") ~ "Cones",
                               str_detect(parm, "wB") ~ "Temperature",
                               str_detect(parm, 'wC') ~ "Precipitation")) %>%
  mutate(lag = str_sub(parm, 4, (nchar(parm)-1)),
         lag = as.numeric(lag))


start <- c(0.5,4.5)
end <- c(4.5, 7.5)
times <- c("predict", "lag")

box_df <- as.data.frame(start) %>%
  bind_cols(end = end, times = times) %>%
  mutate(times = factor(times, levels = c("predict", "lag")))

cone_weight_plot <- weights_df %>%
  filter(covariate == "Cones") %>%
ggplot(aes(x = lag, y = `50%`)) +
  geom_rect(data = box_df, aes(xmin = start,
                               xmax = end,
                               group = times,
                               ymin = 0, 
                               ymax = 1,
                               fill = times),
            inherit.aes = F, alpha = 0.6) +
  scale_fill_manual(values = c('#b3cde3',
                               '#ccebc5')) +
  geom_hline(yintercept = 1/7, linetype = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`,
                    ymax = `97.5%`),
                width = 0.2) +
  facet_wrap(~covariate, scales = "free_x") +
  scale_x_continuous(breaks = c(1:7),
                     labels = c("3", "2", "1", "0", "1", 
                                "2", "3")) +
  annotate(geom = "text", x = 1.5, y = 0.9, label = "Years before cones") +
  annotate(geom = "text", x = 5.5, y = 0.9, label = "Years after cones") +
  labs(x = "",
       y = "Importance weight \n (median and 95% BCI)") 

# ggsave(here('pictures',
#             'cone_weights.pdf'),
#        width = 4.5, 
#        height = 3.5, 
#        units = 'in')
#green - lag
#ccebc5
#blue - predictive
#b3cde3
#yellow - immediate
#ffffcc

#Season 1: 2-4: breeding
#Season 2: 5-6: feeding dependent young
#Season 3: 7: unknown, just chillin?
#Season 4: 8-1: winter foraging, potentially seeking food elsewhere

clim_weights_plot <- weights_df %>%
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Partial plots potentially -----------------------------------------------

# 
# monsoon <- read.csv(here('data',
#                          'spatial_data',
#                          'cleaned_data',
#                          'monsoon_data_df.csv')) %>%
#   mutate(monsoon_scaled = scale(PRISM_ppt_30yr_normal_800mM4_07_bil))
# 
# cones <- read.csv(here('data',
#                        'spatial_data',
#                        'cleaned_data',
#                        'cone_masting_df.csv')) %>%
#   dplyr::select(cell, X2010:X2024) %>%
#   pivot_longer(X2010:X2024,
#                names_to = "year",
#                values_to = "cones") %>%
#   mutate(scale_cones = scale(cones))
# 
# 
# beta_wide <- beta_df %>%
#   dplyr::select(covariate, `50%`) %>%
#   pivot_wider(names_from = covariate,
#               values_from = `50%`)
# 
# cones <- seq(from= -2.5, to = 2.5, by = 0.2)
# monsoon <- seq(from = -2.5, to = 2.5, by = 0.1)
# 
# df <- as.data.frame(expand.grid(cones = cones, 
#                                 monsoon = monsoon)) %>%
#   mutate(a0 = beta_wide$a0,
#          aCones = beta_wide$Cones,
#          aMons = beta_wide$Monsoon,
#          aConexMons = beta_wide$ConexMonsoon) %>%
#   rowwise() %>%
#   mutate(loglambda = a0 + aCones*cones +
#            aMons*monsoon + aConexMons*cones*monsoon,
#          lambda = exp(loglambda)) %>%
#   ungroup() 
# 
# df %>%
#   filter(monsoon %in% c(-2.5, 0, 2.5)) %>%
#   mutate(monsoon_lev = case_when(monsoon == -2.5 ~ "low",
#                                  monsoon == 0.0 ~ "average",
#                                  monsoon == 2.5 ~ "high")) %>%
#   mutate(monsoon_lev = factor(monsoon_lev, levels = c("low", "average",
#                               "high")) )%>%
#   ggplot(aes(x = cones, y = lambda, group = monsoon_lev, color = monsoon_lev)) +
#   geom_line(linewidth = 1)+ 
#   scale_color_manual(values = c("#7bccc4", '#2b8cbe',
#                                 '#084081')) +
#   labs(x = "Cone availability",
#        y = bquote('Jays'/~meter^2)) +
#   guides(color=guide_legend(title="Monsoonality")) +
#   xlim(-2, 2.3) +
#   scale_y_continuous(limits = c(0, 0.0000045),
#                      labels = scales::comma) 
#   
# 
# 
# 
# 
# 
