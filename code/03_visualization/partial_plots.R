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
                      'ebird_abund_model_run1nolags_summary.RDS'))

monsoon <- read.csv(here('data',
                         'spatial_data',
                         'cleaned_data',
                         'monsoon_data_df.csv')) %>%
  mutate(monsoon_scaled = scale(PRISM_ppt_30yr_normal_800mM4_07_bil))

cones <- read.csv(here('data',
                       'spatial_data',
                       'cleaned_data',
                       'cone_masting_df.csv')) %>%
  dplyr::select(cell, X2010:X2024) %>%
  pivot_longer(X2010:X2024,
               names_to = "year",
               values_to = "cones") %>%
  mutate(scale_cones = scale(cones))

# Pull out beta DF --------------------------------------------------------
covariates <- c('a0', "Cones", "Tmax", "PPT", "Monsoon",
                "PinyonBA", "ConexTmax", "ConexPPT",
                'ConexMonsoon', "ConexBA")
beta_df <- as.data.frame(betas$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "a")) %>%
  filter(!parm %in% c("deviance")) %>%
  bind_cols(covariates) %>%
  rename(covariate = `...7`)

beta_wide <- beta_df %>%
  dplyr::select(covariate, `50%`) %>%
  pivot_wider(names_from = covariate,
              values_from = `50%`)

cones <- seq(from= -2.5, to = 2.5, by = 0.2)
monsoon <- seq(from = -2.5, to = 2.5, by = 0.1)

df <- as.data.frame(expand.grid(cones = cones, 
                                monsoon = monsoon)) %>%
  mutate(a0 = beta_wide$a0,
         aCones = beta_wide$Cones,
         aMons = beta_wide$Monsoon,
         aConexMons = beta_wide$ConexMonsoon) %>%
  rowwise() %>%
  mutate(loglambda = a0 + aCones*cones +
           aMons*monsoon + aConexMons*cones*monsoon,
         lambda = exp(loglambda)) %>%
  ungroup() 

df %>%
  filter(monsoon %in% c(-2.5, 0, 2.5)) %>%
  mutate(monsoon_lev = case_when(monsoon == -2.5 ~ "low",
                                 monsoon == 0.0 ~ "average",
                                 monsoon == 2.5 ~ "high")) %>%
  mutate(monsoon_lev = factor(monsoon_lev, levels = c("low", "average",
                              "high")) )%>%
  ggplot(aes(x = cones, y = lambda, group = monsoon_lev, color = monsoon_lev)) +
  geom_line(linewidth = 1)+ 
  scale_color_manual(values = c("#7bccc4", '#2b8cbe',
                                '#084081')) +
  labs(x = "Cone availability",
       y = bquote('Jays'/~meter^2)) +
  guides(color=guide_legend(title="Monsoonality")) +
  xlim(-2, 2.3) +
  scale_y_continuous(limits = c(0, 0.0000045),
                     labels = scales::comma) 
  

ggsave(here('pictures',
            'jay_cone_monsoon.pdf'),
       height = 3,
       width = 4.5,
       units = "in")





