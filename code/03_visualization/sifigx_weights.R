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

# Weights -----------------------------------------------------------------

weights_df <- as.data.frame(betas$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "w")) %>%
  mutate(covariate = case_when(str_detect(parm, "wA") ~ "Cones",
                               str_detect(parm, "wB") ~ "Temperature",
                               str_detect(parm, 'wC') ~ "Precipitation")) %>%
  mutate(lag = str_sub(parm, 4, (nchar(parm)-1)),
         lag = as.numeric(lag)) 

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

