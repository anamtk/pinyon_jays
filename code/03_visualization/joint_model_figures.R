#Model summary figures
#August 6, 2024
#Ana Miller-ter Kuile

#exploreing results of the BBS-EBIRD joint likelihood abundance model

# Load packages -----------------------------------------------------------

package.list <- c('tidyverse',
                  'here', 'patchwork') 

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())
# Load model summary ------------------------------------------------------

sum <- readRDS(here('monsoon',
                    'outputs',
                    'ebird_bbs_joint_abund_model_summary.RDS'))


# Pull out results --------------------------------------------------------

med <- as.data.frame(sum$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "a|wA|wB")) %>%
  filter(!parm %in% c("deviance", "a0"))

med_effects <- med %>%
  filter(parm %in% c("a[1]", 'a[2]')) %>%
  mutate(parm = case_when(parm == "a[1]" ~ "Cones",
                          parm == "a[2]" ~ "VPD"))

weights <- med %>%
  filter(str_detect(parm, "w")) %>%
  mutate(lag = case_when(parm == "wA[1]" ~ "3 years ago",
                         parm == "wA[2]" ~ "2 years ago", 
                         parm == "wA[3]" ~ "1 year ago",
                         parm == "wA[4]" ~ "This year",
                         parm == "wA[5]" ~ "1 year in future",
                         parm == "wA[6]" ~ "2 years in future",
                         parm == "wA[7]" ~ "3 years in future",
                         parm == "wB[1]" ~ "This year",
                         parm == "wB[2]" ~ "1 year ago", 
                         parm == "wB[3]" ~ "2 years ago",
                         parm == "wB[4]" ~ "3 years ago",
                         parm == "wB[5]" ~ "4 years ago",
                         parm == "wB[6]" ~ "5 years ago",
                         parm == "wB[7]" ~ "6 years ago")) %>%
  mutate(covariate = case_when(str_detect(parm, "wA") ~ "Cones",
                               TRUE ~ "VPD"))

(a <- ggplot(med_effects) +
    geom_pointrange(aes(x = `50%`, y = parm, xmin = `2.5%`, xmax = `97.5%`)) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(x = "Covariate effect \n (Median and 95% BCI)",
         y = "")) 

(b <- weights %>%
    filter(covariate == "Cones") %>%
    mutate(lag = factor(lag, levels = c("3 years ago", "2 years ago",
                                        "1 year ago", "This year", 
                                        "1 year in future", "2 years in future", 
                                        "3 years in future"))) %>%
    ggplot() +
    geom_hline(yintercept = 1/7, linetype = 2) +
    geom_vline(xintercept = 'This year', linetype = 2, color = "grey") +
    geom_pointrange(aes(x = lag, y = `50%`, ymin = `2.5%`, ymax = `97.5%`)) +
    labs(x = "Time period", 
         y = "Lag effect \n (Median and 95% BCI)",
         title = "Cone weights") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))


(c <- weights %>%
    filter(covariate == "VPD") %>%
    mutate(lag = factor(lag, levels = c("6 years ago", "5 years ago",
                                        "4 years ago", "3 years ago",
                                        "2 years ago","1 year ago", 
                                        "This year"))) %>%
    ggplot() +
    geom_hline(yintercept = 1/7, linetype = 2) +
    geom_vline(xintercept = 'This year', linetype = 2, color = "grey") +
    geom_pointrange(aes(x = lag, y = `50%`, ymin = `2.5%`, ymax = `97.5%`)) +
    labs(x = "Time period", 
         y = "Lag effect \n (Median and 95% BCI)",
         title = "VPD weights") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))

a + (b/c) +
  plot_annotation(tag_levels = "A")

ggsave(here('pictures',
            'bbs_ebird_joint_results.jpg'),
            height = 5,
            width = 7, 
            units = "in")
