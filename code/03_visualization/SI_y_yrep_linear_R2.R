
# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 'data.table',
                  'corrplot',
                  'sf', 'coda', 'patchwork')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# -------------------------------------------------------------------------

data <- readRDS(here('data',
                      '01_ebird_data',
                      'cleaned_data',
                      'ebird_check_blob_yr_ids.RDS')) %>%
  dplyr::select(numID, yrID, obsrvtn_c, checkID)

yrep_sum <- readRDS(here('monsoon',
                           'ebird',
                           'nospuncert',
                           'outputs',
                           'ebird_abund_model_yrep.RDS'))

yrep_sum2 <- as.data.frame(yrep_sum$quantiles) %>%
  rownames_to_column(var = "parm")%>%
  filter(parm != "deviance") %>%
  separate(parm, 
           into = c("yrID", "numID", "checkID"),
           sep = ",") %>%
  mutate(yrID = str_sub(yrID, 17, nchar(yrID)),
         checkID =str_sub(checkID, 1, (nchar(checkID)-1)),
         yrID =as.numeric(yrID),
         numID = as.numeric(numID),
         checkID = as.numeric(checkID))

yyrepr2 <- readRDS(here('monsoon',
                        'ebird',
                        'nospuncert',
                        'outputs',
                        'ebird_abund_model_yyrepr2.RDS'))

# y ~ yrep graph ----------------------------------------------------------

y_yrep_med <- yrep_sum2 %>%
  left_join(data, by = c("yrID", "numID", "checkID"))

linear <- ggplot(y_yrep_med, aes(x = obsrvtn_c, y = `50%`)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point() +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  labs(x = "Observed count", y = "Predicted count")

# R2 from samples ---------------------------------------------------------

mean_r2 <- yyrepr2 %>%
  summarise(mean = mean(V1, na.rm = T))

yyrepr2 %>%
  summarise(mean = mean(V1, na.rm = T),
            sd = sd(V1, na.rm = T),
            total = n(),
            se = sd/sqrt(total))

r2 <- ggplot(yyrepr2) +
  geom_histogram(aes(x = V1)) +
  geom_vline(xintercept = mean_r2$mean, linetype = 2) +
  labs(x = expression(paste(R^2)), y = "Count")

linear + r2 +
  plot_layout(widths = c(2,1)) +
  plot_annotation(tag_levels = "a",
                  tag_suffix = ")")

ggsave(here('pictures',
            'final',
            'y_yrep_linearr2.jpg'),
       width = 5,
       height = 2.5,
       units = 'in')       
