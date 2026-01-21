
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

# data <- readRDS(here('data',
#                       '01_ebird_data',
#                       'cleaned_data',
#                       'ebird_check_blob_yr_ids.RDS')) %>%
#   dplyr::select(numID, yrID, obsrvtn_c, checkID)

# yrep_sum <- readRDS(here('monsoon',
#                            'ebird',
#                            'nospuncert',
#                            'outputs',
#                            'ebird_abund_model_yrep.RDS'))

# yrep_sum2 <- as.data.frame(yrep_sum$quantiles) %>%
#   rownames_to_column(var = "parm")%>%
#   filter(parm != "deviance") %>%
#   separate(parm, 
#            into = c("yrID", "numID", "checkID"),
#            sep = ",") %>%
#   mutate(yrID = str_sub(yrID, 17, nchar(yrID)),
#          checkID =str_sub(checkID, 1, (nchar(checkID)-1)),
#          yrID =as.numeric(yrID),
#          numID = as.numeric(numID),
#          checkID = as.numeric(checkID))

r2 <- readRDS(here('monsoon',
                        'outputs',
                        'ebird_abund_model_yyrepr2.RDS')) %>%
  mutate(model = "full") 

r2_nocone <- readRDS(here('monsoon',
                        'outputs',
                        'ebird_abund_model_yyrepr2_nocone.RDS'))%>%
  mutate(model = "nocone")

r2_noppt <- readRDS(here('monsoon',
                        'outputs',
                        'ebird_abund_model_yyrepr2_noppt.RDS'))%>%
  mutate(model = "noppt") 


# COmbine -----------------------------------------------------------------

r2_all <- r2 %>%
  bind_rows(r2_nocone, r2_noppt)


# Plot --------------------------------------------------------------------

ggplot(r2_all) +
  geom_histogram(aes(x = R2, fill = model), alpha = 0.6)

ggplot(r2_all) +
  geom_boxplot(aes(x = model, 
                   y = R2,
                   fill = model), alpha = 0.6)

# Get proportions ---------------------------------------------------------

R2_full <- r2 %>%
  summarise(Mean = mean(R2)) %>%
  as_vector()

noconeR2 <- r2_nocone%>%
  summarise(Mean = mean(R2)) %>%
  as_vector()

nopptR2 <- r2_noppt%>%
  summarise(Mean = mean(R2)) %>%
  as_vector()

cone_prop <- (R2_full-noconeR2)/R2_full
ppt_prop <- (R2_full-nopptR2)/R2_full
(noconeR2)/R2_full
(nopptR2)/R2_full
(R2_full-noconeR2)/R2_full
(R2_full-nopptR2)/R2_full


# Stacked bar -------------------------------------------------------------

r2_sum <- r2_all %>%
  group_by(model) %>%
  summarise(mean_r2 = mean(R2)) %>%
  ungroup() 

ggplot(r2_sum) +
  geom_bar(aes(x = 1, 
               y = mean_r2, 
               fill = model),
           stat = 'identity',
           position = 'stack')


