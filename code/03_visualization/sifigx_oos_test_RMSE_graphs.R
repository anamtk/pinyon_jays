
# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 'data.table',
                  'corrplot',
                  'sf', 'coda', 'patchwork')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data objects -------------------------------------------------------

data <- readRDS(here('data',
                     '03_jags_input_data',
                     'oos',
                     'oos_ebird_data_list_nospuncert.RDS'))

test <- readRDS(here("data",
                     "03_jags_input_data",
                     "ebird_data_list_nospuncert.RDS"))

mean_y <- mean(c(data$ebird.count, test$ebird.count), na.rm = T)

rmse_samples <- readRDS(here('monsoon',
                             'outputs',
                             'ebird_abund_model_RMSE_samples.RDS'))

rmse_df <- bind_rows(as.data.frame(rmse_samples[[1]]),
                     as.data.frame(rmse_samples[[2]]),
                     as.data.frame(rmse_samples[[3]])) %>%
  dplyr::select(RMSE) %>%
  mutate(type = "test")

oos_RMSE_df <- readRDS(here('data',
                          '05_cross_validation',
                          'oos_RMSE.RDS')) %>%
  mutate(type = "oos") %>%
  rename(RMSE = oos_RMSE)

# Plot both RMSE ----------------------------------------------------------

both_RMSE <- oos_RMSE_df %>%
  bind_rows(rmse_df) %>%
  mutate(nRMSE = RMSE/mean_y)
  
theme_set(theme_bw())
a <- ggplot() +
  geom_boxplot(data = both_RMSE, aes(x = type, y = RMSE)) +
  labs(x = "Dataset") +
  scale_x_discrete(labels = c("out-of-sample", "in-sample (test)"))
a

ggsave(here('pictures',
            'final',
            'RMSE.jpg'),
       width = 3,
       height = 3,
       units = 'in')    

b <- ggplot() +
  geom_boxplot(data = both_RMSE, aes(x = type, y = nRMSE))

a+ b

both_RMSE %>%
  group_by(type) %>%
  summarise(mean = mean(RMSE),
            sd = sd(RMSE),
            total = n(),
            se = sd/sqrt(total))
