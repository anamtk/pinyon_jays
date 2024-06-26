# Load packages -----------------------------------------------------------

package.list <- c("tidyverse", 'here', #general packages
                  'jagsUI', #jags wrapper
                  'coda', #gelman.diag() function
                  'mcmcplots') #trace plot function

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_classic())


# Load data ---------------------------------------------------------------

data_list <- readRDS(here('data',
                          'simulated_data',
                          'data_list.RDS'))


# Model path --------------------------------------------------------------


model_file <- (here('code',
               'models',
               'ebird_bbs_joint_abund_JAGS.R'))


# Parameters to save ------------------------------------------------------


parameters <- c('a0',
                'a',
                'wA',
                'wB',
                'wC',
                'b0',
                'b',
                'c0',
                'c1',
                'c')

# Initials ----------------------------------------------------------------

inits <- readRDS(here('data',
                      'simulated_data',
                      'inits_list.RDS'))
# Run model ---------------------------------------------------------------

model <- jagsUI::jags(data = data_list,
                      parameters.to.save = parameters,
                      inits = inits,
                      model.file = model_file,
                      parallel = TRUE,
                      n.chains = 3,
                      n.burnin = 2000,
                      n.iter = 6000,
                      DIC = TRUE)


# Diagnose model ----------------------------------------------------------

mcmcplot(model$samples)


# Plots -------------------------------------------------------------------

sum <- summary(model$samples)

quants <- as.data.frame(sum$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "a0|a|wA|wB|wC")) %>%
  filter(parm != 'deviance') %>%
  mutate(ID = case_when(parm == "a0" ~ "Intercept",
                        parm == "a[1]" ~ "Cone Production",
                        parm == 'a[2]' ~ "Temperature",
                        parm == "a[3]" ~ "Precipitation"))

quants %>%
  filter(parm %in% c('a[1]', 'a[2]', 'a[3]')) %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype =2) +
  geom_pointrange(aes(x = `50%`,
                      xmin = `2.5%`,
                      xmax = `97.5%`,
                      y = ID)) +
  labs(x = "Covariate effect \n (median and 95% BCI)",
       y = "Covariate")

quants %>%
  filter(str_detect(parm, "wA")) %>%
  mutate(lag = str_sub(parm, 4, (nchar(parm)-1))) %>%
  ggplot() +
  geom_hline(yintercept = 1/5, linetype = 2) +
  geom_pointrange(aes(y = `50%`,
                      ymax = `97.5%`,
                      ymin = `2.5%`,
                      x = lag)) +
  labs(x = "Time into past",
       y = "Importance weight \n (median and 95% BCI)",
       title = "Cone production antecedent effects")


# Pull out N from model ---------------------------------------------------

parms2 <- c("N.tot")

mod2 <- update(model,
               parameters.to.save = parms2,
               n.iter = 375)

sum2 <- summary(mod2$samples)

pops <- as.data.frame(sum2$quantiles) %>%
  rownames_to_column(var = 'parm') %>%
  filter(parm != "deviance") %>%
  mutate(year = str_sub(parm, 7, (nchar(parm)-1)))

ggplot(pops, aes(x = year, y = `50%`, group = 1)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  geom_line() +
  labs(x = "Year", y = "Total relative population size \n (median and 95% BCI)") +
  scale_x_discrete(labels = c('2010', '2011', '2012', '2013', '2014'))
