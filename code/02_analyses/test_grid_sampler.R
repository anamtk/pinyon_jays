# Testing the covariate paraemter sampler
# Ana Miller-ter Kuile
# March 2, 2023

# this script tests the sampler Kiona generated that lets me
# sample from the posterior samples of each covariate parameter
# from statistical models for the full IPM

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "jagsUI")


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Generate data -----------------------------------------------------------

set.seed(1)
n.years <- 2
#n.ebird.grids is yearly (grids[t])
n.ebird.grids <- c(5,2)
#background grids (indexed by site and year [t,i])
n.background.grids <- matrix(data = rpois(n.years*max(n.ebird.grids), lambda = 4),
                             nrow = n.years,
                             ncol = max(n.ebird.grids))

grid_df <- as.data.frame(n.background.grids) %>%
  rownames_to_column(var = 'year') %>%
  pivot_longer(V1:V5,
               names_to = 'ebird_grid',
               values_to = "background_grid_count") %>%
  mutate(ebird_grid = str_sub(ebird_grid, 2, nchar(ebird_grid)),
         ebird_grid = as.numeric(ebird_grid)) %>%
  filter(!(year == 2 & ebird_grid > 2)) %>%
  mutate(background_grid = factor(background_grid_count,
                                  levels = c('1', '2', '3',
                                             '4', '5','6',
                                             '7', '8', '9', '10'))) %>%
  group_by(year, ebird_grid) %>%
  complete(background_grid) %>%
  fill(background_grid_count, .direction = "updown") %>%
  mutate(background_grid = as.numeric(background_grid),
         year = as.numeric(year)) %>%
  filter(background_grid <= background_grid_count) %>%
  dplyr::select(-background_grid_count) %>%
  mutate(
    weight = runif(n()),
    weight = weight / sum(weight)) %>%
  ungroup() %>%
  group_by(year) %>%
  mutate(ebird_grid = as.factor(ebird_grid)) %>%
  complete(ebird_grid) %>%
  ungroup() %>%
  mutate(background_grid = as.factor(background_grid)) %>%
  group_by(year, ebird_grid) %>%
  complete(background_grid) %>%
  mutate(background_grid = as.numeric(background_grid)) %>%
  replace_na(list(weight = 0)) %>%
  filter(!is.na(background_grid))

#now, generate IDs for the for loop where
# we will populate the array
gridID <- grid_df$ebird_grid
yearID <- grid_df$year
backgroundID <- grid_df$background_grid

#the "data" to randomly sample from, which is year, site, background grid [t,i,g]
#make a blank array
grid.array <- array(NA, dim = c(n.years, max(n.ebird.grids), max(n.background.grids)))

#fill taht array based on the values in those columns
for(i in 1:dim(grid_df)[1]){ #dim[1] = n.rows
  #fill grid array with the values for the background points
  #for that year and ebird grid
  grid.array[yearID[i], gridID[i],  backgroundID[i]] <- as.numeric(grid_df[i,3])
}


#pi is a vector of probabilities for a "background grid"
#for each year and ebird grid

#pi is indexed year, ebirdgrid, background grid pi[t,i,g]
pi <- array(NA, dim = c(n.years, max(n.ebird.grids), max(n.background.grids)))

#fill taht array based on the values in those columns
for(i in 1:dim(grid_df)[1]){ #dim[1] = n.rows
  #fill grid array with the values for the background points
  #for that year and ebird grid
  pi[yearID[i], gridID[i],  backgroundID[i]] <- as.numeric(grid_df[i,4])
}


data <- list(n.years = n.years,
             n.ebird.grids = n.ebird.grids,
             n.background.grids = n.background.grids,
             grid.array = grid.array,
             pi = pi)

# Parameters to save ------------------------------------------------------

parms <- c('grid.index',
           'ebird.grid')


# Run model ---------------------------------------------------------------

model <- here("code", 
              'jags_models',
              'grid_sampler.R')


mod <- jagsUI::jags(data = data,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = parms,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 4000,
                    DIC = TRUE)

sum <- summary(mod$samples)

View(sum$quantiles)
