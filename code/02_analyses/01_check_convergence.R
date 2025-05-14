# Load packages -----------------------------------------------------------

package.list <- c("tidyverse", 'here' #general packages
                  )

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

#theme_set(theme_classic())


# Load Rhat ---------------------------------------------------------------

mod_rhat <- readRDS(here('monsoon',
                         'ebird',
                         'nospuncert',
                         'outputs',
                         'ebird_abund_model2_Rhat.RDS'))


# Graph RHat per parameter ------------------------------------------------

rhat_graph_fun <- function(list){
  
  #this creates a dtaaframe out of the Rhat values from the model
  df <- data.frame(id = names(list),
                   Rhat = unlist(lapply(list, paste, collapse = ","))) %>%
    #splits Rhat by , when that list element had more than one value
    mutate(Rhat = str_split(Rhat, ",")) %>%
    #unnests - so makes each Rhat a new row in the df
    unnest(c(Rhat)) %>%
    #make sure Rhat is a numeric
    mutate(Rhat = as.numeric(Rhat)) 
  
  #plot histogram and make sure all below 1.1
  plot <- ggplot(df, aes(x = Rhat)) +
    geom_histogram() +
    geom_vline(xintercept = 1.1, linetype = 2) +
    theme_bw() +
    scale_y_sqrt() +
    #this is facetted by the parameter ID so you
    # can find problematic parameters
    facet_wrap(~ id)
  
  return(plot)
}

rhat_graph_fun(list = mod_rhat)

ggsave(here('pictures',
            'final',
            'rhat.jpg'),
       width = 5,
       height = 4,
       units = 'in')       

