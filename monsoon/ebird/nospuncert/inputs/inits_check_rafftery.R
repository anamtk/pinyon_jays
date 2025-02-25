
# Load packages -----------------------------------------------------------

(start.time <- Sys.time())

package.list <- c("jagsUI", "coda",
                  'dplyr', 'stringr',
                  'magrittr', 'tidyr',
                  'tibble', 'purrr',
                  'mcmcplots','ggplot2') 

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load data ---------------------------------------------------------------

mod <- readRDS("/scratch/atm234/pinyon_jays/ebird/nospuncert/outputs/ebird_abund_model.RDS")

# Get initials from previous model ----------------------------------------

#get the MCMC chains
samples <- mod$samples

#function to make each chain a dataframe
df_fun <- function(chain){
  df <- as.data.frame(chain) %>%
    rownames_to_column(var = "iteration")
  return(df)
}

#use that function on all list elements
samp_dfs <- lapply(samples, df_fun)

#make into one dataframe
samp_df <- bind_rows(samp_dfs, .id = "chain")

#get values for all parameters from the last iteration of the
#chain with the lowest deviance
samp_df2 <- samp_df %>%
  group_by(chain) %>%
  #get mean deviance by chain
  mutate(mean_dev = mean(deviance, na.rm = T)) %>%
  ungroup() %>%
  #get only the chain with the minimum average deviance
  filter(mean_dev == min(mean_dev)) %>%
  #pull out the final iteration from that chain
  filter(iteration == max(iteration)) %>%
  dplyr::select(-chain, -iteration,
                -deviance, -mean_dev) 

samp_df2

#root nodes
#a0
a0 <- as.vector(samp_df2$a0)
#a
a <- samp_df2 %>%
  dplyr::select('a[1]':'a[9]') %>%
  pivot_longer('a[1]':'a[9]',
               names_to = 'a',
               values_to = "value") %>%
  dplyr::select(value) %>%
  as_vector()

#cs
c0 <- as.vector(samp_df2$c0)

c1 <- samp_df2 %>%
  dplyr::select("c1[2]") %>%
  as_vector()

c1 <- c(NA, c1)

c <- samp_df2 %>%
  dplyr::select('c[2]':'c[5]') %>%
  pivot_longer('c[2]':'c[5]',
               names_to = 'c',
               values_to = "value") %>%
  dplyr::select(value) %>%
  as_vector()

# Initials ----------------------------------------------------------------

inits_list <- readRDS('/scratch/atm234/pinyon_jays/ebird/nospuncert/inputs/ebird_init_list_nospuncert.RDS')

N <- inits_list[[1]]$N

inits_list2 <- list(list(N = N,
                         a0 = a0,
                         a = a,
                         c0 = c0,
                         c1 = c1,
                         c = c),
                    list(N = N,
                         a0 = a0 + 0.5,
                         a = a+ 0.25,
                         c0 = c0 +0.05 ,
                         c1[2] = c1[2] + 0.05,
                         c = c + 0.05),
                    list(N = N,
                         a0 = a0 - 0.5,
                         a = a- 0.25,
                         c0 = c0 -0.05 ,
                         c1[2] = c1[2] - 0.05,
                         c = c - 0.05))

inits_list2

# Raftery -----------------------------------------------------------------

raf <- raftery.diag(mod$samples)

names <- rownames(raf[[1]]$resmatrix)
ch1 <- raf[[1]]$resmatrix[,2]
ch2 <- raf[[2]]$resmatrix[,2]
ch3 <- raf[[3]]$resmatrix[,2]

raf_all <- as.data.frame(cbind(names, 
                               ch1, ch2, ch3)) %>%
  mutate(ch1 = as.numeric(ch1),
         ch2 = as.numeric(ch2),
         ch3 = as.numeric(ch3)) %>%
  pivot_longer(ch1:ch3,
               names_to = "chain",
               values_to = 'iterations') 

raf_all %>%
  summarise(iterations_90 = quantile(iterations, 
                                     probs = 0.9, 
                                     na.rm = T)/3,
            iterations_95 = quantile(iterations,
                                     probs = 0.95,
                                     na.rm = T)/3,
            max = max(iterations, 
                      na.rm = T)/3)

bu1 <- raf[[1]]$resmatrix[,1]
bu2 <- raf[[2]]$resmatrix[,1]
bu3 <- raf[[3]]$resmatrix[,1]

burn <- as.data.frame(cbind(names, bu1, bu2, bu3)) %>%
  mutate(bu1 = as.numeric(bu1),
         bu2 = as.numeric(bu2),
         bu3 = as.numeric(bu3)) %>%
  filter(!str_detect(names, "z")) %>%
  pivot_longer(bu1:bu3,
               names_to = "chain",
               values_to = 'iterations') 

burn %>%
  summarise(max(iterations, na.rm = T))



