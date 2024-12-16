#Get the list of possible grid cell IDs per observation
#based on creating a buffer polygon around each 
#(5 miles for all BBS; 1/2 distance traveled for eBIRD)
#get the proportional coverage per each for each observation
#as a way to weight them

model{
  
  for(t in 1:n.years){ #each year
    for(i in 1:n.ebird.grids[t]){ #has diff # grids
      #random index for the grid array - which background point to collect
      grid.index[t,i] ~ dcat(pi[t,i,])
      #pulling ut the grid ID for that one
      ebird.grid[t,i] <- grid.array[t,i, grid.index[t,i]]
    }
  }
}