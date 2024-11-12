data{
  
  for(t in 1:n.years){ #could make this different for ebird than other loops if the data has more time
    for(i in 1:n.ebird.check[t]){ #number of checklists in year t
      
      #Spatial uncertainty in eBIRD locations for populating
      #the covariates
      #pi and grid.array are "data" 
      #pi is the proportion of the likely survey in a given grid
      #cell, grid.array is the gridIDs that correspond to each of those
      #same probabilities and from which a sample grid ID is drawn
      #random index for the grid array - which background point to collect
      ebird.grid.index[t,i] ~ dcat(ebird.pi[t,i,1:n.ebird.cells[t,i]])
      #pulling ut the grid ID for that one
      ebird.grid[t,i] <- ebird.grid.array[t,i,ebird.grid.index[t,i]]  
      
    } 
    
    n.grids[t] <- max(ebird.grid[t,])
    
  } 
}
model{
  for(t in 1:n.years){
    N[t] ~ dpois(n.grids[t])
  }
}