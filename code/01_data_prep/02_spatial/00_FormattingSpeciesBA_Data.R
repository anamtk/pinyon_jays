################################################################################
## This script creates annual BA data for pinyon pine from 2010 to 2022, and modifies
  # pinyon cone predictions for analysis

## Code by Kyle C. Rodman, Ecological Restoration Institute. 
# 10/18/2024

################################################################################
### Bring in necessary packages
package.list <- c("here", "tidyverse", "sf", "terra", "USA.state.boundaries")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

################################################################################
### Get pinyon BA for Southwestern states

# ## Get pinyon pine basal area maps for contiguous US
# pinyonBA <- rast(here("Data", "SpeciesDominance", "PinyonBA_sqftPerAc.tif"))
# 
# ## Crop to extent of Four Corners States and save output in smaller area
# fourCorners <- USA.state.boundaries::state_boundaries_wgs84 %>%
#   filter(NAME == "Arizona"|NAME == "Colorado"|NAME=="New Mexico"|NAME == "Utah") %>%
#   st_transform(crs(pinyonBA))
# pinyonBA %>%
#   crop(ext(st_bbox(st_buffer(fourCorners, dist = 4000)))) %>%
#   writeRaster(here("Data", "SpeciesDominance", "PinyonBA_sqftPerAc_fourCorners.tif"))

################################################################################
### Read in the remaining/new data

## Get 4-km cone prediction data
pinyonCones <- rast(here("Data", "TimeSeries", "pinyon_jay_quantile_predictions.tif"))

## Get initial BA data
pinyonBA <- rast(here("Data", "SpeciesDominance", "PinyonBA_sqftPerAc_fourCorners.tif"))

## Get annual CC data
annualCC1 <- rast(here("Data", "SpeciesDominance", "fourCornersCC-0000000000-0000013824.tif"))
annualCC2 <- rast(here("Data", "SpeciesDominance", "fourCornersCC-0000000000-0000000000.tif"))

################################################################################
### Modify Pinyon Cone Predictions Using BA maps

## Threshold of < 0.1 sq m/ha
baProj <- terra::project(pinyonBA, pinyonCones[[1]], method = "average")
pinyonCones[baProj < 0.43560] <- 0 ## Any pixels with < 0.1 sq m/ha of pinyon BA are assumed to have zero cones
plot(pinyonCones[[1]])
writeRaster(pinyonCones, here("Data", "AnalysisReady", "ConePredictions_final.tif"), overwrite = T)

################################################################################
### Convert annual CC maps and initial pinyon BA maps to annual BA at 4km -
##  aligning with cone prediction grids

## Get baseline CC from 2000-2009, aligning with years and locations of pinyon BA
annualCC1_baseline <- app(annualCC1[[1:10]], mean)
annualCC2_baseline <- app(annualCC2[[1:10]], mean)
annualCC_baseline <- terra::merge(annualCC1_baseline, annualCC2_baseline)
annualCC_baseline <- project(annualCC_baseline, pinyonBA, method = "bilinear", gdal = T)

## With this, figure out the amount of BA represented by each % CC in each pixel
cellAdjustFactor <- pinyonBA/annualCC_baseline
cellAdjustFactor <- clamp(cellAdjustFactor, lower = 0, upper = 2.5, values = T)

## Now, loop through each year, figure out the change in CC from baseline period,
# and use that to adjust pinyon BA annually
for(i in 11:23){## Outer loop is the index of band in CC data. 11 = 2010 CC, 23 = 2022 CC
  ## First, merge and project CC data for that year
  tempCC1 <- annualCC1[[i]]
  tempCC2 <- annualCC2[[i]]
  annualCC <- terra::merge(tempCC1, tempCC2)
  annualCC <- project(annualCC, cellAdjustFactor, method = "bilinear", gdal = T)
  
  ## Next, get difference in CC from baseline period
  ccDiff <- annualCC-annualCC_baseline
  
  ## Use cell-specific BA to CC relationship to get change in pinyon BA
  baDiff <- ccDiff*cellAdjustFactor
  tempPinyon <- pinyonBA+baDiff
  tempPinyon <- clamp(tempPinyon, lower = 0, upper = 120, values = T)
  
  ## Project and align to same grid as Andreas' cone data, saving to output folder as 4k grid for corresponding year
  writeRaster(tempPinyon, here("Data", "SpeciesDominance", "PinyonBA_250m", paste("PinyonBA_",1999+i, ".tif", sep = "")), overwrite = T)
  finalPinyon <- terra::project(tempPinyon, pinyonCones, method = "average", gdal = T)
  writeRaster(finalPinyon, here("Data", "SpeciesDominance", "PinyonBA_4km", paste("PinyonBA_",1999+i, ".tif", sep = "")), overwrite = T)
}

## Merging all of them into one file and exporting - 250m res
allGrids <- rast(list.files(here("Data", "SpeciesDominance", "PinyonBA_250m"), full.names = T))
names(allGrids) <- paste("PinyonBA_sqftPerAc_", 2010:2022, sep = "") 
writeRaster(allGrids, here("Data", "SpeciesDominance", "PinyonBA_250m_sqftPerAc.tif"), overwrite = T)

## Merging all of them into one file and exporting - 4km res
allGrids <- rast(list.files(here("Data", "SpeciesDominance", "PinyonBA_4km"), full.names = T))
names(allGrids) <- paste("PinyonBA_sqftPerAc_", 2010:2022, sep = "") 
writeRaster(allGrids, here("Data", "SpeciesDominance", "PinyonBA_4km_sqftPerAc.tif"), overwrite = T)
writeRaster(allGrids/4.356, here("Data", "AnalysisReady", "PinyonBA_4km_sqmPerHa.tif"), overwrite = T)
