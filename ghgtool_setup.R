## ---------------------------
##
## Script name: ghgtool_setup.R
##
## Purpose of script: setup the directory architecture for the ghg tool, 
##                    with a specific focus on the outputs. Also, populate
##                    the data_in directory
##
## Run before: all other ghg scripts
##             ghg_tables_create.R 
##
## Specific numbered tasks:
## 1 - create directories
## 2 - populate data_in directory
## 3 - create project packages folder (proj_packages.R)
##
## Author: Dr. Paul M. Evans
##
## Date Created: 2023-01-13
##
## Copyright (c) Paul M. Evans, 2023
## Email: paueva@ceh.ac.uk
## Git: https://github.com/pevans13

## set working directory 
print(getwd())      # Working directory 
# setwd("C:/wd")    # Set working directory (PC)

## ---------------------------
options(scipen = 6, digits = 4) # for non-scientific notation
## ---------------------------

# remove objects
rm(list = ls())

# #### 0 - local or not ####
local <- T
if(local){
  localPath <- getwd()
}
rm(local)

#### 0 - which parts to run? ####
## ------------ Notes --------------  ##
## here you can select parts which you would like to run
## To do this, change the right side of the relevant arrow to a 'T'
## ------------ ----- --------------  ##
part1 <- T # create architecture
part2 <- T # create grid that will be used for the whole project
part3 <- T # imports data from the original place they were stored
part4 <- T # create an initial script that contains all the R packages required for the project

#### 0 - load libraries ####
library(raster)
library(dplyr)
library(sf)
library(tictoc)
library(data.table)

if(part1){
  #### 1 - create directories ####
  tic("creating directories")
  # rscripts (for r scripts)
  dir.create("r_scripts/", showWarnings = F, recursive = T)
  # data inputs
  dir.create("data_in/", showWarnings = F, recursive = T)
  # outputs
  dir.create("results/", showWarnings = F, recursive = T)
  ## specific categories
  catList <- c("land_use", "margins", "hedgerow", "perennial", "arable", "animals")
  for (i in 1:length(catList)){
    # datasets
    dir.create(paste0("results/", catList[i],"/"), showWarnings = F, recursive = T)
    # images/ figures
    dir.create(paste0("images/", catList[i],"/"), showWarnings = F, recursive = T)
  }
  
  # data inputs
  catList <- c("crops", "land_cover", "animals")
  for (i in 1:length(catList)){
    # datasets
    dir.create(paste0("data_in/", catList[i],"/"), showWarnings = F, recursive = T)
  }
  
  # tidy
  rm(catList, i)
  toc(log = T)
}

#### part 2 - create grid that will be used for the whole project ####
if(part2){
  #### 2a - grid ####
  tic("creating grid")
  ## ------------ Notes --------------  ##
  ## A grid can be created from the a replicate raser grid based on the initial 
  ## archetype data grid, which was 1 km.
  ## This will set the basis for all future calculations.
  ## ------------ ----- --------------  ##
  rasterPlan <- raster(ncol = 625, nrow = 1225
                       , xmn = 50000, xmx = 675000
                       , ymn = 0, ymx = 1225000)
  crs(rasterPlan) <- CRS('+init=EPSG:27700')
  print(rasterPlan)
  
  # create grid
  ## get extent
  e <- extent(rasterPlan)
  ## create bbox from it
  ebb <- st_bbox(c(xmin = xmin(e), xmax = xmax(e), ymax = ymax(e), ymin = ymin(e))
                 , crs = st_crs(27700))
  
  # create grid
  eGrid <- st_make_grid(ebb
                        , cellsize = c(1000, 1000)
                        , what = "polygons"
                        , square = TRUE) %>%
    st_as_sf(crs = st_crs(27700)) %>%
    rename(geometry = x)
  st_geometry(eGrid)
  str(eGrid)
  
  # get a reference for each cell, using centroid points in EPSG: 27700
  ## get centroid
  eGrid2Centroid <- eGrid %>%
    st_centroid()
  st_geometry(eGrid2Centroid)
  st_coordinates(eGrid2Centroid)
  head(eGrid2Centroid)
  ## get x and y from coords
  eGrid2Centroid <- eGrid2Centroid %>%
    mutate(x = st_coordinates(.)[, 1]
           , y = st_coordinates(.)[, 2]
    )
  head(eGrid2Centroid)
  # combine to get a rcFid ref
  eGrid2Centroid <- eGrid2Centroid %>%
    mutate(rcFid_1km = paste(round(x,0)
                             , round(y, 0)
                             , sep = "_"))
  head(eGrid2Centroid)
  
  # join back to original polygons
  eGrid3 <- eGrid %>%
    st_join(eGrid2Centroid)
  # remove archetype references
  eGrid3 <- eGrid3 %>%
    dplyr::select(rcFid_1km)
  class(eGrid3)
  head(eGrid3)
  str(eGrid3)
  
  # save 
  st_write(eGrid3, "data_in/agland_grid_ESW.gpkg", append = F)
  
  # save a version of this grid in all the data_in directories
  x <- list.dirs("data_in/")
  for(i in 2:length(x)){ file.copy("data_in/agland_grid_ESW.gpkg"
                                   , x[i]) }
  toc(log = T)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path("data_in", "readme.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2023-27-01
Last update:",  format(Sys.Date()), "

The files in this directory contain the data used to create the outputs of the GHG tool. There are three directories:
  1. 'land_cover' = data input relating to land covers, including hedgerows and field margins
  2. 'crops' = data relating to crops
  3. 'animals' = data relating to animals

The other file in this directory is 'agland_grid_ESW.gpkg', which contains that initial 1 km2 grid that will be used in all future analysis for this GHGtool
The only column is 'rcFid_1km', which is the centroid point of 1 km reference id (EPSG: 27700) [x and y coordinates are seperated by '_']. This column and the values in it will be used to enusre all data is in the same km
Data:
    Spatial resolution: 1 km2
    Spatial extent: GB
    Temporal coverage: N/A
    Native projection: 27700")
  
  sink(file = NULL)
}

#### part 3 - imports data from the original place they were stored ####
if(part3){
  #### 3a - land cover ####
  tic("loaded land cover")
  
  # copy the 1 km resolution land cover
  file.copy("N:\\Data\\UK\\Land_cover\\LCM2015/Targets21/lcm2015_gb_1km_dominant_target_class.img"
            , "data_in/land_cover/")
  # copy the 25 m resolution land cover 2015
  file.copy("N:\\Data\\UK\\Land_cover\\LCM2015\\Targets21/LCM2015_25/LCM2015_GB.tif"
            , "data_in/land_cover/")
  # copy the 25 m resolution land cover 2019 - this is for animal calculations
  file.copy("N:/Data/UK/Land_cover/LCM2019/data/gb2019lcm25m.tif"
            , "data_in/land_cover/")
  toc()
  
  #### 3b - Crops ####
  cat("copyng crops...\n")
  file.copy("N:\\Data\\UK\\Land_cover/LCMPlusCrops/LCM2015Pluscrops/LCM2015PlusCrops.tif"
            , "data_in/land_cover/")
  
  #### 3c - animals ####
  # agcensus - numeric
  agList <- list.files("N:\\Data\\UK\\agcensus/", pattern = "agcensus_5km"
                       , full.names = T)
  lapply(agList, file.copy, "data_in/animals/")
  
  cat("copyng animals...\n")
  # agcensus - names
  file.copy("N:\\Data\\UK\\agcensus/England_2010_5k.csv"
            , "data_in/animals/")
  
  # cts data
  cat("copyng CTS...\n")
  file.copy("N:/Data/UK/CTS/BEEF_2020.csv"
            , "data_in/animals/")
  file.copy("N:/Data/UK/CTS/DAIRY_2020.csv"
            , "data_in/animals/")
  }

if(part4){
  #### 4 - create project packages folder (proj_packages.R) ####
  
  # creare preamble
  text1 <- "#### 0 - load libraries for the GHGtool project ####\n## automatic install of packages if they are not installed already"
  text2 <- "list.of.packages <- c("
  # create list of packages required
  packageList <- c('raster', 'stars'
                   ,'sf'
                   ,'dplyr'
                   ,'data.table'
                   ,'tictoc'
                   ,'terra'
                   ,'tidyr'
                   , 'nngeo'
                   , 'pbapply'
                   , 'stringr'
                   , 'readxl'
                   , 'tidyverse'
                   , 'xml2'
                   , 'rvest'
                   , 'spatialEco' 
                   , 'parallel'
                   , 'doParallel'
                   , 'janitor')
  for(i in 1:length(packageList)){
    packageList[i] <- paste0("'", packageList[i], "'")
    if(i > 1){
      packageList[i] <- paste0(",", packageList[i])
    }
  }
  text3 <- ")" 
  text4 <- 'new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}
# loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(library(package.i, character.only = TRUE)
  )
}
# tidy
rm(list.of.packages, new.packages, package.i)'

writeLines(c(text1, text2, packageList, text3, text4), "r_scripts/proj_packages.R")
}
