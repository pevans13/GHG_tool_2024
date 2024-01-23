## ---------------------------
##
## Script name: total_ghg_script.R
##
## Purpose of script: to combine all of the different elements of greenhouse gas
##                    emissions calculation into one script##                    
##                    
## Run after: ghgtool_setup.R
##
## ------------ Notes --------------  ##
## The script encompasses many different process / emissions 
## It is separated into different sections based on these different processes
## and sources of emissions

## The user can select which emissions they are interested in i.e. want to calculate
## by selecting the appropriate section in the 'run which parts?' below. To run
## a specific part, the user will have to place a 'T' to the right of the arrow
## for that section. The ones that do not need to be run should have an 'F' to 
## the right.

## Each section should be self-contained, meaning that you should be able to run
## a single section without the others. However, this does assume that the user
## has run all of the sections previously, and created the outputs that feed 
## into the next sections.
## ------------ ----- --------------  ##
##
## Specific numbered tasks:
## 1 - 
## 2 - 
## 3 - 
## 4 - 
##
## list of final outputs:
##    data_in/ghgGridGB.gpkg -> 1 km2 grid of inland GB that contains 'rcFid_1km', which is a unique cell identifier
##    data_in/land_cover_table.gpkg -> 
##    data_in/animals/CTSCowsGrid.gpkg ->
##    data_in/land_cover/land_cover_table_19.gpkg ->
##    data_in/animals/cow_grid_2020_1km.gpkg ->
##    data_in/animals/nonBovine_grid_2015_1km.gpkg ->
##
## Author: Dr. Paul M. Evans
##
## Date Created: 2024-01-16
##
## Copyright (c) Paul M. Evans, 2024
## Email: paueva@ceh.ac.uk
## Git: https://github.com/pevans13

## ---------------------------
options(scipen = 6, digits = 4) # for non-scientific notation
## ---------------------------

# remove objects
rm(list = ls())

#### 0 - local or not ####
## ------------ Notes --------------  ##
## for this section, the user specifies where the model is being run. 
## This section should be changed by the user to meet their set-up needs
## The important thing is that 'home' is set to the root directory where the 
## outputs are required to go and the inputs can be found
## ------------ ----- --------------  ##

# determine whether the script is run locally, or on a HPC
local <- T

if(local){
  home <- getwd()
} else {
  # set this to a path that will be used
  home <- setwd(file.path(getwd(), "/ghgtool/"))
}
rm(local)

#### 0 - load libraries ####
## ------------ Notes --------------  ##
## This sections uses 'proj_packages.R' to install and read in all library
## packages neccessary for this project. Should anymore need to be added,
## they can be added by altering part 4 of 'ghgtool_setup.R'. 
## ------------ ----- --------------  ##
source("./r_scripts/proj_packages.R")

#### 0 - paths ####

#### 0 - run which parts? ####
# select which bits need running
## ------------ Notes --------------  ##
## This is the section where the user chooses which section they want to run by
## putting a 'T' to the right of the specific arrows
## ------------ ----- --------------  ##
createGrid <- F

# part 1 creates a grid that contains the composition of 1 km2 in terms of land 
# cover. The land cover information comes from CEH's LCM2015PlusCrops (LCM2015) map
part1LandCoverGrid <- F

# part 2 creates a grid that contains the number of animals present in a 1 km2
# cell. It uses Defra CTS data and AgCensus data to derive the animal numbers
part2AnimalGrid <- F

# part 3 creates attribute tables for each type of animal
part3AnimalAttributes <- F

# part 4 calculates the CO2e values from enteric fermentation for individual animals
part4Enteric <- F

# part 5 calculates the amount of manure per animal and the emissions from its management 
part5ManureManagement <- F

# part 6 calculates the amount of excretions produced per animal
part6excretions <- T

# part 7 calculates the amount of N used as fertilisers, based on RB209, and then the
# emissions created by that N use. 
## Specific numbered tasks:
## 1. Determine pH, and thus how much lime may be needed
## 2. determine rainfall, crop type, and Soil Nitrogen Supply (SNS)
## 3. determine N use per km2 - derived from individual crops
## 4. calculate emissions based purely on assumed fertiliser
part7fertiliserUse <- F

# part 8 determines the emission and sequestration values based on non-agricultural land use
part8landuse <- F

# part 9 determines the emissions from crop residues (i.e. the non-harvested part of the plant)
# It does this by assume that residues are managed in a certain way
part9residue <- F

# part 10 determines on-farm energy in terms of CO2e emissions  
# this includes on-farm fuel use, and potato storage
part10onfarmenergy <- F

part11combineemissions <- F

#### 0 - select animal scenario ####
## ------------ Notes --------------  ##
## This sections allows you to choose what the set-up of animals will be. This
## relates to the quality of food, for both the inside portion (6 months) and
## the outside portion (6 months). The range should be between 60 and 85% 
## digestible energy (DE)
## ------------ ----- --------------  ##
DEinside <- 85
DEoutside <- 65

#### 0 - load the GB grid ####
## ------------ Notes --------------  ##
## The GB grid should have been created in the 'ghgtool_setup.R' script, and
## saved as 'agland_grid_ESW.gpkg'. That grid needs to be loaded. The script 
## below modifies it as to only keep the 1 km2 of inland GB.
## ------------ ----- --------------  ##

if(createGrid){
  ghgGridPath <- "data_in/agland_grid_ESW.gpkg"
  ghgGrid <- st_read(ghgGridPath)
  rm(ghgGridPath)
  
  # set location of files that relate to land cover data
  currPath <- file.path("./data_in", "land_cover")
  
  # load 25 m LCM raster
  lcmCrops25 <- rast(file.path(currPath, "LCM2015PlusCrops.tif"))
  
  ##### create a GB-only 1 km2 grid #####
  # convert the initial ghgGrid to points
  ghgGridPoints <- ghgGrid %>% 
    st_centroid()
  head(ghgGridPoints)
  ## extract raster values of 25 m LCM to points, to determine just GB pixels
  ghgGridPointsLCM <- terra::extract(lcmCrops25, ghgGridPoints) %>%
    ### where first layer is NA, remove
    filter(!is.na(Band_1))
  head(ghgGridPointsLCM)
  ### extract those ID points from the original ghgGrid
  ghgGridGB <- ghgGrid %>%
    slice(c(unique(ghgGridPointsLCM$ID)
            # other points need including
            , 361464, 123034, 123037, 123038, 123041, 122419, 122420, 115536, 115532
            , 454677, 454052, 455301, 298407, 298419, 296532, 295907, 294030, 294656
            , 349009))
  ### set the crs, which should be 27700
  st_crs(ghgGridGB) <- 27700
  ghgGridGB
  ### save the new grid
  st_write(ghgGridGB, "data_in/ghgGridGB.gpkg", append = F)
  # ghgGridGB <- st_read("data_in/ghgGridGB.gpkg")
  
  # tidy
  rm(ghgGridPoints, ghgGridPointsLCM, ghgGrid)
  
} else {
  
  # load in, if already created
  ghgGridGB <- st_read("data_in/ghgGridGB.gpkg")
  
}
rm(createGrid)

#### 1 - part1LandCoverGrid ####
if(part1LandCoverGrid){
  
  ## ------------ Notes --------------  ##
  ## The data used for crops in this tools comes from the LCM2015PlusCrops, which
  ## matches the LCM 2015 data. 
  ## ------------ ----- --------------  ##
  
  ##### 1a - set grid and import crops #####
  # set location of files that relate to land cover data
  currPath <- file.path("./data_in", "land_cover")
  
  # load 25 m LCM raster
  lcmCrops25 <- rast(file.path(currPath, "LCM2015PlusCrops.tif"))
  
  ##### 1b - extract land cover information #####
  tic("extract land cover information")
  ## ------------ Notes --------------  ##
  ## using the LCM 2015 with crops, the area, in hectares, of each land cover  
  ## will be obtained. 
  ## ------------ ----- --------------  ##
  
  # get nrows in GB polygon grid
  GBrows <- nrow(ghgGridGB)
  
  # extract land covers from underlying 25 m2 land cover raster
  # put it inside a loop, to be able to track progress
  ## storage list
  landCoversFor1km <- list()
  
  ### loop 
  #### create a sequence, in order to track progress of the loop
  lcmSeq <- round(c(seq(1, GBrows, GBrows/10), GBrows), 0)
  lcmSeq
  
  ## ------------ Notes --------------  ##
  ## here, the land cover types under each km2 are extracted. It should take about 25 mins
  ## ------------ ----- --------------  ##
  tic("total extraction time")
  landCovers <- pblapply(1:(length(lcmSeq) - 1), function(x){
    # landCovers <- pblapply(1:2, function(x){
    
    cat("start =", xStart <- lcmSeq[x], "| end = "
        , xEnd <- ifelse(x == (length(lcmSeq)-1)
                         , (lcmSeq[x+1]+1), lcmSeq[x+1])-1
        , "| at", format(Sys.time(), "%b %d %X"),  "\n")
    
    tic("extracted one tenth")
    landCoversFor1km[[x]] <- as.data.frame(table(terra::extract(lcmCrops25, ghgGridGB[xStart:xEnd, ])))
    # get list of grid row numbers
    gridrowNums <- ghgGridGB[xStart:xEnd, ]$rcFid_1km
    landCoversFor1km[[x]]$rcFid_1km <- gridrowNums
    
    toc(log = T)
    return(landCoversFor1km[[x]])
  })
  toc(log = T)
  
  #### bind them all together
  lc1km <- bind_rows(landCovers)
  head(lc1km)
  #### get unique values of ID - should match nrows
  stopifnot(identical(length(unique(lc1km$rcFid_1km)), GBrows))
  #### make wide, to be able to determine composition
  lc1kmWide <- lc1km %>%
    pivot_wider(names_from = Band_1, values_from = Freq) %>%
    # remove ID
    dplyr::select(-ID)
  names(lc1kmWide)[2:ncol(lc1kmWide)] <- paste0("lc_", names(lc1kmWide)[2:ncol(lc1kmWide)])
  head(lc1kmWide)
  #### convert names to land covers
  lc1kmWide <- lc1kmWide %>%
    # rename the columns for crops
    rename(winterwheat = lc_302
           , winterbarley = lc_303
           , springwheat = lc_304
           , oats = lc_305
           , maize = lc_306
           , rapeseed = lc_307
           , springbarley = lc_308
           , potato = lc_309
           , fieldbeans = lc_310
           , sugarbeet = lc_311
           # improved grassland
           , improved_grass = lc_401
           # others
           , broadleaf = lc_1
           , conifer = lc_2
           , neutral_grass = lc_5
           , calc_grass = lc_6
           , acid_grass = lc_7
           , fen_marsh = lc_8
           , heather = lc_9
           , heather_grass = lc_10
           , bog = lc_11
           , inland_rock = lc_12
           , saltwater = lc_13
           , freshwater = lc_14
           , sup_lit_rock = lc_15
           , sup_lit_sed = lc_16
           , lit_rock = lc_17
           , lit_sed = lc_18
           , saltmarsh = lc_19
           , urban = lc_20
           , suburban = lc_21) %>%
    # convert to hectares
    mutate(across(2:ncol(.), ~ (. * 625) / 10000))
  # change colnames to reflect that they are in hectares
  names(lc1kmWide)[2:ncol(lc1kmWide)] <- paste0(names(lc1kmWide)[2:ncol(lc1kmWide)], "_ha")
  head(lc1kmWide)
  #### make spatial
  stopifnot(length(unique(lc1kmWide$rcFid_1km)) == length(unique(ghgGridGB$rcFid_1km)))
  stopifnot(identical(sort(unique(lc1kmWide$rcFid_1km)), sort(unique(ghgGridGB$rcFid_1km))))
  lc1kmSpatial <- lc1kmWide %>% 
    merge(ghgGridGB
          , .
          , by = "rcFid_1km"
          , all = T)
  head(lc1kmSpatial)
  class(lc1kmSpatial)
  #### save
  st_write(lc1kmSpatial, file.path(currPath, "land_cover_table.gpkg"), append = F)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(currPath, "readme_land_cover_table.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-08-10
  Last update:",  format(Sys.Date()), "
  
  Description of 'land_cover_table.gpkg':
  This is a gridded dataset that contains the area, in ha, of different land classes in a 1 km2 cell. 
  The areas of land cover were derived from 25 m2 land cover map 2015 (Rowland et al., 2017)
  
  Citation for the original data: 
  Rowland, C. S., Morton, R. D., Carrasco, L., McShane, G., O’Neil, A. W., & Wood, C. M. (2017). Land Cover Map 2015 (25m raster, GB) [Data set]. NERC Environmental Information Data Centre. https://doi.org/10.5285/bb15e200-9349-403c-bda9-b430093807c7
  
  Columns of 'land_cover_table.gpkg':
    'rcFid_1km' = centroid point of 1 km reference id (EPSG: 27700) [x and y coordinates are seperated by '_']
    '[land_cover]_ha' = area of the respective land covers within the km, in ha
    
  Data: Land cover area
  units: ha
  Spatial resolution: 1 km2
  Spatial extent: GB
  Temporal coverage: 2015
  Native projection: 27700")
  sink(file = NULL)
  
} # end of part1LandCoverGrid

# tidy
rm(part1LandCoverGrid)

#### 2 - part2AnimalGrid ####
if(part2AnimalGrid){
  
  ## ------------ Notes --------------  ##
  ## The data used for animal numbers comes from a combination of CTS (defra data)
  ## and AgCensus data. The Livestock numbers comes from CTS, while the other animals
  ## are inferred based on their temporal difference between 2010 and 2020.
  ## ------------ ----- --------------  ##
  
  # set location of files that relate to animal data
  currPath <- file.path("./data_in", "animals")
  
  ##### 2a - load in the bovine data #####
  # load CTS (Defra 2020) data
  beefCows <- fread(file.path(currPath, "BEEF_2020.csv"))
  dairyCows <- fread(file.path(currPath, "DAIRY_2020.csv"))
  # see unique categories
  print(unique(beefCows$BREED_DESC))
  print(unique(beefCows$ROLE_DESC))
  print(unique(dairyCows$BREED_DESC))
  print(unique(dairyCows$ROLE_DESC))
  
  # sort data into one spatial table
  beefCows2 <- beefCows %>%
    tidyr::unite("Category", which(names(beefCows) == "BREED_DESC") : which(names(beefCows) == "ROLE_DESC")
                 , remove = T)
  dairyCows2 <- dairyCows %>%
    tidyr::unite("Category", which(names(dairyCows) == "BREED_DESC") : which(names(dairyCows) == "ROLE_DESC")
                 , remove = T)
  ## merge
  ctsCows <- merge(beefCows2, dairyCows2
                   , by = c("EASTING", "NORTHING", "Category", "TOTAL")
                   , all = T)
  ## make wide i.e. one row per cell
  ctsWide <- ctsCows %>% 
    pivot_wider(names_from = Category, values_from = TOTAL)
  ## fill in NAs as 0
  ctsWide[is.na(ctsWide)] <- 0
  ## make spatial
  ctsWide <- st_as_sf(ctsWide, coords = c("EASTING", "NORTHING")
                      , crs = 27700)
  plot(ctsWide[1])
  ## get the bounding box of a 5000 m (half of 10,000 m (1 km)) buffered circle around the point
  box <- st_buffer(ctsWide, 15000) %>%
    st_bbox() %>%
    st_as_sfc()
  ### make a 10x10 grid of the bounding box
  grid <- st_make_grid(box, cellsize = c(10000, 10000)) %>%
    st_as_sf()
  ### make spatial in polygon squares
  cowsGrid <- st_join(grid, ctsWide)
  ### remove rows that a all NA
  cowsGrid <- cowsGrid[rowSums(is.na(cowsGrid)) == 0, ]
  ### make fids
  cowsGrid$fid10km <- paste0("fid10km_", 1:nrow(cowsGrid))
  head(cowsGrid)
  plot(cowsGrid[1])
  #### write shapefile
  st_write(cowsGrid, file.path(currPath, "CTSCowsGrid.gpkg"), append = F)
  cat("cow grid at 10 km saved:", file.path(currPath, "CTSCowsGrid.gpkg"), "\n")
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(currPath, "readme_CTSCowsGrid.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-08-10
  Last update:",  format(Sys.Date()), "
  
  Description of 'CTSCowsGrid.gpkg':
  This is a gridded dataset that contains the quantity of different farmed bovine animals. 
  The bovines are split up into different 'purposes' and life stages e.g. 1-year-old dairy cow
  These data were obtained from Defra's Cattle Tracing System for the year 2020.
  
  Citation for the original data: 
  
  
  Columns of 'CTSCowsGrid.gpkg':
    'fid10km' =10 km reference id based on its row position within the 10 km dataset 
    '[animal_type]' = numbers of animals of that type present within a 10 km2 cell for 2020
    
  Data: Bovine numbers
  units: quantity
  Spatial resolution: 10 km2
  Spatial extent: GB
  Temporal coverage: 2020
  Native projection: 27700")
  sink(file = NULL)
  
  # tidy
  rm(beefCows, beefCows2, dairyCows, dairyCows2)
  
  ##### 2b - scale down bovine data from 10 km to 1km #####
  ## ------------ Notes --------------  ##
  ## The scaling down will be done based on the composition of grasslands 
  ## within each 1 km. The grassland data will come from the lcm2019 map,
  ## as it is the closest available data to 2020 cow data
  ## ------------ ----- --------------  ##  
  
  ###### 2b1 - create land cover map dataframe for 2019 ######
  
  # load if it exists, otherwise run through
  if(file.exists(file.path("./data_in", "land_cover", "land_cover_table_19.gpkg"))){
    grass19Grid <- st_read(file.path("./data_in", "land_cover", "land_cover_table_19.gpkg"))
  } else {
    # import grasslands from 2019 - closest year to cow data we have
    grass2019 <- rast(file.path("./data_in", "land_cover", "gb2019lcm25m.tif"))
    # create empty df to store land cover data, just for grass (i.e. 4-7 and 10)
    grass2019Df <- data.frame(grasscode = c(4:7, 10))
    
    # extract land covers from underlying 25 m2 land cover raster
    # put it inside a loop, to be able to track progress
    ## storage list
    landCoversFor1km <- list()
    
    # get nrows in GB polygon grid
    GBrows <- nrow(ghgGridGB)
    
    ### loop 
    #### create a sequence, in order to track progress of the loop
    lcmSeq <- round(c(seq(1, GBrows, GBrows/10), GBrows), 0)
    lcmSeq
    
    tic("total extraction time")
    landCovers2b1 <- pblapply(1:(length(lcmSeq) - 1), function(x){
      # landCovers <- pblapply(1:2, function(x){
      
      cat("start =", xStart <- lcmSeq[x], "| end = "
          , xEnd <- ifelse(x == (length(lcmSeq)-1)
                           , (lcmSeq[x+1]+1), lcmSeq[x+1])-1
          , "| at", format(Sys.time(), "%b %d %X"),  "\n")
      
      tic("extracted one tenth")
      landCoversFor1km[[x]] <- as.data.frame(table(terra::extract(grass2019[[1]], ghgGridGB[xStart:xEnd, ])))
      # get list of grid row numbers
      gridrowNums <- ghgGridGB[xStart:xEnd, ]$rcFid_1km
      landCoversFor1km[[x]]$rcFid_1km <- gridrowNums
      
      toc(log = T)
      return(landCoversFor1km[[x]])
    })
    toc(log = T)
    
    #### bind them all together
    lc1km <- bind_rows(landCovers2b1)
    head(lc1km)
    
    #### get unique values of ID - should match nrows
    stopifnot(identical(length(unique(lc1km$rcFid_1km)), GBrows))
    #### make wide, to be able to determine composition
    lc1kmWide <- lc1km %>%
      pivot_wider(names_from = gb2019lcm25m_1, values_from = Freq) %>%
      # remove ID
      dplyr::select(-ID)
    names(lc1kmWide)[2:ncol(lc1kmWide)] <- paste0("lc_", names(lc1kmWide)[2:ncol(lc1kmWide)])
    head(lc1kmWide)
    #### convert names to land covers
    lc1kmWide <- lc1kmWide %>%
      # rename the columns for crops
      rename(improved_grass = lc_4 # improved grassland
             # others
             , broadleaf = lc_1
             , conifer = lc_2
             , arable = lc_3
             , neutral_grass = lc_5
             , calc_grass = lc_6
             , acid_grass = lc_7
             , fen_marsh = lc_8
             , heather = lc_9
             , heather_grass = lc_10
             , bog = lc_11
             , inland_rock = lc_12
             , saltwater = lc_13
             , freshwater = lc_14
             , sup_lit_rock = lc_15
             , sup_lit_sed = lc_16
             , lit_rock = lc_17
             , lit_sed = lc_18
             , saltmarsh = lc_19
             , urban = lc_20
             , suburban = lc_21) %>%
      # convert to hectares
      mutate(across(2:ncol(.), ~ (. * 625) / 10000))
    # change colnames to reflect that they are in hectares
    names(lc1kmWide)[2:ncol(lc1kmWide)] <- paste0(names(lc1kmWide)[2:ncol(lc1kmWide)], "_ha")
    head(lc1kmWide)
    #### make spatial
    stopifnot(length(unique(lc1kmWide$rcFid_1km)) == length(unique(ghgGridGB$rcFid_1km)))
    stopifnot(identical(sort(unique(lc1kmWide$rcFid_1km)), sort(unique(ghgGridGB$rcFid_1km))))
    lc1kmSpatial <- lc1kmWide %>% 
      merge(ghgGridGB
            , .
            , by = "rcFid_1km"
            , all = T)
    head(lc1kmSpatial)
    class(lc1kmSpatial)
    #### save
    st_write(lc1kmSpatial, file.path("./data_in", "land_cover", "land_cover_table_19.gpkg"), append = F)
    
    # write readme
    sink(file = NULL)
    sink(file = file.path("./data_in", "land_cover", "readme_land_cover_table_19.md"))
    cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-08-10
  Last update:",  format(Sys.Date()), "
  
  Description of 'land_cover_table_19.gpkg':
  This is a gridded dataset that contains the area, in ha, of different land classes in a 1 km2 cell. 
  The areas of land cover were derived from 25 m2 land cover map 2019 (Morton et al., 2020)
  
  Citation for the original data: 
  Morton, R. D., Marston, C. G., O’Neil, A. W., & Rowland, C. S. (2020). Land Cover Map 2019 (25m rasterised land parcels, GB) [Data set]. NERC Environmental Information Data Centre. https://doi.org/10.5285/F15289DA-6424-4A5E-BD92-48C4D9C830CC
  
  Columns of 'land_cover_table.gpkg':
    'rcFid_1km' = centroid point of 1 km reference id (EPSG: 27700) [x and y coordinates are seperated by '_']
    '[land_cover]_ha' = area of the respective land covers within the km, in ha
    
  Data: Land cover area
  units: ha
  Spatial resolution: 1 km2
  Spatial extent: GB
  Temporal coverage: 2019
  Native projection: 27700")
    sink(file = NULL)
    
    # load
    grass19Grid <- st_read(file.path("./data_in", "land_cover", "land_cover_table_19.gpkg"))
  } # end of creating the 2019 land cover (for grasses, mostly)
  
  ###### 2b2 - extract grass per 10 km2 cells ######
  # create if it does not already exist
  if(file.exists(file.path(currPath, "cow_grid_2020_1km.gpkg"))){
    cows1kmGB <- st_read(file.path(currPath, "cow_grid_2020_1km.gpkg"))
  } else {
    # create centroids first
    grass19Cent <- grass19Grid %>% st_centroid() %>%
      # just keep grasses
      dplyr::select(rcFid_1km, improved_grass_ha
                    , neutral_grass_ha, calc_grass_ha, acid_grass_ha, heather_grass_ha) %>%
      # get total sum of grass
      mutate(grass_total_ha = rowSums(.[,2:(ncol(.)-1),drop=TRUE], na.rm = TRUE)) %>%
      dplyr::select(rcFid_1km, grass_total_ha)
    head(grass19Cent)
    
    tic("total extraction time for grass")
    # Set up parallel processing (adjust the number of cores as needed)
    registerDoParallel(cores = 4)
    
    mt <- list() # save list
    # Use foreach with %dopar%
    cAndG <- foreach(i = 1:nrow(cowsGrid)
                     , .combine = "rbind"
                     , .packages = c("tictoc", "sf", "dplyr", "tidyr")) %dopar% {
                       tic("one row")
                       # extract row
                       cowRow <- cowsGrid[i, ]
                       
                       # determine which 1 km cells each 10 km cell contains - use st_intersection
                       grassRow <- st_intersection(grass19Cent, cowRow)
                       cat("i =", i, "| number grassRow: ", nrow(grassRow), "\n")
                       if(nrow(grassRow) > 0){
                         # get grass in 1 km as a percentage of grass in 10 km2
                         grassRow$grassPC10km <- grassRow$grass_total_ha / sum(grassRow$grass_total_ha)
                         
                         # duplicate rows of the cow grid to equal the 1 km underlying it
                         cowRow2 <- rbind(cowRow %>%
                                            st_drop_geometry(), cowRow[rep(1, (nrow(grassRow)-1)), ] %>%
                                            st_drop_geometry()) %>% dplyr::select(-fid10km)
                         
                         # multiply, to get the quantity of each bovine category in each 1 km2
                         cowsAndGrass <- grassRow$grassPC10km * cowRow2
                         cowsAndGrass <- bind_cols(grassRow %>% st_drop_geometry() %>% dplyr::select(rcFid_1km), cowsAndGrass)
                         head(cowsAndGrass)
                         mt[[i]] <- cowsAndGrass
                         toc()
                         
                         # Return the result for each iteration
                         if (nrow(grassRow) > 0) {
                           return(cowsAndGrass)
                         } else {
                           return(NULL)
                         }
                       }
                     } 
    
    # Stop parallel processing
    stopImplicitCluster()
    toc(log = T)
    head(cAndG)
    
    # merge with GB grid
    cows1kmGB <- merge(ghgGridGB, cows1km, by = "rcFid_1km", all = T)
    cows1kmGB[is.na(cows1kmGB)] <- 0
    head(cows1kmGB[, 1:20])
    class(cows1kmGB)
    #### save
    st_write(cows1kmGB, file.path(currPath, "cow_grid_2020_1km.gpkg"), append = F)
    
    # write readme
    sink(file = NULL)
    sink(file = file.path(currPath, "readme_cow_grid_2020_1km.md"))
    cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-08-10
  Last update:",  format(Sys.Date()), "
  
  Description of 'cow_grid_2020_1km.gpkg':
  This is a gridded dataset that contains the quantity of different farmed bovine animals, based on the underlying grass area. 
  The bovines are split up into different 'purposes' and life stages e.g. 1-year-old dairy cow
  These data were obtained from Defra's Cattle Tracing System for the year 2020, which was at the 10 km2 resolution
  Using the amount of grass in each 1 km2 underlying the each 10 km2, the bovine numbers in a 10 km2 cell where proportied to each 1 km2 based on proportion of grass
  
  Citation for the original data: 
  
  
  Columns of 'cow_grid_2020_1km.gpkg':
    'rcFid_1km' = centroid point of 1 km reference id (EPSG: 27700) [x and y coordinates are seperated by '_']
    '[animal_type]' = derived (using proportion of grass in 10 km2) numbers of animals of that type present within a 1 km2 cell for 2020
    
  Data: Bovine numbers
  units: quantity
  Spatial resolution: 1 km2
  Spatial extent: GB
  Temporal coverage: 2020
  Native projection: 27700")
    sink(file = NULL)
    
  } # end of cow grid creation at 1 km2
  
  # tidy
  rm(cowsAndGrass, cAndG, cowsGrid, ctsCows, ctsWide)
  
  ###### 2b3 - determine non-bovine animals per 5 km ######
  
  ## ------------ Notes --------------  ##
  ## We could only get bovine data from Defra. Data on other animals came from
  ## AgCensus, of which data is available at 5 km2 resolution. These data include
  ## numbers of sheep, pigs, poultry, and horses.
  ## Assuming these animals are always on grasslands, a similar method will be
  ## used for these animals as it was for the bovines: they were calculated per
  ## grassland proportion. As the AgCensus data were from 2010, the numbers 
  ## needed to be updated to 2020 i.e. to match the cattle data
  ## This was done using the total number of pigs, sheep, and poultry between
  ## 2010 and 2020, and adjusting the proportion appropriately. In practice, this
  ## meant that each animal had more space / less space per individual depending on whether
  ## that animal type increased or decreased.
  
  ## Pigs decreased by 26% between 2010 and 2020, while sheep decreased 21%
  ## and poultry increased by 7%
  ## ------------ ----- --------------  ##
  
  # create if it does not already exist
  if(file.exists(file.path(currPath, "nonBovine_grid_2015_1km.gpkg"))){
    nonBov1kmGB <- st_read(file.path(currPath, "nonBovine_grid_2015_1km.gpkg"))
  } else {
    # load agcensus data
    AgCensus5kmOG <- st_read(file.path(currPath, "agcensus_5km.shp"))
    # adjust AgCensus so slightly that it does not count the saem underlying 1 km twice
    ## Define the amount to move the polygon (in metres)
    moveDistance <- c(10, 10)
    ## Move the polygon by adding the move_distance to each vertex
    AgGeom <- st_geometry(AgCensus5kmOG) + moveDistance
    st_geometry(AgCensus5kmOG) <- AgGeom
    st_crs(AgCensus5kmOG) <- 27700
    st_write(AgCensus5kmOG, "data_in/animals/AgCensus5kmOGm.gpkg", append = F)
    
    head(AgCensus5kmOG)
    ## only keep info of interest i.e. pig, horses, sheep, poultry, spatial refs
    agNames <- fread(file.path(currPath, "England_2010_5k.csv"))
    # get the full names
    xNames <- which(grepl(paste("id|pig|sheep|ram|ewe|non-bree|lamb|shee|fattening|poul|horse"
                                , "cattle", "beef", "dairy", sep = "|")
                          , names(agNames)))
    # as AgCensus has the additional id cols, remove those first
    AgCensus5km <- AgCensus5kmOG %>%
      dplyr::select(-c("left","top","right","bottom")) %>%
      # select those found above
      dplyr::select(c(id, all_of(xNames)))
    # assign colnames based on xNames
    colnames(AgCensus5km)[2:(ncol(AgCensus5km)-1)] <- (c(names(agNames)[xNames]))
    # make fids
    AgCensus5km$fid5km <- paste0("fid5km_", 1:nrow(AgCensus5km))
    head(AgCensus5km)
    # load in lcm plus crops for 2015 - the shapefile made earlier in the script
    gridName <- paste0("land_cover_table" ,".gpkg")
    # create centroids first
    grass2015 <- st_read(file.path("data_in", "land_cover", gridName)) %>% st_centroid() %>%
      # just keep grasses
      dplyr::select(rcFid_1km, improved_grass_ha
                    , neutral_grass_ha, calc_grass_ha, acid_grass_ha, heather_grass_ha) %>%
      # get total sum of grass
      mutate(grass_total_ha = rowSums(.[,2:(ncol(.)-1),drop=TRUE], na.rm = TRUE)) %>%
      dplyr::select(rcFid_1km, grass_total_ha)
    head(grass2015)
    
    grassyOut <- list() # save list
    
    if(file.exists(file.path(currPath, "nonBovine_grid_2015_1km.gpkg"))){
      nonBov1kmGB <- st_read(file.path(currPath, "nonBovine_grid_2015_1km.gpkg"))
    } else {
      ## Parallel code
      ## ------------ Notes --------------  ##
      ## should take 1h 20 mins
      ## ------------ ----- --------------  ##
      tic("total time for non-bovine grid")
      registerDoParallel(cl <- makeCluster(5))
      results_list <- foreach(i = 1:nrow(AgCensus5km)
                              , .combine = "rbind"           # Combine results
                              # Self-load
                              , .packages = c("sf", "dplyr")) %dopar% {
                                
                                # extract row
                                nonBovRow <- AgCensus5km[i, ]
                                # determine which 1 km cells each 5 km cell contains
                                grassRow <- st_filter(grass2015, nonBovRow)
                                # st_write(grassRow, "data_in/animals/grassRow.gpkg", append = F)
                                if(nrow(grassRow) > 0){
                                  # get grass in 1 km as a percentage of grass in 5 km2
                                  grassRow$grassPC5km <- grassRow$grass_total_ha / sum(grassRow$grass_total_ha)
                                  
                                  # duplicate rows of the cow grid to equal the 1 km underlying it
                                  rowAgain <- rbind(nonBovRow %>%
                                                      st_drop_geometry(), nonBovRow[rep(1, (nrow(grassRow)-1)), ]
                                                    %>% st_drop_geometry()) %>%
                                    dplyr::select(-c(fid5km, id))
                                  head(rowAgain)
                                  
                                  # multiply, to get the quantity of each bovine category in each 1 km2
                                  nonBovGrass <- grassRow$grassPC5km * rowAgain
                                  nonBovGrass <- bind_cols(grassRow %>% st_drop_geometry() %>% dplyr::select(rcFid_1km), nonBovGrass)
                                  
                                }
                              }
      stopCluster(cl)
      toc()
      
      ### bind them
      nonBov1km <- bind_rows(results_list)
      # merge with GB grid
      nonBov1kmGB <- merge(ghgGridGB, nonBov1km, by = "rcFid_1km", all = T)
      nonBov1kmGB[is.na(nonBov1kmGB)] <- 0
      head(nonBov1kmGB)
      class(nonBov1kmGB)
      
      #### save
      st_write(nonBov1kmGB, file.path(currPath, "nonBovine_grid_2015_1km.gpkg"), append = F)
      
      # write readme
      sink(file = NULL)
      sink(file = file.path(currPath, "readme_nonBovine_grid_2015_1km.md"))
      cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-08-10
  Last update:",  format(Sys.Date()), "
  
  Description of 'nonBovine_grid_2015_1km.gpkg':
  This is a gridded dataset that contains the quantity of different farmed non-bovine animals, based on the underlying grass area. 
  The animals are split up into different species and life stages e.g. breeding ewes intended for further breeding
  These data were obtained from 2015 AgCensus, which was at the 5 km2 resolution and featured data collected in 2010
  Using the amount of grass in each 1 km2 underlying the each 5 km2, the non-bovine numbers in a 5 km2 cell where proportied to each 1 km2 based on proportion of grass
  
  Citation for the original data: 
  
  
  Columns of 'nonBovine_grid_2015_1km.gpkg':
    'rcFid_1km' = centroid point of 1 km reference id (EPSG: 27700) [x and y coordinates are seperated by '_']
    '[animal_type]' = derived (using proportion of grass in 5 km2) numbers of animals of that type present within a 1 km2 cell for 2010
    
  Data: non-bovine numbers
  units: quantity
  Spatial resolution: 1 km2
  Spatial extent: GB
  Temporal coverage: 2015
  Native projection: 27700")
      sink(file = NULL)
    }
    
    ###### 2b4 - determine which km is in which region of England ######
    ## ---------- notes ---------- # 
    ## for the animals, emissions differ by region. Therefore, each cell needs to
    ## be assigned to a specific region. However, we only have that data for England
    ## so, all others will be labelled 'other' and matched up with the 'all' region
    ## category for the emissions 
    ## --------------------------- # 
    
    # load in regions of England
    regions <- st_read("N:/Data/UK/Boundary/region/RGN_DEC_2021_EN_BUC.shp")
    # get centroids of animal grids
    gridCentres <- nonBov1kmGB %>%
      st_centroid()
    # load in regions for Wales and Scotland
    regions.WS <- st_read("N:/Data/UK/Boundary/data/bdline_gb.gpkg"
                          , layer='scotland_and_wales_region')
    plot(regions.WS)
    # combine
    regions.UK <- bind_rows(regions, regions.WS) %>%
      mutate(fid = 1:n()) %>% relocate(fid) %>%
      dplyr::select(fid, RGN21NM, Name, Area_Code) 
    ## derive names for countries and counties
    regions.UK <- regions.UK %>%
      mutate(area = ifelse(is.na(Name), RGN21NM, Name)) %>%
      mutate(area = gsub(" PER", "", area))
    plot(regions.UK[6])
    
    # loop through the regions, extracting cell centroids, and then bind these back
    # together
    for(i in unique(regions.UK$area)){
      
      tic("using region")
      nm <- i
      # print
      cat("i =", i, "| name =", nm, "\n")
      
      goodPoints <- st_filter(gridCentres, regions.UK %>% filter(area == i)) %>%
        dplyr::select(rcFid_1km) %>%
        mutate(region := nm) %>%
        relocate(rcFid_1km, region)
      head(goodPoints)
      
      # if first iteration, create new df
      # else, bind with previously-created df
      if(i == unique(regions.UK$area)[1]){
        cellRegions <- goodPoints
      } else {
        cellRegions <- rbind(cellRegions, goodPoints)
      }
      toc()
    }
    head(cellRegions)
    plot(cellRegions[2])
    
    ###### 2b5 - proportional relationship of non-bovine animals between 2010 and 2020 ######
    ## ------------ Notes --------------  ##
    ## in this section, proportional relationships between non-bovine animals recorded in
    ## 2010 (i.e. same time as AgCensus data) and 2020 will be used to adjust the 
    ## original 2010 AgCensus numbers. This is only for England. For Scotland and Wales, 
    ## see the next sections (2b5b, 2b5c)
    
    ## Pigs decreased by 26% between 2010 and 2020, while sheep decreased 21%
    ## and poultry increased by 7%
    ## ------------ ----- --------------  ##
    
    ###### 2b5a - proportional relationship of non-bovine animals - England
    nonBov1kmGB2 <- merge(nonBov1kmGB 
                          , cellRegions %>% st_drop_geometry()
                          , all = T)  %>%
      # select only non-bovine
      dplyr::select(contains(c("Fid", "region"
                               # extract pigs
                               , "pigs"
                               # extract poultry
                               , "poultry"
                               # extract sheep
                               , "ewe", "sheep", "rams"
                               # extract horse
                               , "horse")))
    names(nonBov1kmGB2)
    head(nonBov1kmGB2)
    str(nonBov1kmGB2)
    plot(nonBov1kmGB2[2])        
    
    # for each animal multiply it by its proportion difference.
    # Note: horse does not change
    ## pigs
    nonBov1kmGB.prop2 <- nonBov1kmGB2 %>%
      mutate(across(contains("pigs"), ~.*0.74)
             ## sheep
             , across(contains(c("ewe", "sheep", "rams")), ~.*0.79)
             ## sheep
             , across(contains("poultry"), ~.*1.07))
    head(nonBov1kmGB.prop2)
    
    st_write(nonBov1kmGB.prop2, "nonBov1kmGBprop2.gpkg")
    
    ###### 2b5b - proportional relationship of non-bovine animals - Scotland and Wales
    ## ------------ Notes --------------  ##
    ## livestock data for Wales can be obtained in two ways: total livestock
    ## numbers in 2017, by region (https://statswales.gov.wales/Catalogue/Agriculture/Agricultural-Survey/Area-Survey-Results/total-livestock-in-wales-by-area)
    ## and 2020 livestock numbers, by animal type and 'role' (https://www.gov.wales/survey-agriculture-and-horticulture-june-2022)
    ## From the 2017 data, the proportion of animals in different regions of Wales could have been calculated in terms
    ## of percentage. However, to make the methodology consistent, they were calculated the same way as Scotland. See next 'Notes' section
    ## ------------ ----- --------------  ##
    
    ## ------------ Notes --------------  ##
    ## livestock data for Scotland were obtained by using 1 km2 published maps (in image form)
    ## from UK Government data. The original categories were ranges of:
    ## <=1, 1-5, 5 - 20, 20-80, 80-200, 200-430 for pigs
    ## <=1, 1-30, 30-60, 60-120, 120-240, 240-460 for sheep
    ## the top number of each bracket was used 
    ## These data were geoprocessed from the original, polygonised.
    ## Below, these polygons will rasterised 
    ## ------------ ----- --------------  ##
    stop("grass roots")
    
    # read in polygons
    pigNumb.poly <- st_read(file.path(currPath, "pig_numbers.gpkg"))
    pigNumb.poly
    st_bbox(pigNumb.poly); st_crs(pigNumb.poly)
    # convert to 27700
    pigNumb.poly27700 <- st_transform(pigNumb.poly, 27700)
    pigNumb.rast <- st_rasterize(pigNumb.poly27700 %>% dplyr::select(pigs_km)
                                 , dy = 1000, dx = 1000) %>%
      rast()
    pigNumb.rast
    plot(pigNumb.rast)
    # save
    writeRaster(pigNumb.rast, file.path(currPath, "pig_numbers.tif")
                , overwrite=TRUE)
    
    sheepNumb.poly <- st_read(file.path(currPath, "sheep_numbers.gpkg"))
    sheepNumb.poly
    st_bbox(sheepNumb.poly); st_crs(sheepNumb.poly)
    # convert to 27700
    sheepNumb.poly27700 <- st_transform(sheepNumb.poly, 27700)
    sheepNumb.rast <- st_rasterize(sheepNumb.poly27700 %>% 
                                     rename(sheep_per_km = sheep_numb) %>%
                                     dplyr::select(sheep_per_km)
                                   , dy = 1000, dx = 1000) %>%
      rast()
    sheepNumb.rast
    plot(sheepNumb.rast)
    # save
    writeRaster(sheepNumb.rast, file.path(currPath, "sheep_numbers.tif")
                , overwrite=TRUE)
    
    # extract to polygon
    nonBov1kmGBprop2 <- st_read("nonBov1kmGBprop2.gpkg")
    sheepNumb.point <- raster(sheepNumb.rast) %>% rasterToPoints() %>%
      # Convert SpatialPointsDataFrame to sf object
      as.data.frame() %>%
      st_as_sf(coords = c("x", "y"))
    st_crs(sheepNumb.point) <- st_crs(nonBov1kmGBprop2)
    sheepNumb.point <- st_intersection(nonBov1kmGBprop2, sheepNumb.point)
    plot(sheepNumb.point[1])
    sheepNumb.point <- sheepNumb.point %>% st_drop_geometry() %>% 
      dplyr::select(rcFid_1km, sheep_per_km) %>%
      filter(!is.na(rcFid_1km))
    nonBov1km <- dplyr::full_join(nonBov1kmGBprop2, sheepNumb.point
                                  , relationship = "many-to-many") %>%
      filter(!is.na(region))
    
    pigNumb.point <- raster(pigNumb.rast) %>% rasterToPoints() %>%
      # Convert SpatialPointsDataFrame to sf object
      as.data.frame() %>%
      st_as_sf(coords = c("x", "y"))
    st_crs(pigNumb.point) <- st_crs(nonBov1kmGBprop2)
    pigNumb.point <- st_intersection(nonBov1kmGBprop2, pigNumb.point)
    plot(pigNumb.point[1])
    pigNumb.point <- pigNumb.point %>% st_drop_geometry() %>% 
      dplyr::select(rcFid_1km, pigs_km) %>%
      rename(pigs_per_km = pigs_km) %>%
      filter(!is.na(rcFid_1km))
    nonBov1km <- dplyr::full_join(nonBov1km, pigNumb.point
                                  , relationship = "many-to-many") %>%
      filter(!is.na(region))
    
    ## use the proportions of different roles of different pigs and sheep to determine those in scotland
    nonBov1km.spl <- nonBov1kmGBprop2 %>% 
      st_drop_geometry() %>%
      filter(total.pigs..c47. > 0) %>%
      dplyr::select(contains(c("breeding.pig", "fattening.pig"))) %>%
      ### get totals, and then proportion of each caetgories
      mutate(total = rowSums(., na.rm = T))
    nonBov1km.s1 <- mean(nonBov1km.spl$breeding.pigs..incl..gilts.and.boars...c45. / nonBov1km.spl$total)
    ### breeding pigs make up 0.91 of all pigs
    
    ## use the proportions of different roles of different pigs and sheep to determine those in scotland
    nonBov1km.spl <- nonBov1kmGBprop2 %>% 
      st_drop_geometry() %>%
      filter(total.pigs..c47. > 0) %>%
      dplyr::select(contains(c("ewe", "sheep", "lamb", "ram"))) %>%
      ### get totals, and then proportion of each caetgories
      mutate(total = rowSums(., na.rm = T))
    nonBov1km.s1 <- mean(nonBov1km.spl$breeding.ewes.intended.for.further.breeding..c38. / nonBov1km.spl$total)
    nonBov1km.s2 <- mean(nonBov1km.spl$breeding.ewes.intended.for.slaughter..c39. / nonBov1km.spl$total)
    nonBov1km.s3 <- mean(nonBov1km.spl$ewes.intended.for.first.time.breeding..c40. / nonBov1km.spl$total)
    nonBov1km.s4 <- mean(nonBov1km.spl$non.breeding.sheep..1yr..male.and.female...c42. / nonBov1km.spl$total)
    nonBov1km.s5 <- mean(nonBov1km.spl$rams.for.service..c41. / nonBov1km.spl$total)
    nonBov1km.sum <- sum(nonBov1km.s1, nonBov1km.s2, nonBov1km.s3, nonBov1km.s4, nonBov1km.s5)
    nonBov1km.s1 <- nonBov1km.s1 / nonBov1km.sum
    nonBov1km.s2 <- nonBov1km.s2 / nonBov1km.sum
    nonBov1km.s3 <- nonBov1km.s3 / nonBov1km.sum
    nonBov1km.s4 <- nonBov1km.s4 / nonBov1km.sum
    nonBov1km.s5 <- nonBov1km.s5 / nonBov1km.sum
    ### breeding ewes make up 0.54 and 0.36, first time breeders are 0.03, non-breeding sheep are 0.01
    ### and rams for service are 0.06 
    
    nonBov1km.part2 <- nonBov1km %>%
      # deal with pig proportions
      mutate(pigs_per_km = ifelse(is.na(pigs_per_km), 0, pigs_per_km)) %>%
      mutate(breed_pigs_gilts_boars = ifelse(pigs_per_km > 1
                                             , pigs_per_km * 0.91
                                             , breeding.pigs..incl..gilts.and.boars...c45.)
             , fattening_pigs_piglets = ifelse(pigs_per_km > 1
                                               , pigs_per_km * 0.09
                                               , fattening.pigs.and.piglets..includes.barren.sows...c46.)
             , breed_pigs_gilts_boars = ifelse(pigs_per_km == 1, 1, breed_pigs_gilts_boars)) %>%
      # deal with sheep proportions
      mutate(sheep_per_km = ifelse(is.na(sheep_per_km), 0, sheep_per_km)) %>%
      mutate(breed_ewes_extra = ifelse(sheep_per_km > 1
                                       , sheep_per_km * 0.54
                                       , breeding.ewes.intended.for.further.breeding..c38.)
             , breed_ewes_slaughter = ifelse(sheep_per_km > 1
                                             , sheep_per_km * 0.36
                                             , breeding.ewes.intended.for.slaughter..c39.)
             , first_time_ewes = ifelse(sheep_per_km > 1
                                        , sheep_per_km * 0.03
                                        , ewes.intended.for.first.time.breeding..c40.)
             , non_breed_sheep = ifelse(sheep_per_km > 1
                                        , sheep_per_km * 0.01
                                        , non.breeding.sheep..1yr..male.and.female...c42.)
             , rams = ifelse(sheep_per_km > 1
                             , sheep_per_km * 0.06
                             , rams.for.service..c41.)
             , breed_ewes_extra = ifelse(sheep_per_km == 1, 1, breed_ewes_extra))
    
    st_write(nonBov1km.part2, "nonBov1km2.gpkg", append=TRUE)
    
    ###### 2b5c - proportional relationship of non-bovine animals - Poultry
    ## ------------ Notes --------------  ##
    ## according to Scotland Gov agricultural survey results [https://www.gov.scot/publications/results-december-2020-agricultural-survey/]
    ## there were 6.68 million poultry for meat production in 2020, and 6.6 million for egg production,
    ## which equals 13.28 together.
    
    ## for Wales, there were 2,438,761,  1,787,890, and 5,326,057 poultry for 
    ## North Wales, South Wales, Mid and West Wales, respectively. 
    
    ## like sheep and pigs, poultry was calculated using georeferenced maps
    ## ------------ ----- --------------  ##
    
    poultryNumb.poly <- st_read(file.path(currPath, "poultry_numbers.gpkg")) %>%
      dplyr::select(poul_numb) %>%
      mutate(poultry_per_km = as.numeric(poul_numb))
    poultryNumb.poly
    st_bbox(poultryNumb.poly); st_crs(poultryNumb.poly)
    # convert to 27700
    poultryNumb.poly27700 <- st_transform(poultryNumb.poly, 27700)
    poultryNumb.rast <- st_rasterize(poultryNumb.poly27700 %>% 
                                       dplyr::select(poultry_per_km)
                                     , dy = 1000, dx = 1000) %>%
      rast()
    poultryNumb.rast
    plot(poultryNumb.rast)
    # save
    writeRaster(poultryNumb.rast, file.path(currPath, "poultry_numbers.tif")
                , overwrite=TRUE)
    
    # convert to point data
    poultryNumb.point <- raster(poultryNumb.rast) %>% rasterToPoints() %>%
      # Convert SpatialPointsDataFrame to sf object
      as.data.frame() %>%
      st_as_sf(coords = c("x", "y"))
    st_crs(poultryNumb.point) <- st_crs(nonBov1kmGBprop2)
    poultryNumb.point <- st_intersection(nonBov1kmGBprop2, poultryNumb.point)
    plot(poultryNumb.point[1])
    poultryNumb.point <- poultryNumb.point %>% st_drop_geometry() %>% 
      dplyr::select(rcFid_1km, poultry_per_km) %>%
      filter(!is.na(rcFid_1km))
    nonBov1km <- dplyr::full_join(nonBov1km.part2, poultryNumb.point
                                  , relationship = "many-to-many") %>%
      filter(!is.na(region))
    
    nonBov1km.part3 <- nonBov1km %>%
      # deal with poultry
      mutate(poultry_per_km = ifelse(is.na(poultry_per_km), total.poultry..c48., poultry_per_km)) %>%
      # limit to 2500
      mutate(poultry_per_km = ifelse(poultry_per_km > 2500, 2501, poultry_per_km))
    
    st_write(nonBov1km.part3, "nonBov1km3.gpkg", append=TRUE)
    
  } # end of non-bovine grid
  
  ##### 2b6 - add non-bovine animals back in #####
  animalsFinal <- cows1kmGB %>% 
    merge(., nonBov1km.part3 %>% dplyr::select(rcFid_1km, region, breed_pigs_gilts_boars:poultry_per_km) %>%
            st_drop_geometry()
          , by = "rcFid_1km") %>%
    relocate(rcFid_1km, region)
  names(animalsFinal)
  
  # save 
  st_write(animalsFinal, file.path(currPath, "all_animals_1km.gpkg"), append = F)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(currPath, "readme_all_animals_1km.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-08-10
  Last update:",  format(Sys.Date()), "
  
  Description of 'all_animals_1km.gpkg':
  This is a gridded dataset that contains the quantity of different farmed bovine and non-bovine animals.
  The non-bovine animals were derived from spatially-aggregated government agricultural data from 2020
  The bovine animals are split up into different species and life stages. This data was not available for non-bovine animals
  
  Citation for the original data: 
  
  
  Columns of 'all_animals_1km.gpkg':
    'rcFid_1km' = centroid point of 1 km reference id (EPSG: 27700) [x and y coordinates are seperated by '_']
    '[animal_type]' = derived numbers of animals of that type present within a 1 km2 cell for 2020
    'region' = region of GB (sub-divided into differnt UK regions)
    
  Data: animal (bovine and non-bovine) numbers
  units: quantity
  Spatial resolution: 1 km2
  Spatial extent: GB
  Temporal coverage: 2019/2020
  Native projection: 27700")
  sink(file = NULL)
  
  # tidy
  rm(landCoversFor1km, lc1km, lc1kmSpatial, lc1kmWide, i, part2AnimalGrid, xEnd, xStart
     , animalsFinal, cowsGrid, goodPoints)
  
} # end of 'part2AnimalGrid'

#### 3 - part3AnimalAttributes ####
if(part3AnimalAttributes){
  
  # set location of files that relate to animal data
  currPath <- file.path("./data_in", "animals")
  
  # load in animal numbers
  animalNumbers <- st_read(file.path(currPath, "all_animals_1km.gpkg"))
  animalNames <- names(animalNumbers) <- sub("_sum", "", names(animalNumbers))
  head(animalNames)
  tail(animalNames)
  rm(animalNumbers)
  
  ##### 3a - categorise the animals #####
  ## ------------ Notes --------------  ##
  ## To match the structure of CFT, herd structure has 5 elements:
  ## Calves, heifers, dairy cows, nursing or suckling cows, and other cattle
  ## CFT split cows into breeds, however, that info is not available in CTS
  ## ------------ ----- --------------  ##
  # get names, and assign to appropriate category
  animalNms <- animalNames[2:(length(animalNames)-1)]
  animalNms
  # just cows
  cowNm <- as.data.frame(cbind(animal = animalNms[which(animalNms == "Continental_0...3.months_Breeding.Heifers"): which(animalNms == "Upland_9...12.months_Breeding.bulls")]
                               , category = ""))
  # sort each cow into the relevant CFT category
  cowNm$category <- ifelse(grepl("eifer", cowNm$animal), "heifer"
                           # calves are considered to be below 8 months old
                           , ifelse(grepl("0...3.months|3...6.months", cowNm$animal), "calf"
                                    , ifelse(grepl("0...3.months|3...6.months|calves", cowNm$animal), "calf"
                                             # other cattle
                                             , ifelse(grepl("bull|Bull|teer", cowNm$animal), "other_cattle"
                                                      # according to RSCPA (https://www.rspca.org.uk/adviceandwelfare/farm/dairy/farming)
                                                      # cows are impregnated around 15 months, to around 24 months
                                                      , ifelse(grepl("15...18|18...21|21...24|In.calf", cowNm$animal), "pregnant"
                                                               # cows are impregnated around 15 months, to around 24 months
                                                               , ifelse(grepl("24...27|27...30|30...33|airy", cowNm$animal), "lactating_dairy"
                                                                        , "other_cattle"))))))
  
  # non-bovine animal names
  nonBovineNm <- animalNames[!animalNames %in% c(cowNm$animal
                                                 , "rcFid_1km"
                                                 , "geom")]
  nonBovineNm <- nonBovineNm[!grepl("region", nonBovineNm)]
  cat(nonBovineNm, sep = "\n")
  # not cows
  nonBovineNm <- as.data.frame(cbind(animal = nonBovineNm
                                     , category = ""))
  nonBovineNm$category <- ifelse(grepl("ram", nonBovineNm$animal), "rams"
                                 , ifelse(grepl("ewe", nonBovineNm$animal), "ewes"
                                          , ifelse(grepl("lamb", nonBovineNm$animal), "lambs"
                                                   , ifelse(grepl("oultry|aying|urkey|hicken", nonBovineNm$animal), "poultry"
                                                            , ifelse(grepl("orse", nonBovineNm$animal), "horse"
                                                                     , ifelse(grepl("pig", nonBovineNm$animal), "pigs"
                                                                              , "other_sheep"))))))
  print(nonBovineNm)
  # combine, to get all animals
  names(cowNm) <- names(nonBovineNm)
  animalDf <- bind_rows(cowNm, nonBovineNm)
  rm(animalNms, cowNm, nonBovineNm)
  gc()
  
  ##### 3b - create animal attributes table #####
  ###### 3b1 - animal weights ######
  ## ------------ Notes --------------  ##
  ## weights come from Analysis of Farm Business Survey (FBS)
  ## module data on livestock productivity, practices and 
  ## greenhouse gas emissions. - DO0148 (Parson and Williams (2015))
  ## https://randd.defra.gov.uk/ProjectDetails?ProjectID=19397&FromSearch=Y&Status=3&Publisher=1&SearchText=DO0148&SortString=ProjectCode&SortOrder=Asc&Paging=10#Description
  ## ------------ ----- --------------  ##
  
  ## dairy cows
  animalWeights <- animalDf[animalDf$category == "lactating_dairy", ]
  animalUnique <- c(unique(animalWeights$animal))
  
  # weights significantly vary by region for dairy (Table 35 in PW2015)
  # for medium cows, use average, else use average -/+ 1/2 SD
  for(i in 1:length(animalUnique)){
    if(i == 1){
      # create empty list
      dairyWeight <- list()
    }
    
    if(grepl("edium", animalUnique[[i]])){
      dairyWeight[[i]] <- as.data.frame(
        cbind(animal = animalUnique[[i]]
              , East_Midlands = 587
              , East_of_England = 603
              , North_East = 595
              , North_West = 583
              , South_East = 541
              , South_West = 551
              , West_Midlands = 610
              , Yorkshire_and_the_Humber = 648)) %>% rowwise() %>%
        mutate(other = as.character(mean(as.numeric(.[1, ]), na.rm = T)))
    } else if(grepl("arge", animalUnique[[i]])){
      dairyWeight[[i]] <- as.data.frame(
        cbind(animal = animalUnique[[i]]
              , East_Midlands = 587 + (62/2)
              , East_of_England = 603 + (62/2)
              , North_East = 595 + (47/2)
              , North_West = 583 + (59/2)
              , South_East = 541 + (64/2)
              , South_West = 551 + (54/2)
              , West_Midlands = 610 + (69/2)
              , Yorkshire_and_the_Humber = 648 + (53/2))) %>% rowwise() %>%
        mutate(other = as.character(mean(as.numeric(.[1, ]), na.rm = T)))
    } else {
      dairyWeight[[i]] <- as.data.frame(
        cbind(animal = animalUnique[[i]]
              , East_Midlands = 587 - (62/2)
              , East_of_England = 603 - (62/2)
              , North_East = 595 - (47/2)
              , North_West = 583 - (59/2)
              , South_East = 541 - (64/2)
              , South_West = 551 - (54/2)
              , West_Midlands = 610 - (69/2)
              , Yorkshire_and_the_Humber = 648 - (53/2))) %>% rowwise() %>%
        mutate(other = as.character(mean(as.numeric(.[1, ]), na.rm = T)))
    }
    
    if(i == length(animalUnique)){
      # combine outputs
      dairyWeight <- bind_rows(dairyWeight) %>%
        pivot_longer(cols = !animal, names_to = "UK.Region", values_to = "Weight_Kg")
      # tidy 
      rm(animalUnique, animalWeights, i)
      head(dairyWeight)
    }
  }
  
  ## dairy heifers
  animalWeights <- animalDf[grepl("In.calf.heifers", animalDf$animal)| animalDf$category == "pregnant", ]
  animalUnique <- c(unique(animalWeights$animal))
  # weights significantly vary by region for dairy (Table 37 in PW2015)
  # for medium cows, use average, else use average -/+ 1/2 SD
  for(i in 1:length(animalUnique)){
    if(i == 1){
      # create empty list
      dairyHefWeight <- list()
    }
    
    if(grepl("edium", animalUnique[[i]])){
      dairyHefWeight[[i]] <- as.data.frame(
        cbind(animal = animalUnique[[i]]
              , East_Midlands = 514
              , East_of_England = 522
              , North_East = 525
              , North_West = 535
              , South_East = 497
              , South_West = 512
              , West_Midlands = 522
              , Yorkshire_and_the_Humber = 594)) %>% rowwise() %>%
        mutate(other = as.character(mean(as.numeric(.[1, ]), na.rm = T)))
    } else if(grepl("arge", animalUnique[[i]])){
      dairyHefWeight[[i]] <- as.data.frame(
        cbind(animal = animalUnique[[i]]
              , East_Midlands = 514 + (74/2)
              , East_of_England = 522 + (71/2)
              , North_East = 525 + (62/2)
              , North_West = 535 + (56/2)
              , South_East = 497 + (72/2)
              , South_West = 512 + (54/2)
              , West_Midlands = 522 + (60/2)
              , Yorkshire_and_the_Humber = 594 + (52/2))) %>% rowwise() %>%
        mutate(other = as.character(mean(as.numeric(.[1, ]), na.rm = T)))
    } else {
      dairyHefWeight[[i]] <- as.data.frame(
        cbind(animal = animalUnique[[i]]
              , East_Midlands = 514 - (74/2)
              , East_of_England = 522 - (71/2)
              , North_East = 525 - (62/2)
              , North_West = 535 - (56/2)
              , South_East = 497 - (72/2)
              , South_West = 512 - (54/2)
              , West_Midlands = 522 - (60/2)
              , Yorkshire_and_the_Humber = 594 - (52/2))) %>% rowwise() %>%
        mutate(other = as.character(mean(as.numeric(.[1, ]), na.rm = T)))
    }
    
    if(i == length(animalUnique)){
      # combine outputs
      dairyHefWeight <- bind_rows(dairyHefWeight) %>%
        pivot_longer(cols = !animal, names_to = "UK.Region", values_to = "Weight_Kg")
      # tidy 
      rm(animalUnique, animalWeights, i)
    }
  }
  head(dairyHefWeight)
  
  ## dairy calves
  animalWeights <- animalDf[animalDf$category == "calf" & grepl("airy", animalDf$animal), ]
  # Table 38 in PW2015
  # for medium cows, use average, else use average -/+ 1/2 SD
  animalWeights$Weight_Kg <- ifelse(grepl("edium", animalWeights$animal), 50
                                    , ifelse(grepl("arge", animalWeights$animal), 50 + (8.5/2)
                                             , 50 - (8.5/2)))
  # set regions
  animalWeights$UK.Region <- "All"
  # correct order
  dairyCalfWeight <- animalWeights %>%
    select(colnames(dairyWeight))
  # tidy
  rm(animalWeights)
  
  ## 3a4 beef cattle - LFA (i.e. non-lowland)
  animalWeights <- animalDf[animalDf$category == "other_cattle"
                            & !grepl("Lowland", animalDf$animal) 
                            & grepl("6...9|9...12|12...15|15...18|18...21|21...24", animalDf$animal),  ]
  # according to PW2015, beef weights did not significantly vary per region, 
  # so one value can be used for each
  # does separate LFA and lowland however
  # only animals under two years
  # Table 40 in PW2015
  # for non-lowland cows use average, else use average -/+ 1/2 SD based on ages
  animalWeights$Weight_Kg <- ifelse(grepl("6...9", animalWeights$animal), 587 - (66/2)
                                    , ifelse(grepl("21...24", animalWeights$animal), 587 + (66/2)
                                             , 587))
  # set regions
  animalWeights$UK.Region <- "All"
  # correct order
  beefLFAWeight <- animalWeights %>%
    select(colnames(dairyWeight))
  # tidy
  rm(animalWeights)
  
  ## beef cattle - lowland
  animalWeights <- animalDf[animalDf$category == "other_cattle"
                            & grepl("Lowland", animalDf$animal) 
                            & grepl("6...9|9...12|12...15|15...18|18...21|21...24", animalDf$animal),  ]
  # Table 41 in PW2015
  # for lowland cows use average, else use average -/+ 1/2 SD based on ages
  animalWeights$Weight_Kg <- ifelse(grepl("6...9", animalWeights$animal), 594 - (86/2)
                                    , ifelse(grepl("21...24", animalWeights$animal), 594 + (86/2)
                                             , 594))
  # set regions
  animalWeights$UK.Region <- "All"
  # correct order
  beefLowlandWeight <- animalWeights %>%
    select(colnames(dairyWeight))
  # tidy
  rm(animalWeights)
  
  ## beef cattle - lowland
  animalWeights <- animalDf[animalDf$category == "heifer" 
                            & !grepl("airy", animalDf$animal),]
  
  # Table 42 in PW2015
  # for beef heifer cows use average, else use average -/+ 1/2 SD based on ages (i.e.
  # assuming young ones are not in-calf)
  animalWeights$Weight_Kg <- ifelse(grepl("0...3|3...6", animalWeights$animal), 554 - (78/2)
                                    , ifelse(grepl("24...27|27...30|30...33|33...36", animalWeights$animal), 554 + (78/2)
                                             , 554))
  # set regions
  animalWeights$UK.Region <- "All"
  # correct order
  beefHeiferWeight <- animalWeights %>%
    select(colnames(dairyWeight))
  # tidy
  rm(animalWeights)
  
  ## beef cattle - 2+ years, male
  animalWeights <- animalDf[animalDf$category == "other_cattle"
                            & grepl("ull|teer", animalDf$animal) 
                            & grepl("24...27|27...30|30...33|33...36|36...48|48...60|240", animalDf$animal),  ]
  # Table 43 in PW2015
  # for male 2+ cows use average, else use average -/+ 1/2 SD based on ages 
  animalWeights$Weight_Kg <- ifelse(grepl("24...27", animalWeights$animal), 596 - (91/2)
                                    , ifelse(grepl("48...60", animalWeights$animal), 596 + (91/2)
                                             , 596))
  # set regions
  animalWeights$UK.Region <- "All"
  # correct order
  beef2YrMaleWeight <- animalWeights %>%
    select(colnames(dairyWeight))
  # tidy
  rm(animalWeights)
  
  ## beef cattle - 2+ years, female
  animalWeights <- animalDf[animalDf$category == "other_cattle"
                            & grepl("cow|emale", animalDf$animal) 
                            & grepl("24...27|27...30|30...33|33...36|36...48|48...60|240", animalDf$animal),  ]
  animalUnique <- c(unique(animalWeights$animal))
  # female weights significantly vary by region (Table 45 in PW2015)
  # for medium cows, use average, else use average -/+ 1/2 SD
  for(i in 1:length(animalUnique)){
    if(i == 1){
      # create empty list
      beef2YrFemaleWeight <- list()
    }
    
    if(grepl("27...30|30...33|33...36|36...48", animalUnique[[i]])){
      beef2YrFemaleWeight[[i]] <- as.data.frame(
        cbind(animal = animalUnique[[i]]
              , East_Midlands = 438
              , East_of_England = 451
              , North_East = 455
              , North_West = 447
              , South_East = 480
              , South_West = 443
              , West_Midlands = 465
              , Yorkshire_and_the_Humber = 511)) %>% rowwise() %>%
        mutate(other = as.character(mean(as.numeric(.[1, ]), na.rm = T)))
    } else if(grepl("48...60", animalUnique[[i]])){
      beef2YrFemaleWeight[[i]] <- as.data.frame(
        cbind(animal = animalUnique[[i]]
              , East_Midlands = 438 + (84/2)
              , East_of_England = 451 + (104/2)
              , North_East = 455 + (99/2)
              , North_West = 447 + (79/2)
              , South_East = 480 + (86/2)
              , South_West = 443 + (82/2)
              , West_Midlands = 465 + (72/2)
              , Yorkshire_and_the_Humber = 511 + (85/2))) %>% rowwise() %>%
        mutate(other = as.character(mean(as.numeric(.[1, ]), na.rm = T)))
    } else {
      beef2YrFemaleWeight[[i]] <- as.data.frame(
        cbind(animal = animalUnique[[i]]
              , East_Midlands = 438 - (84/2)
              , East_of_England = 451 - (104/2)
              , North_East = 455 - (99/2)
              , North_West = 447 - (79/2)
              , South_East = 480 - (86/2)
              , South_West = 443 - (82/2)
              , West_Midlands = 465 - (72/2)
              , Yorkshire_and_the_Humber = 511 - (85/2))) %>% rowwise() %>%
        mutate(other = as.character(mean(as.numeric(.[1, ]), na.rm = T)))
    }
    
    if(i == length(animalUnique)){
      # combine outputs
      beef2YrFemaleWeight <- bind_rows(beef2YrFemaleWeight) %>%
        pivot_longer(cols = !animal, names_to = "UK.Region", values_to = "Weight_Kg")
      # tidy 
      rm(animalUnique, animalWeights, i)
    }
  }
  
  ## other calves
  animalWeights <- animalDf[animalDf$category == "calf"
                            & !grepl("airy", animalDf$animal),  ]
  animalUnique <- c(unique(animalWeights$animal))
  # calf weights significantly vary by region (Table 51 in PW2015)
  for(i in 1:length(animalUnique)){
    if(i == 1){
      # create empty list
      beefCalfWeight <- list()
    }
    
    if(grepl("3...6.months", animalUnique[[i]])){
      beefCalfWeight[[i]] <- as.data.frame(
        cbind(animal = animalUnique[[i]]
              , East_Midlands = 53.3
              , East_of_England = 45.5
              , North_East = 47.5
              , North_West = 48.5
              , South_East = 49
              , South_West = 51.1
              , West_Midlands = 53.4
              , Yorkshire_and_the_Humber = 42.3)) %>% rowwise() %>%
        mutate(other = as.character(mean(as.numeric(.[1, ]), na.rm = T)))
    } else {
      beefCalfWeight[[i]] <- as.data.frame(
        cbind(animal = animalUnique[[i]]
              , East_Midlands = 53.3 - (7.1/2)
              , East_of_England = 45.5 - (6.8/2)
              , North_East = 47.5 - (3.8/2)
              , North_West = 48.5 - (7.2/2)
              , South_East = 49 - (6.3/2)
              , South_West = 51.1 - (8.5/2)
              , West_Midlands = 53.4 - (11.3/2)
              , Yorkshire_and_the_Humber = 42.3 - (8/2))) %>% rowwise() %>%
        mutate(other = as.character(mean(as.numeric(.[1, ]), na.rm = T)))
    }
    
    if(i == length(animalUnique)){
      # combine outputs
      beefCalfWeight <- bind_rows(beefCalfWeight) %>%
        pivot_longer(cols = !animal, names_to = "UK.Region", values_to = "Weight_Kg")
      # tidy 
      rm(animalUnique, animalWeights, i)
    }
  }
  
  # combine
  x <- unique(rbind(beef2YrFemaleWeight, beef2YrMaleWeight, beefHeiferWeight, beefLFAWeight
                    , beefLowlandWeight, dairyCalfWeight, dairyHefWeight, dairyWeight, beefCalfWeight)$animal)
  # combine
  cowWeights <- rbind(beef2YrFemaleWeight, beef2YrMaleWeight, beefHeiferWeight, beefLFAWeight
                      , beefLowlandWeight, dairyCalfWeight, dairyHefWeight, dairyWeight, beefCalfWeight)
  head(cowWeights)
  # tidy
  rm(beef2YrFemaleWeight, beef2YrMaleWeight, beefHeiferWeight, beefLFAWeight
     , beefLowlandWeight, dairyCalfWeight, dairyHefWeight, dairyWeight, beefCalfWeight
     , x)
  
  ## non-bovine weights ##
  ## rams
  animalWeights <- animalDf[animalDf$category == "rams", ]
  animalUnique <- c(unique(animalWeights$animal))
  # ram weights significantly vary by region (Table 55 in PW2015)
  for(i in 1:length(animalUnique)){
    ramWeight <- as.data.frame(
      cbind(animal = animalUnique[[i]]
            , East_Midlands = 71
            , East_of_England = 57
            , North_East = 66
            , North_West = 68
            , South_East = 55
            , South_West = 72
            , West_Midlands = 72
            , Yorkshire_and_the_Humber = 83)) %>% rowwise() %>%
      mutate(other = as.character(mean(as.numeric(.[1, ]), na.rm = T)))
    if(i == length(animalUnique)){
      # combine outputs
      ramWeight <- bind_rows(ramWeight) %>%
        pivot_longer(cols = !animal, names_to = "UK.Region", values_to = "Weight_Kg")
    }
  }
  # tidy 
  rm(animalUnique, i)
  head(ramWeight)
  
  ## ewes
  animalWeights <- animalDf[animalDf$category == "ewes", ]
  animalUnique <- c(unique(animalWeights$animal))
  # ewe weights significantly vary by region
  # (first time Table 61 in PW2015) and others Table 59
  for(i in 1:length(animalUnique)){
    if(i == 1){
      eweWeight <- list()
    }
    
    if(grepl("older", animalUnique[[i]])){
      eweWeight[[i]] <- as.data.frame(
        cbind(animal = animalUnique[[i]]
              , East_Midlands = 63
              , East_of_England = 46
              , North_East = 52
              , North_West = 60
              , South_East = 57
              , South_West = 63
              , West_Midlands = 62
              , Yorkshire_and_the_Humber = 69)) %>% rowwise() %>%
        mutate(other = as.character(mean(as.numeric(.[1, ]), na.rm = T)))
    } else if(grepl("first", animalUnique[[i]])){
      eweWeight[[i]] <- as.data.frame(
        cbind(animal = animalUnique[[i]]
              , East_Midlands = 48
              , East_of_England = 46
              , North_East = 41
              , North_West = 44
              , South_East = 46
              , South_West = 54
              , West_Midlands = 54
              , Yorkshire_and_the_Humber = 52)) %>% rowwise() %>%
        mutate(other = as.character(mean(as.numeric(.[1, ]), na.rm = T)))
    } else {
      eweWeight[[i]] <- as.data.frame(
        cbind(animal = animalUnique[[i]]
              , East_Midlands = 63 + (11/2)
              , East_of_England = 46 + (8/2)
              , North_East = 52 + (6/2)
              , North_West = 60 + (13/2)
              , South_East = 57 + (8/2)
              , South_West = 63 + (12/2)
              , West_Midlands = 62 + (12/2)
              , Yorkshire_and_the_Humber = 69 + (11/2))) %>% rowwise() %>%
        mutate(other = as.character(mean(as.numeric(.[1, ]), na.rm = T)))
    }
    
    if(i == length(animalUnique)){
      # combine outputs
      eweWeight <- bind_rows(eweWeight) %>%
        pivot_longer(cols = !animal, names_to = "UK.Region", values_to = "Weight_Kg")
      # tidy 
      rm(animalUnique, i)
    }
  }
  
  ## lambs ('other_saheep')
  unique(animalDf$category)
  animalWeights <- animalDf[animalDf$category == "other_sheep", ]
  animalUnique <- c(unique(animalWeights$animal))
  # lamb weights significantly vary by region (Table 68 in PW2015)
  for(i in 1:length(animalUnique)){
    if(i == 1){
      lambWeight <- list()
    }
    lambWeight[[i]] <- as.data.frame(
      cbind(animal = animalUnique[[i]]
            , East_Midlands = 38
            , East_of_England = 25
            , North_East = 37
            , North_West = 36
            , South_East = 32
            , South_West = 35
            , West_Midlands = 26
            , Yorkshire_and_the_Humber = 37)) %>% rowwise() %>%
      mutate(other = as.character(mean(as.numeric(.[1, ]), na.rm = T)))
    print(lambWeight[[i]])
    if(i == length(animalUnique)){
      # combine outputs
      lambWeight <- bind_rows(lambWeight) %>%
        pivot_longer(cols = !animal, names_to = "UK.Region", values_to = "Weight_Kg")
      # tidy 
      rm(animalUnique, i)
    }
  }
  
  ## pigs
  animalWeights <- animalDf[animalDf$category == "pigs", ]
  animalUnique <- c(unique(animalWeights$animal))
  # pig weight data from UK average clean pig carcase weights
  pigWeight <- as.data.frame(rbind(cbind(animal = animalUnique[[1]], UK.Region = "All", Weight_Kg = 83.3)
                                   , cbind(animal = animalUnique[[2]], UK.Region = "All", Weight_Kg = 83.3)
  ))
  
  ## poultry
  animalWeights <- animalDf[animalDf$category == "poultry", ]
  animalUnique <- c(unique(animalWeights$animal))
  # liveweights obtained from Defra (from Oct 2021)
  # https://www.gov.uk/government/statistics/historical-statistics-notices-on-poultry-and-poultry-meat-production-2021/united-kingdom-poultry-and-poultry-meat-statistics-december-2021
  # all poultry generally are considered broilers
  poultryWeight <- as.data.frame(rbind(cbind(animal = animalUnique[[1]], UK.Region = "All", Weight_Kg = 2.3)))
  
  ###### 3b2 - combine all animal weights ######
  # combine
  # merge all
  cspWeights <- Reduce(function(x, y) merge(x, y, all=TRUE)
                       , list(eweWeight, lambWeight
                              # , otherSheepWeight
                              , ramWeight, pigWeight
                              , poultryWeight)) %>%
    # bind with cows
    bind_rows(cowWeights, .) %>%
    mutate(Weight_Kg = as.numeric(Weight_Kg))
  
  # merge to get category names as well
  cspWeights <- rbind(merge(cspWeights, animalDf
                            , by = "animal"
                            , all = F)) %>%
    mutate(Weight_Kg = as.numeric(Weight_Kg))
  head(cspWeights)
  
  # tidy
  rm(cowWeights, eweWeight, lambWeight, ramWeight, pigWeight
     , animalUnique, animalWeights, poultryWeight)
  
  ###### 3b3 - Net and gross energies (NE / GE) ######
  ###### Net energy maintenance (NEm) ######
  ## ------------ Notes --------------  ##
  ## Cfi is the net energy maintenance coefficient
  ## Values taken from IPCC Table 10.4 for cows and sheep.
  ## There are no values for pigs, horses, or poultry. 
  ## Therefore, an assumption was made for horses and poultry that they all 
  ## had the average (diff between highest and lowest) of the
  ## other animals: 0.236, 0.315, 0.217, 0.37, 0.386, 0.322
  ## average: 0.2935
  ## Pigs, due to their similar size, were considered to be the same
  ## as adult sheep
  ## ------------ ----- --------------  ##
  # Net energy maintenance coefficient (cfi)
  cspParameters <- cspWeights %>%
    mutate(Cfi = if_else(grepl("lamb|non_breed_sheep", animal), 0.236
                         , if_else(grepl("sheep|ewe|pig|ram", animal), 0.217
                                   , if_else(grepl("ull", animal), 0.370
                                             , if_else(category == "lactating_dairy", 0.386
                                                       , if_else(grepl("teer|ther_|emale|eifer|calv|calf|cow", animal), 0.322, 0.2935)))))) %>%
    dplyr::relocate(animal, category, UK.Region)
  head(cspParameters)
  
  ## ------------ Notes --------------  ##
  ## NEm only requires weight and an animal coefficient,
  ## both of which are in the parameter table (cspParameters)
  ## average weights were derived as above
  ## eq. 10.3 IPCC 
  ## unit: MJ day-1
  ## cfi is a coefficient to determine net energy for maintenance
  ## ------------ ----- --------------  ##
  cspOutputs <- cspParameters %>%
    mutate(NEm = Weight_Kg ^ 0.75 * Cfi) %>%
    select(-c(Weight_Kg, Cfi)) %>%
    dplyr::relocate(animal, category, UK.Region)
  head(cspOutputs)
  
  # ensure both parameter and output dataset are in the same order
  cspOutputs <- cspOutputs[with(cspOutputs, order(animal, category, UK.Region)),]
  cspParameters <- cspParameters[with(cspParameters, order(animal, category, UK.Region)),]
  stopifnot(cspOutputs[, c(1:3)] == cspParameters[, c(1:3)])
  
  ###### Net energy activity (NEa) #####
  ## ------------ Notes --------------  ##
  ## Ca = coefficient corresponding to animal’s feeding situation
  ## lambs are always housed
  ## Table 10.5 in IPCC
  ## NEa unit: MJ day-1
  ## cattle and sheep use different equations: eq. 10.4/5 IPCC (cow/ sheep)
  ## different values for animal situation - so get all values
  ## IPCC report assumes that ‘[f]or poultry and swine, the feeding situation is
  ## assumed to be under confinement conditions and consequently the activity
  ## coefficient (Ca) is assumed to be zero as under these conditions very little
  ## energy is expended in acquiring feed.’
  ## ------------ Notes --------------  ##
  cspParameters <- cspParameters %>%
    mutate(Ca_housed = if_else(grepl("lamb|non_breed_sheep", animal), 0.0067
                               , if_else(grepl("oultry", animal), 0
                                         , if_else(grepl("pig", category), 0
                                                   , if_else(grepl("sheep|ewe|ram", animal), 0.0096
                                                             , if_else(grepl("other_cattle|heifer|calf|pregnant|dairy", category), 0 
                                                                       , -999)))))
           , Ca_flat_sml_past = if_else(grepl("lamb|non_breed_sheep", animal), -999
                                        , if_else(grepl("other_cattle|heifer|calf|pregnant|dairy", category), 0.17    
                                                  , if_else(grepl("sheep|ewe|ram", category), 0.0107
                                                            , if_else(grepl("poultry|pig", category), -999
                                                                      , -999))))
           , Ca_hill_lrg_past = if_else(grepl("lamb|non_breed_sheep", animal), -999
                                        , if_else(grepl("other_cattle|heifer|calf|pregnant|dairy", category), 0.36    
                                                  , if_else(grepl("sheep|ewe|ram", category), 0.024
                                                            , if_else(grepl("poultry|pig", category), -999
                                                                      , -999))))
    )
  # convert -999 to NA
  cspParameters[cspParameters == -999] <- NA
  head(cspParameters)
  
  ## NEa equations ##
  cspOutputs$NEa_housed <- cspOutputs$NEa_flat_sml <- cspOutputs$NEa_hill_lrg <- 0
  
  ###### Net energy growth (NEg) #####
  ## ------------ Notes --------------  ##
  ## eq. 10.6/7 IPCC (cows/sheep)
  ## NEg unit: MJ day-1
  ## three additional parameters are required here: 
  ## average daily weight gain (kg day-1)
  ## mature body weight (kg). 
  ## (For non-adult animals, target weight related to stage of growth is mature weight)
  ## Coefficient based on sex of animals (Csex)
  ## ------------ Notes --------------  ##
  
  ###### mature animal weight
  cspParameters <- cspParameters %>%
    # for cattle (adults, heifers, calves), target weights were taken from https://projectblue.blob.core.windows.net/media/Default/Beef%20&%20Lamb/BR_FeedingGrowingFinishingCattle-WEB.pdf
    mutate(mature_Kg = if_else(grepl("heifer", category), 685
                               , if_else(grepl("other_cattle|pregnant|dairy", category), 685
                                         # for calves, target weights can be split into month stages
                                         # https://ahdb.org.uk/knowledge-library/monitoring-dairy-calf-growth
                                         , if_else(grepl("calf", category) & grepl("0...3.|Small", animal), 116 
                                                   , if_else(grepl("calf", category) & grepl("3...6.|Med", animal), 185 
                                                             , if_else(grepl("calf", category), 200 
                                                                       # ewes taken from average: https://ahdb.org.uk/knowledge-library/selecting-ewe-lambs-for-breeding
                                                                       , if_else(grepl("ewe", animal), mean(c(60,70,80,65,80,65))
                                                                                 , if_else(grepl("lamb|non_breed_sheep", animal), 40
                                                                                           , if_else(grepl("ram", animal), 90
                                                                                                     # pigs based on average https://ahdb.org.uk/knowledge-library/feeding-growing-and-finishing-pigs
                                                                                                     , if_else(grepl("pig", animal), mean(c(93, 78))
                                                                                                               , if_else(grepl("oultry", animal), 2.5
                                                                                                                         , if_else(grepl("orse", animal), 500
                                                                                                                                   , -999))))))))))))
  head(cspParameters)
  head(cspParameters[order(cspParameters$mature_Kg, decreasing = F), ])
  
  ###### daily weight gain (kg) 
  cspParameters <- cspParameters %>%
    # for lamb, Target growth rates section in https://projectblue.blob.core.windows.net/media/Default/Beef%20&%20Lamb/GrowingAndFinishingLambsForBR3340_200415_WEB-1.pdf
    mutate(wGain_Kg = if_else(grepl("lamb|non_breed_sheep", animal), 0.25
                              # https://projectblue.blob.core.windows.net/media/Default/Beef%20&%20Lamb/BR_FeedingGrowingFinishingCattle-WEB.pdf
                              # Table 3
                              , if_else(grepl("heifer|other_cattle|pregnant|dairy", category), 1
                                        # for calves, target weights can be split into month stages
                                        , if_else(grepl("calf", category) & grepl("0...3.|Small", animal), .83 
                                                  , if_else(grepl("calf", category) & grepl("3...6.|Med", animal), .77
                                                            , if_else(grepl("calf", category), .80
                                                                      # https://projectblue.blob.core.windows.net/media/Default/Pork/Pork%20MI%20files/uk-pig-facts-and-figures_2506_190507_web.pdf
                                                                      # AHDB’s Pig Facts (Table 4.4) - assuming breeding pigs are ‘rearing’ and fattening pigs are ‘finishing’
                                                                      , if_else(grepl("pig", animal), c(mean(.866, .469))
                                                                                # Gibbs 2002 says sufficient to only include lamb weight gain
                                                                                , if_else(grepl("ewe", animal), .25
                                                                                          , if_else(grepl("sheep|ram", animal), .115
                                                                                                    # https://www.rspca.org.uk/webContent/staticImages/BroilerCampaign/EatSitSufferRepeat.pdf
                                                                                                    , if_else(grepl("oultry", animal), 0.063
                                                                                                              # DOI:10.5713/ajas.2006.86
                                                                                                              , 0))))))))))
  head(cspParameters)
  head(cspParameters[order(cspParameters$wGain_Kg, decreasing = F), ])
  
  for(i in 1:nrow(cspOutputs)){
    # print(i)
    
    # for cattle
    # Eq 10.4 IPCC: NEa = Ca * NEm
    if(grepl("other_cattle|heifer|calf|pregnant|dairy", cspOutputs[i, ]$category)){  
      cspOutputs[i, ]$NEa_housed <- cspParameters[i, ]$Ca_housed * cspOutputs[i, ]$NEm 
      cspOutputs[i, ]$NEa_flat_sml <- cspParameters[i, ]$Ca_flat_sml_past * cspOutputs[i, ]$NEm 
      cspOutputs[i, ]$NEa_hill_lrg <- cspParameters[i, ]$Ca_hill_lrg_past * cspOutputs[i, ]$NEm 
      # for sheep/ goats
      # Eq 10.5 IPCC: NEa = Ca * weight
    } else if (grepl("sheep|goat|ewe|ram", cspOutputs[i, ]$category)){
      cspOutputs[i, ]$NEa_housed <- cspParameters[i, ]$Ca_housed * cspParameters[i, ]$Weight_Kg 
      cspOutputs[i, ]$NEa_flat_sml <- cspParameters[i, ]$Ca_flat_sml_past * cspParameters[i, ]$Weight_Kg 
      cspOutputs[i, ]$NEa_hill_lrg <- cspParameters[i, ]$Ca_hill_lrg_past * cspParameters[i, ]$Weight_Kg 
      # IPCC report assumes that ‘[f]or poultry and swine, the feeding situation is
      # assumed to be under confinement conditions and consequently the activity
      # coefficient (Ca) is assumed to be zero as under these conditions very little
      # energy is expended in acquiring feed.’
    } else {
      cspOutputs[i, ]$NEa_housed <- 0
      cspOutputs[i, ]$NEa_flat_sml <- 0
      cspOutputs[i, ]$NEa_hill_lrg <- 0
    }
  }
  head(cspOutputs)
  
  ###### sex coefficients 
  # Csex. 0.8 for females, 1.2 for bulls (only required for cattle)
  # note: unsure which would be castrated, but they would have a coefficient of 1
  cspCsex <- cspParameters %>%
    filter(category %in% c("other_cattle", "calf", "heifer", "pregnant", "lactating_dairy")) %>%
    mutate(Csex = if_else(grepl("ull|teer", animal), 1.2, 0.8))
  
  # The NEg equation used for sheep includes two empirical constants (a and b) that vary by animal species/category
  # for sheep (create list)
  # Table 10.6
  # unsure of lamb sex, so go mean
  cspCsex2 <- cspParameters %>%
    filter(category %in% c("lambs")) %>%
    mutate(CsexA = if_else(grepl("lamb", animal)
                           , mean(c(2.1, 25))
                           , 0)
           , CsexB = if_else(grepl("lamb", animal)
                             , mean(c(0.35, .45))
                             , 0)
    )
  # merge back
  cspParameters <- merge(merge(cspParameters
                               , cspCsex[, c("animal", "category", "UK.Region", "Csex")]
                               , by = c("animal", "category", "UK.Region")
                               , all = T)
                         , cspCsex2[, c("animal", "category", "UK.Region", "CsexA", "CsexB")]
                         , by = c("animal", "category", "UK.Region")
                         , all = T)
  # tidy
  rm(cspCsex, cspCsex2)
  head(cspParameters)
  
  ### ensure both parameter and output dataset are in the same order
  cspOutputs <- cspOutputs[with(cspOutputs, order(animal, category, UK.Region)),]
  cspParameters <- cspParameters[with(cspParameters, order(animal, category, UK.Region)),]
  stopifnot(cspOutputs[, c(1:3)] == cspParameters[, c(1:3)])
  
  ## NEg equations ##
  cspOutputs$NEg <- 0
  
  for(i in 1:nrow(cspOutputs)){
    print(i)
    print(cspOutputs[i, ]$category)
    
    # for cattle
    # Eq 10.6 IPCC: NEg = 22.02 * (bw/ Csex * mw)^0.75 * dwg^1.097
    if(grepl("other_cattle|heifer|calf|pregnant|dairy", cspOutputs[i, ]$category)){  
      cspOutputs[i, ]$NEg <- 22.02 * (cspParameters[i, ]$Weight_Kg/ (cspParameters[i, ]$Csex * cspParameters[i, ]$mature_Kg))^0.75 * cspParameters[i, ]$wGain_Kg^1.097 
      # for sheep (only applies to immature sheep i.e. lambs)
      # Gibbs 2002 says sufficient to only include lamb weight gain
      # Eq 10.5 IPCC: NEg = (wg(yearly) * (a+0.5b(bw at weaning + bw when slaughtered)))/365
    } else if (grepl("lambs", cspOutputs[i, ]$animal)){
      # https://projectblue.blob.core.windows.net/media/Default/Beef%20&%20Lamb/GrowingAndFinishingLambsForBR3340_200415_WEB-1.pdf
      # states that an aim for lambs should be 18-21 kg when weaned
      # assumed same for goats
      cspOutputs[i, ]$NEg <- (cspParameters[i, ]$mature_Kg - 21) * (cspParameters[i, ]$CsexA + (0.5*cspParameters[i, ]$CsexA*(21 + cspParameters[i, ]$mature_Kg)))/365
    } else {
      # NEg for other animals not considered important in terms of emissions - therefore considered 0
      # i.e. horses/ poultry/ sheep
      cspOutputs[i, ]$NEg <- 0
    }
    print(cspOutputs[i, ]$NEg)
  }
  head(cspOutputs[sort(cspOutputs$NEg, decreasing = T), ])
  head(cspOutputs[sort(cspOutputs$NEg, decreasing = F), ])
  
  ###### Net energy lactation (NEl) #####
  ## ------------ Notes --------------  ##
  ## eq. 10.8/10 cows/sheep & goats (IPCC 2019)
  ## NEl unit: MJ day-1
  ## three additional parameters are required here: 
  ## amount of milk produced, kg of milk day-1
  ## fat content of milk, percent of total milk weight.
  ## the net energy required to produce 1 kg of milk.
  ## ------------ Notes --------------  ##
  
  ## milk - volume, fat, protein ##
  # for cows, 1 litre of milk is 1.033 kg - the litre values are found in John Nix,
  # but kg values need to be input for calculation
  # Friesians produced an average of 8,000 litres a year (Nix; p47)
  # which is equiv to 22 l a day, so 22.6411 kg day-1
  # (similar average to all British cows: https://ahdb.org.uk/dairy/uk-milk-yield)
  kgDayMilk <- 22.6411
  
  # fat and protein in Table 3.1 CFT (for average breed)
  fat <- (3.9+3.9+4.2+8+4.8+3.6+3.6+5+4.5+3.8+4.2+4.5+3.7)/14
  prot <- (3.8+3.3+3.7+3.5+3.5+3.6+3.2+3.2+4.2+4+3.3+3.5+4+3.2)/14
  
  # check
  stopifnot(cspOutputs[, c(1:3)] == cspParameters[, c(1:3)])
  
  # net energy production of milk
  # from IPCC (2019):
  # the energy required to produce 1 kg of milk, MJ kg-1. A default EVmilk value of 4.6 MJ/kg 
  # (sheep) (AFRC 1993; AFRC 1995) ... can be used which
  # corresponds to a milk fat content of 7 percent ... by weight for sheep ...,
  # respectively.
  cspParameters <- cspParameters %>%
    mutate(EVmilk = if_else(grepl("ewes", animal), 4.6, 0))
  head(cspParameters)
  
  ## NEl equations ##
  cspOutputs$NEl <- 0
  
  for(i in 1:nrow(cspOutputs)){
    # print(i)
    
    # for cattle
    # Eq 10.8 IPCC: NEl = milk production (kg day-1) * (1.47+0.4*fat content (%))
    if(grepl("lactating_dairy", cspOutputs[i, ]$category)){  
      cspOutputs[i, ]$NEl <- kgDayMilk * (1.47 + 0.4 * fat) 
      # for sheep, as the milk production is unknown the weight gain between birth and weaning 
      # important, so use 21 kg, as earlier -> 21-4 = 17 kg gained
      # Eq. 10.10 (IPCC, 2019): ((5*wg_bir_wean)/365) * EVmilk
    } else if (grepl("ewes", cspOutputs[i, ]$animal)){
      cspOutputs[i, ]$NEl <- ((5*17)/365) * cspParameters[i, ]$EVmilk
      # NEl for other animals not considered important in terms of emissions
    } else {
      cspOutputs[i, ]$NEl <- 0
    }
  }
  
  # tidy
  rm(kgDayMilk)
  str(cspOutputs)
  head(cspOutputs[order(cspOutputs$NEl, decreasing = TRUE), ])
  head(cspOutputs[order(cspOutputs$NEl, decreasing = F), ])
  
  ###### Net energy wool (NEw) #####
  ## ------------ Notes --------------  ##
  ## the net energy required to produce 1 kg of wool.
  ## eq. 10.12 sheep (IPCC 2019).
  ## goats can be considered but are not here due to goat wool not grown in uk
  ## NEw unit: MJ day-1
  ## ------------ Notes --------------  ##
  
  EVwool <- 24
  # annual wool production per sheep, kg yr-1
  # according to https://www.britishwool.org.uk/ksupload/userfiles/Shearing/Best_Management_Practices_27092017.pdf
  # between 2 - 5 kg 
  PRwool <- mean(c(2, 5))
  # NEw = ((EVwool * PRwool)/365)
  cspOutputs <- cspOutputs %>%
    mutate(NEw = if_else(grepl("sheep|ewes|ram", category) & !grepl("lamb|other_shee", animal)
                         , (EVwool * PRwool)/365
                         , 0))
  # tidy
  rm(EVwool, PRwool)
  head(cspOutputs[order(cspOutputs$NEw, decreasing = TRUE), ], 10)
  head(cspOutputs[order(cspOutputs$NEw, decreasing = F), ], 10)
  
  ###### Net energy pregnancy (NEp) #####
  ## ------------ Notes --------------  ##
  ## Eq 10.13 IPCC
  ## NEp unit: MJ day-1
  ## pregnancy - for sheep coefficient is between single birth and double birth
  ## Table 10.7 (IPCC 2019)
  ## ------------ Notes --------------  ##
  cspOutputs <- cspOutputs %>%
    mutate(NEp = if_else(grepl("pregnant", category) , NEm * 0.1
                         , if_else(grepl("ewes", animal), NEm * mean(c(0.077, 0.126))
                                   , 0)))
  
  # tidy
  rm(cspWeights)
  
  ###### 3b4 - calculate dry matter intake #####
  ## ---------- notes ---------- # 
  ## dry matter intake (DMI) is the amount of feed an animal consumes 
  ## on a moisture-free basis per day. DMI can be calculated
  ## with 3 additional bits of information:
  ##    estimated dietary net energy concentration of the feed or diet (NEmf)
  ##    body weight
  ##    fat-corrected milk (FCM)
  
  ## However, when using Eqs. 10.17 and 10.18 in IPCC (2019) for calves and
  ## growing cattle, respectively, the values of DMI per day were way too high
  ## (20% of BW, and 170% BW(!))
  ## Because of this, the approach for these two groups of animals were derived
  ## from the bodyweight alone. It is recommended ~2-2.5% of a cow’s body weight
  ## should be fed to them in the form of dry matter per day, so this will be
  ## used...
  
  ## it should also be noted that Brown et al. (2021) used a UK-specific way of deriving
  ## DMI. See section '5.3.2.1.1 Feed intake' in that publication.
  ## --------------------------- # 
  
  DMIinput <- cspParameters
  
  # FCM is determined based on volume and fat content of milk (EQ. 10.18b IPCC)
  kgDayMilk <- 22.6411 # check above for why this amount
  # fat and protein come from Table 3.1 in CFT - it is a mean value of all breeds
  fat <- mean(c(3.9,3.9,4.2,8,4.8,3.6,3.6,5,4.5,3.8,4.2,4.5,3.7))
  prot <- mean(c(3.8,3.3,3.7,3.5,3.5,3.6,3.2,3.2,4.2,4,3.3,3.5,4,3.2))
  kgDayMilkFat <- (fat/100) * kgDayMilk
  FCM <- (0.4324 * kgDayMilk) + (16.216 * kgDayMilkFat)
  # tidy
  rm(kgDayMilkFat, prot)
  
  ### ----- final used methodology ---- ###
  # calculate DMI - based on specific equations
  DMIinput$DMIkgDay <- ifelse(grepl("calf", DMIinput$category)|grepl("lamb", DMIinput$animal), DMIinput$Weight_Kg * 0.025
                              # heifers (EQUATION 10.18A)
                              , ifelse(DMIinput$category == "heifer", 3.184 + 0.01536 * DMIinput$Weight_Kg * 0.96
                                       # bulls/steers (EQUATION 10.18A)
                                       , ifelse(grepl("teer|ull", DMIinput$animal), 3.83 + 0.0143 * DMIinput$Weight_Kg * 0.96
                                                # lactating dairy cows (EQUATION 10.18B)
                                                , ifelse(grepl("lactating|pregn", DMIinput$category), 0.0185 * DMIinput$Weight_Kg + 0.305  * FCM
                                                         # EQUATION 10.18
                                                         , ifelse(grepl("other_cat|sheep|pig", DMIinput$category), DMIinput$Weight_Kg * 0.02
                                                                  # 2.5% bodyweight assumed for poultry DMI
                                                                  , 0.025 * DMIinput$Weight_Kg)))))
  # check the %s of BW
  pcWeight <- (DMIinput$DMIkgDay/DMIinput$Weight_Kg)*100
  # tidy
  rm(pcWeight)
  
  # calculate DMI - based on body weight percentage (for cows)
  # use these aims: https://ahdb.org.uk/knowledge-library/calculating-dry-matter-intakes-for-rotational-grazing-of-cattle
  DMIinput$DMIkgDayBW <- ifelse(grepl("calf", DMIinput$category), DMIinput$Weight_Kg * 0.0275
                                # lactating dairy cows (mean of 2 and 2.5%)
                                , ifelse(grepl("lactating|pregn", DMIinput$category), DMIinput$Weight_Kg * 0.0225
                                         # heifers (2.5%)
                                         , ifelse(DMIinput$category == "heifer", DMIinput$Weight_Kg * 0.025
                                                  # beef cattle is divided by weight
                                                  ## 300 - 600 kg (2.5% BW)
                                                  , ifelse(grepl("other_cat", DMIinput$category) & grepl("eef", DMIinput$animal) & DMIinput$Weight_Kg >= 300, DMIinput$Weight_Kg * 0.025
                                                           ## 100 - 200 kg (3% BW)
                                                           , ifelse(grepl("other_cat", DMIinput$category) & grepl("eef", DMIinput$animal) & DMIinput$Weight_Kg < 300, DMIinput$Weight_Kg * 0.03
                                                                    # other 'dry' cattle (1.5%)
                                                                    , ifelse(grepl("other_cat", DMIinput$category), DMIinput$Weight_Kg * 0.015
                                                                             , DMIinput$DMIkgDay))))))
  
  
  head(DMIinput)
  
  # combine the two tables and save
  parametersEnergy <- merge(DMIinput
                            , cspOutputs
                            , by = c("animal", "category", "UK.Region")
                            , all = T)
  
  ##### 3c - save tables ####
  fwrite(parametersEnergy, file.path(currPath, "animal_coefs_net_energy.csv")
         , row.names = F)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(currPath, "readme_animal_coefs_net_energy.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2023-05-03
Last update:",  format(Sys.Date()), "

Files:
  animal_coefs_net_energy.csv
  
animal_coefs_net_energy.csv
  This file contains categories and parameters and net energies for different animal groups, differentiated by region. 
    The columns:
          'animal' = the category an animal was assigned to, based on either CTS or AgCensus data
          'category' = the category each animal type was assigned to in this work
          'UK.Region' = different regions of the UK, with 'other' representing non-English regions
          'Weight_Kg' = the average weight, in kg, of a typical animal in this category and region combination
          'cfi' = Net energy maintenance coefficient of a typical animal in this category and region combination
          'Ca_[x]' = coefficient corresponding to animal’s feeding situation (x) for a typical animal in this category and region combination, with 'x' being 'housed' (housed all the time), 'flat_sml_past' (able to move across a small, flat pasture), and 'hill_lrg_past' (able to roam across a large, hilly, pasture)
          'mature_Kg' = the mature weight of an animal, in kg
          'wGain_Kg' = daily weight gain of a typical animal in this category and region combination, in kg
          'Csex' = sex coefficients. 0.8 for females, 1.2 for bulls (only required for cattle)
          'EVmilk' = the energy required to produce 1 kg of milk, in MJ kg-1
          'DMIkgDay' = the amount of dry matter intake (DMI) an animal requires in kg per day, based on specific equations
          'DMIkgDayBW' = the amount of dry matter intake (DMI) an animal requires in kg per day, based on the weight of a typical animal in this category and region combination
          'NEm' = Net energy required for maintenance, in MJ day-1 
          'NEa_[x]' = net energy for activity based on an animal’s feeding situation (x) for a typical animal in this category and region combination, with 'x' being 'housed' (housed all the time), 'flat_sml_past' (able to move across a small, flat pasture), and 'hill_lrg_past' (able to roam across a large, hilly, pasture)
          'NEg' = net energy required for growth, in MJ day-1 
          'NEl' = net energy required for lacatation, in MJ day-1 
          'NEw' = net energy required for wool production, in MJ day-1 
          'NEp' = net energy required for pregnancy, in MJ day-1
")
  
  sink(file = NULL)
  
} # end of 'part3AnimalAttributes'

#### 4 - part4Enteric ####
if(part4Enteric){
  # stop("part4")
  # set input location
  currPath <- file.path("./data_in", "animals")
  # set output location
  savePath <- file.path("./results", "animals")
  
  ## ------------ Notes --------------  ##
  ## to calculate enteric fermentation, the IPCC (2019) calculations were used. 
  ## This requires the gross energy requirements of an animal to be determined
  ## as well as their diet and situation (e.g. housed / free-roaming)
  
  ## for this run of the model, there is a big assumption: that all animals (with 
  ## the exception of lambs) spend half their time housed (where they are fed a 
  ## high quality (85% digestible energy) diet), and half of their time allowed
  ## to freely roam and graze a small field (where they eat a lower quality
  ## (65% digestible energy) diet). 
  ## ------------ ----- --------------  ##
  
  ##### 4a - gross energy amounts #####
  ## ------------ Notes --------------  ##
  ## gross energy (GE) needs to be determined for each animal to determine its
  ## overall emissions. Use the tables created in part 3.
  ## ------------ ----- --------------  ##
  
  # load in table with Net Energy amounts per animal category
  netEnergyTable <- fread(file.path(currPath, "animal_coefs_net_energy.csv")) %>%
    as.data.frame()
  head(netEnergyTable)
  
  ###### 4a1 - Ratios of net energy available in a diet ######
  ## ------------ Notes --------------  ##
  ## Ratio of net energy available in a diet for maintenance (REM) and growth (REG)
  ## to digestible energy consumed 
  ## Eq 10.14/10.15 IPCC (REM/REG)
  ## an additional parameter is required here: 
  ## digestible energy (DE) of food consumed
  
  ## using Table 3.7 in CFT, digestible energy (DE) can range from about
  ## 0.5 - 0.9 of GE based on food choice
  ## REM = 1.123 - (4.092*10^-3*DE) + (1.126*10^-5*(DE)^2) - (25.4/DE)
  ## REG = 1.164 - (5.16*10^-3*DE) + (1.308*10^-5*(DE)^2) - (37.4/DE)
  ## important note: IPCC (2019) says to use fractions, but Gibbs et al. (2002) says to use %
  ## the latter was used as it gives non-negative results, which makes sense
  
  ## for this work, we used the assumption of either 65% DE in a diet (when outside)
  ## or 85% when inside (i.e. housed)
  ## ------------ ----- --------------  ##
  
  # loop through possible values of DE, and calculate resulting REM and REG
  for(de in seq(65, 85, 1)){
    print(de)
    REM <- 1.123 - (4.092*10^-3*de) + (1.126*10^-5*(de^2)) - (25.4/de)
    print(REM)
    REG <- 1.164 - (5.16*10^-3*de) + (1.308*10^-5*(de^2)) - (37.4/de)
    print(REG)
    # calculate GE based on Equation 10.16 in IPCC (2019)
    # activity can be split into 3 categories
    GEremHoused <- (netEnergyTable$NEm + netEnergyTable$NEa_housed + netEnergyTable$NEl + netEnergyTable$NEp)/REM
    GEremSml <- (netEnergyTable$NEm + netEnergyTable$NEa_flat_sml + netEnergyTable$NEl + netEnergyTable$NEp)/REM
    GEremLrg <- (netEnergyTable$NEm + netEnergyTable$NEa_hill_lrg + netEnergyTable$NEl + netEnergyTable$NEp)/REM
    GEreg <- (netEnergyTable$NEg + netEnergyTable$NEw)/REG
    GEremregHoused <- GEremHoused+GEreg
    GEremregSml <- GEremSml+GEreg
    GEremregLrg <- GEremLrg+GEreg
    MJdayGEhoused <- GEremregHoused/(de/100)
    MJdayGEsml <- GEremregSml/(de/100)
    MJdayGElrg <- GEremregLrg/(de/100)
    
    if(de == 65){
      DEtable <- as.data.frame(cbind(animal = netEnergyTable$animal
                                     , category = netEnergyTable$category
                                     , UK.Region = netEnergyTable$UK.Region
                                     , DE = de
                                     , MJdayGEhoused = MJdayGEhoused
                                     , MJdayGEsml = MJdayGEsml
                                     , MJdayGElrg = MJdayGElrg))
    } else {
      DEtable <- rbind(DEtable
                       , as.data.frame(cbind(animal = netEnergyTable$animal
                                             , category = netEnergyTable$category
                                             , UK.Region = netEnergyTable$UK.Region
                                             , DE = de
                                             , MJdayGEhoused = MJdayGEhoused
                                             , MJdayGEsml = MJdayGEsml
                                             , MJdayGElrg = MJdayGElrg)))
    }
  }
  head(DEtable)
  
  ###### 4a2 - gross energy requirements based on REG / REM ######
  # pivot with regards to DE%
  GEtableOut <- DEtable %>%
    group_by(DE) %>%
    mutate(row = row_number()) %>%
    pivot_wider(names_from = "DE"
                , values_from = c("MJdayGEhoused", "MJdayGEsml", "MJdayGElrg")) %>%
    dplyr::select(-row) %>%
    mutate(across(c(4:ncol(.)), as.numeric))
  head(GEtableOut)
  
  # save
  fwrite(GEtableOut, file.path(savePath, "mjDayGE.csv"), row.names = F)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(savePath, "readme_mjDayGE.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-05-19
  Last update:",  format(Sys.Date()), "
  
  Description of files in directory:
  The file contains the gross energy (GE) requirements of an animal per day, in MJ. This is based on the net energy requirements of the animals, which originally come from type, properties and growth stage of the animals
  
  Files:
    mjDayGE.csv
  
  Columns of mjDayGE.csv:
        'animal' = the category an animal was assigned to, based on either CTS or AgCensus data
        'category' = the category each animal type was assigned to in this work
        'UK.Region' = different regions of the UK, with 'other' representing non-English regions
        'MJdayGE[freedom]_[DE]' = the GE required for an animal per day, based on its freedom and DE
          where:
            'freedom' = either 'housed' (indoors), 'sml' (is able to roam a small field), or 'lrg' (being able to roam a large, hilly field) 
            'DE' = 65 - 85
    
  Data: GE required per animal category per day per animal
  units: MJ per day
  Spatial extent: UK
  Native projection: 27700")
  
  sink(file = NULL)
  # tidy
  rm(GEremHoused, GEreg, GEremLrg, GEremregHoused, GEremregSml, GEremSml
     , GEremregLrg, de)
  
  # check to ensure that it is about right i.e. 1.5 percent and 3.0 percent of the animal’s weight.
  # to check, (Gibb, 2002) the estimate can be converted in daily intake in kilograms by dividing by 18.45 MJ/kg.
  # check middling column
  GEcheck <- as.data.frame(cbind(animal = netEnergyTable$animal
                                 , GE65 = GEtableOut$MJdayGEsml_65 / 18.45
                                 , GE85 = GEtableOut$MJdayGEsml_85 / 18.45
                                 , weightFrom = netEnergyTable$Weight_Kg * 0.015
                                 , weightTo = netEnergyTable$Weight_Kg * 0.03
                                 , actualWeight = netEnergyTable$Weight_Kg)) %>%
    mutate_at(c("GE65", "GE85",  "weightFrom", "weightTo", "actualWeight"), as.numeric) %>%
    rowwise() %>%
    mutate(fitCriteria65 = if_else(GE65 >= weightFrom & GE65 <= weightTo, "Yes"
                                   , if_else(GE65 < weightFrom, "NO - too low"
                                             , if_else(GE65 > weightTo, "NO - too high", "NO")))) %>%
    mutate(fitCriteria85 = if_else(GE85 >= weightFrom & GE85 <= weightTo, "Yes"
                                   , if_else(GE85 < weightFrom, "NO - too low"
                                             , if_else(GE85 > weightTo, "NO - too high", "NO"))))
  head(GEcheck)
  
  # Using values from GHGi (2019; appendix A3.3.1) (for the UK), check dairy cows have enough GE
  GEcsv <- GEtableOut
  head(GEcsv)
  # keep rows if they contain dairy
  dairyGE <- GEcsv %>%
    filter(grepl("dairy", category)| grepl("dairy", animal) & category != "calf") %>%
    # group by animal and category, to get average regardless of region
    group_by(animal, category) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    rowwise %>% 
    do(as.data.frame(.) %>% { 
      # select just housed 85 and small field 65
      subs <- dplyr::select(., MJdayGEhoused_85: MJdayGEsml_75)
      mutate(., Min = subs %>% min,
             Max = subs %>% max) 
    } ) %>%
    ungroup %>%
    rowwise() %>%
    mutate(fitCriteria = if_else(Min >= 263.80 & Max <= 382.51, "Yes"
                                 , if_else(Min < 263.80, "NO - too low"
                                           , if_else(Max > 382.51, "NO - too high", "NO"))))
  head(dairyGE)
  
  ## ------------ Notes --------------  ##
  ## for this to be seen as accurate, range of daily GE should include from 263.80
  ## to 382.51 GE (MJ/day)
  ## ------------ ----- --------------  ##
  
  UpperLim <- 382.51
  lowerLim <- 263.80
  
  # see which are too high low, specifically
  dairyGEcheck2 <- GEcsv %>%
    filter(grepl("dairy", category)| grepl("dairy", animal) & category != "calf") %>%
    # group by animal and category, to get average regardless of region
    group_by(animal, category) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    mutate(across(where(is.numeric), ~ ifelse(. > lowerLim & . < UpperLim, 1, 0)))
  head(dairyGEcheck2)
  # check the final amount
  colSums(dairyGEcheck2[, 3:ncol(dairyGEcheck2)])
  table(as.matrix(dairyGEcheck2[, 3:ncol(dairyGEcheck2)]))
  
  # some use more/less GE than suggested above (411/ 935: 44%)
  # it appears that the general trend is that young animals e.g. cows < 6 months 
  # need a higher energy (when compared to body weight) use when growing, but this could be justified
  
  ## ------------ Notes --------------  ##
  ## conversion of GE into emissions depends on the additional parameter of
  ## Methane conversion factor (specific for each livestock category) [%] (Ym)
  ## in simple terms, Ym is estimate of CH4 loss per unit of feed
  ## Ym changes with breed, region, stage of life, activity, etc., and it plays
  ## a very big part in determining overall emissions
  ## it is good practice to determine nation-specific Ym values, based highly on milk production
  ## as that cannot be determined for the whole of England, FAO's 
  ## formula is used here (see Eq. 3.2.30 in CFT):
  ## Ym = 9.75 - 0.05 * DE%
  ## Note: IPCC's default value is 6.5
  
  ## to determine final emissions:
  ## Emissions (kg/yr) = [GE Intake (MJ/day) * Ym * (365 days/yr)] / [55.65 MJ/kg of methane]
  ## eq. 10.21 sheep (IPCC 2019)
  ## ------------ ----- --------------  ##
  
  for(de in c(65, 70, 75, 80, 85)){
    print(de)
    # extract just those columns from GEtableOut
    GEDE <- GEtableOut %>%
      dplyr::select(contains(paste0("_",de)))
    # derive appropriate Ym
    Ym <- 9.75 - 0.05 * de
    print(Ym)
    Ym100 <- Ym/100
    # calculate GE * Ym * 365
    kgCH4GE <- (GEDE * Ym100 * 365)/55.65
    # convert to CO2e
    kgCO2GE <- kgCH4GE * 25
    
    if(de == 65){
      DEtable2 <- as.data.frame(cbind(animal = netEnergyTable$animal
                                      , category = netEnergyTable$category
                                      , UK.Region = netEnergyTable$UK.Region
                                      , kgCH4GE = kgCH4GE))
      # save CO2 version
      DEtableCO2 <- as.data.frame(cbind(animal = netEnergyTable$animal
                                        , category = netEnergyTable$category
                                        , UK.Region = netEnergyTable$UK.Region
                                        , kgCO2GE = kgCO2GE))
    } else {
      DEtable2 <- cbind(DEtable2
                        , kgCH4GE = kgCH4GE)
      # save CO2 version
      DEtableCO2 <- cbind(DEtableCO2
                          , kgCO2GE = kgCO2GE)
    }
  }
  head(DEtable2)
  head(DEtableCO2)
  
  # run checks to ensure daily enteric CH4 is consistent with literature
  # keep rows if they contain dairy
  dairyCH4 <- DEtable2 %>%
    filter(grepl("dairy", category)| grepl("dairy", animal) & category != "calf") %>%
    # group by animal and category, to get average regardless of region
    group_by(animal, category) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    rowwise %>% 
    do(as.data.frame(.) %>% { 
      subs <- dplyr::select(., kgCH4GE.MJdayGEhoused_65:kgCH4GE.MJdayGElrg_85)
      mutate(., Min = subs %>% min,
             Max = subs %>% max) 
    } ) %>%
    ungroup %>%
    rowwise() %>%
    mutate(fitCriteria = if_else(Min <= 125.44 & Max >= 125.44, "Yes"
                                 , "NO"))
  head(dairyCH4)
  # for dairy cows, 125.44 kg CH4/head/year (GHGi, 2019) needs to be covered, which it is
  # for dairy heifers, 53.20 kg CH4/head/year (GHGi, 2019) needs to be covered, which it is
  
  # as we are going to use the scenario of 50% time housed with 80% DE and 50% grazing with 85% DE
  # see whether the half and half is ok
  dairyCH4v2 <- dairyCH4 %>%
    dplyr::select(animal, category, kgCH4GE.MJdayGEsml_75, kgCH4GE.MJdayGEhoused_80) %>%
    mutate(kgCH4GE.50pcScen = (kgCH4GE.MJdayGEsml_75/2) + (kgCH4GE.MJdayGEhoused_80/2)) %>%
    rowwise() %>%
    mutate(litDifference = kgCH4GE.50pcScen - 125.44)
  head(dairyCH4v2)
  
  # stop("ouch")
  
  # use this scenario to give a per-animal enteric fermentation emission
  # although lambs always housed
  CH4EntFerm <- DEtable2 %>%
    dplyr::select(animal, category, UK.Region, kgCH4GE.MJdayGEsml_80, kgCH4GE.MJdayGEhoused_85) %>%
    rowwise() %>%
    mutate(kgCH4yr.50pcScen = ifelse(grepl("lamb", category), kgCH4GE.MJdayGEhoused_85
                                     , (kgCH4GE.MJdayGEsml_80/2) + (kgCH4GE.MJdayGEhoused_85/2))) %>%
    # convert to CO2e
    mutate(kgCO2.50pcScen = kgCH4yr.50pcScen * 25)
  head(CH4EntFerm)
  
  # save final per-animal amount of CH4 and CO2e enteric fermentation
  fwrite(CH4EntFerm
         , file.path(savePath, "entericFermCH4CO2e.csv")
         , row.names = F)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(savePath, "readme_entericFermCH4CO2e.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-02-13
  Last update:",  format(Sys.Date()), "

  Files:
    entericFermCH4CO2e.csv
  
  Columns of entericFermCH4CO2e.csv:
      'animal' = the category an animal was assigned to, based on either CTS or AgCensus data
      'category' = the category each animal type was assigned to in this work
      'UK.Region' = different regions of the UK, with 'other' representing non-English regions
      'kgCH4GE.MJdayGE[freedom]_[DE]' = the enteric fermentation emissions, in kg CH4, based on the GE required for an animal per day, based on its freedom and DE
        where:
          'freedom' = either 'housed' (indoors), 'sml' (is able to roam a small field), or 'lrg' (being able to roam a large, hilly field) 
          'DE' = 65 - 85
      'kgCO2.50pcScen' = the final amount of per-animal enteric fermentation emissions, in kg CO2e, based on a scenario in which the animals were housed 50% of the time and fed a 85% digestible energy diet, and could roam and graze 50% of the time, with a 80% DE diet
    
  Data: per-animal enteric fermentation emissions, in kg CO2e or kg cH4
  units: kg CO2e or kg cH4 per year")
  sink(file = NULL)
  
} # end of 'part4Enteric'

#### 5 - part5ManureManagement ####
if(part5ManureManagement){
  
  # stop("part 5")
  
  # set input location
  currPath <- file.path("./data_in", "animals")
  # set output location
  savePath <- file.path("./results", "animals")
  
  # load in table with Net Energy amounts per animal category
  netEnergyTable <- fread(file.path(currPath, "animal_coefs_net_energy.csv")) %>%
    as.data.frame()
  head(netEnergyTable)
  
  ##### manure characteristics #####
  ## ---------- notes ---------- # 
  ## The IPCC Tier 2 method relies on two primary types of inputs that affect the 
  ## calculation of methane emission factors from manure:
  ##  manure characteristics: volatile solids (VS) produced in the manure,
  ##                          and the amount of methane able to be produced from that manure (B0).
  ##  Animal waste management system characteristics (AWMS; the proportion that that management type composed all manure):
  ##                          Includes the types of systems used to manage manure 
  ##                          and a system-specific methane conversion factor (MCF)
  ## The latter reflects the portion of B0 that is achieved
  ## --------------------------- # 
  
  ###### 5a - volatile solids (VS) ######
  ## ---------- notes ---------- # 
  ## Default VS value from IPCC (2019) is 8.4 (1000 kg animal mass-1) day-1 
  
  ## If average daily VS excretion rates are not available, country-specific VS
  ## excretion rates can be estimated from feed intake levels (IPCC, 2019). These
  ## were calculated in the ghg_excretions.
  ## --------------------------- # 
  
  # load in the GE required per day, from which VS can be calculated
  dailyGE <- fread(file.path(savePath, "mjDayGE.csv")) %>%
    as.data.frame()
  head(dailyGE)
  
  # from using EQUATION 10.24 in IPCC (2019), VS can be calculated per animal type
  # VS = [GE * (1 - DE100) + (UE * GE)]*[((1 - ASH)/18.45)]	
  # create VS table with all the necessary columns
  VStable <- dailyGE[, 1:3] %>%
    # assign ASH value (0. 06 for sows (Dämmgen et al. 2011); 0.08 for cattle)
    # with the exception of pigs, assign all others 0.08
    mutate(ash = if_else(category == "pigs", 0.06, 0.08)) %>%
    # calculate right side of equation
    mutate(rightSide = (1-ash)/18.45)
  # add in each column from dailyGE table using GE and DE from the results
  head(VStable)
  
  # add in each column from dailyGE table using GE and DE from the results
  for(i in 4:ncol(dailyGE)){
    # get name
    nm <- names(dailyGE)[i]
    # convert to kg VS per day
    nm <- gsub("MJdayGE", "kgVSday", nm)
    # get de
    de <- as.numeric(gsub(".*_(.+)$", "\\1", nm))
    # get GE results
    GEs <- dailyGE[, i]
    # multiply GE by UE
    UEGE <- GEs * 0.04
    # calc de in terms of proportion
    deProp <- 1 - (de/100)
    # calculate left side
    leftSide <- GEs * deProp + UEGE
    # times left and right
    both <- leftSide * VStable$rightSide
    print(head(leftSide))
    # add to VS table
    VStable <- VStable %>%
      mutate(both = both) %>%
      rename(!!nm := "both")
    
  }
  head(VStable)
  
  # tidy
  rm(de, deProp, both, leftSide, UEGE, GEs, nm)
  
  # get means for all of the same animal, category
  summariseVS <- VStable %>%
    group_by(animal, category) %>%
    dplyr::summarise(across(starts_with("kgVS"), list(mean = ~ mean(., na.rm = TRUE)))) %>%
    # dplyr::select(animal, category, contains(c("housed_80", "sml_80"))) %>%
    filter(grepl("lactat", category))
  head(summariseVS)
  
  # from Appuhamy et al. (2018), the mean VS (of lactating dairy cows) was 5.73 kg/d
  
  # save tables
  # average VS per day
  fwrite(VStable, file.path(savePath, "dailyVS.csv"), row.names = F)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(savePath, "readme_dailyVS.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-03-17
  Last update:",  format(Sys.Date()), "
  
  Description of files in directory:
  The files containing 'dailyVS.csv' show the calculated amount of volatile solids (VS) produced per day per animal, in kilograms. The values are based on assumptions to do with the animals' situations i.e. the quality of their diet, and their freedom of moved (e.g. housed, or allowed to roam)  

  Columns of entericFermCH4CO2e.csv:
      'animal' = the category an animal was assigned to, based on either CTS or AgCensus data
      'category' = the category each animal type was assigned to in this work
      'UK.Region' = different regions of the UK, with 'other' representing non-English regions
      'kgVSday[freedom]_[DE]' = VS produced, in kg, based on the GE required for an animal per day, based on its freedom and DE
        where:
          'freedom' = either 'housed' (indoors), 'sml' (is able to roam a small field), or 'lrg' (being able to roam a large, hilly field) 
          'DE' = 65 - 85
          
  Data: per-animal VS production
  units: kg VS per day per animal")
  
  sink(file = NULL)
  
  ###### 5b - per-animal volatile solids (VS) ######
  ## ------------ Notes --------------  ##
  ## this final output will be based upon the scenario of DE written in the 
  ## '0 - select animal scenario' section above
  ## ------------ ----- --------------  ##
  
  # set names to keep
  insideName <- paste0("kgVSdayhoused_", DEinside)
  outsideName <- paste0("kgVSdaysml_", DEoutside)
  
  VStableScen <- VStable %>%
    dplyr::select(animal, category, UK.Region
                  , (all_of(insideName))
                  , (all_of(outsideName))) %>%
    rename("kgVSdayIn" = insideName
           , "kgVSdayOut" = outsideName) %>%
    # get the 6-month value of each 
    # per 6-month value
    mutate(kgVS6mnIn = kgVSdayIn * (365/2)
           , kgVS6mnOut = kgVSdayOut * (365/2))
  
  head(VStableScen)
  
  # reduce
  VStableScen <- VStableScen %>%
    dplyr::select(-c(kgVSdayIn, kgVSdayOut))
  head(VStableScen)
  
  # tidy
  rm(GEDE)
  
  ###### 5c - B0 (maximum methane producing capacity) ######
  ## ---------- notes ---------- # 
  ## Default B0 values are used from IPCC (2019) Table 10.16
  ## using 'Western Europe' values
  ## Units for B0 are M^3 CH4 KG^-1 VS
  ## --------------------------- # 
  
  # add B0 to animal attributes table
  VStableScen <- VStableScen %>%
    mutate(B0 = if_else(grepl("lactating_dairy|pregnant|heifer", category) , 0.24
                        , if_else(grepl("other_cattle|calf", category), 0.18
                                  , if_else(category == "pigs", 0.45
                                            , if_else(category == "poultry", 0.39
                                                      , if_else(grepl("ewes|rams|lambs|other_sheep", category) , 0.19 
                                                                , if_else(category == "poultry", 0.30
                                                                          , 0)))))))
  
  ###### 5d - AWMS ######
  ## ---------- notes ---------- # 
  ## AWMS(T,S,k) = fraction of livestock category T's manure handled using animal 
  ## waste management system S in climate region k.
  ## k, the climate region, is always considered 'Cool Temperate Moist' as the data 
  ## are for the UK
  ## T and S change
  
  ## values for the AWMS were extracted from Table 10.14 in IPCC (2019)
  ## --------------------------- # 
  
  # determine the type of management that is used for different animals manure
  # this is based on the table of Brown et al. (2021), specifically Table A 3.3.5
  awmsTable <- as.data.frame(cbind(category2 = c("lactating_dairy|pregnant|heifer", "other_cattle|calf", "pigs"
                                                 , "ewes", "rams", "lambs|other_sheep", "horse", "poultry")
                                   , liquid = c(60, 18, 36, rep(0, 5))
                                   , daily_spread = c(8, 12, 14, rep(0, 4), 39)
                                   , FYM_manure = c(9, 21, 36, 8, 1, 1, 30, 50)
                                   , pasture_range_paddock = c(22, 48, 10, 92, 99, 99, 70, 3)
                                   , anaerobic_digest_biogas = c(2, 1, 4, 0, 0, 0, 0, 8))) %>%
    mutate(across(.cols = c(2:6), .fns = as.numeric)) %>%
    # get sum of all parts
    mutate(total = rowSums(across(where(is.numeric)))) %>%
    # convert all columns to fractions of the whole
    mutate(across(.cols = c(2:6), ~ ./total)) %>%
    # remove total
    dplyr::select(-total)
  head(awmsTable)
  
  # manipulate the awms table to create more rows where different categories encompassed
  # different animal groups
  awms2 <- awmsTable %>%
    mutate(divides = str_count(category2, "\\|")) %>%
    rowwise() %>%
    slice(rep(1:n(), each = divides+1)) %>%
    group_by(category2) %>%
    mutate(Code_n = make.unique(as.character(category2))) %>%
    # get unique code as an index
    # filter(grepl("\\.", Code_n)) %>%
    mutate(x2 = gsub(".*\\.", "\\1", Code_n)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(x = list(str_split(Code_n, "\\|", simplify = TRUE)[, (as.numeric(x2) + 1)])) %>%
    mutate(x = gsub("^(.+)\\..*", "\\1", x)) %>%
    # correct the category
    mutate(category2 = if_else(is.na(x), str_remove(category2, "\\|.*")
                               , x)) %>%
    dplyr::select(-c(divides, Code_n, x2, x)) 
  rm(awmsTable)
  print(awms2)
  
  ###### 5e - methane conversion factors (MCFs) ######
  # Table A 3.3.3 in Brown et al. (2021)
  # except 'Anaerobic digestion', which was taken from IPCC (2019)
  # NOTE: these values are specific to the UK
  mcfTable <- rbind(liquid = 17
                    , daily_spread = 0.1
                    , FYM_bovine_porcine = 17
                    , FYM_ovine = 2
                    , pasture_range_paddock = 1
                    , poultry_manure = 1.5
                    , anaerobic_digest_biogas = mean(c(1, 1.41, 3.55, 9.59, 10.85, 12.14)))
  # divide by 100
  mcfTable <- mcfTable / 100
  head(mcfTable)
  
  # convert the specific manure ones to the specific manure groups
  VStableScenMM <- VStableScen %>%
    mutate(mcf_FYM_manure = if_else(grepl("lactating_dairy|pregnant|heifer|other_cattle|calf|pigs", category) , mcfTable[3,1]
                                    , if_else(category == "poultry", mcfTable[6,1]
                                              , if_else(grepl("ewes|rams|lambs|other_sheep|horse", category) , mcfTable[4,1] 
                                                        , 0)))
           # add the others (i.e. non-manure ones)
           , mcf_liquid = mcfTable[1,1]
           , mcf_daily_spread = mcfTable[2,1]
           , mcf_pasture = mcfTable[5,1]
           , mcf_digest = mcfTable[7,1])
  
  head(VStableScenMM)
  
  # ensure both AWMS and mcf tables are in the same order in terms of manure management type
  ## merge using category 
  VStableScenMM <- VStableScenMM %>%
    relocate(mcf_liquid, mcf_daily_spread, mcf_FYM_manure, mcf_pasture, mcf_digest) %>%
    merge(., awms2
          , by.x = "category"
          , by.y = "category2"
          , all = F)
  head(VStableScenMM)
  
  ## extract both mcfs and awms
  ### mcfs
  aaMcf <- VStableScenMM %>%
    dplyr::select(contains("mcf_"))
  ### awms
  aaAWMS <- VStableScenMM %>%
    dplyr::select(liquid:anaerobic_digest_biogas)
  
  ## multiply together
  ## ---------- notes ---------- # 
  ## This next part will give the right side of equation 10.23 in IPCC (2019)
  ## It assumes that the different management systems correspond to Table A 3.3.5 in Brown et al. (2021)
  ## The left side will be determined by using previously-derived VS, which itself depends on DE%
  ## --------------------------- # 
  summedMcfAwms <- as.data.frame(aaMcf * aaAWMS) %>%
    mutate(total = rowSums(.[])) %>%
    # add categories back in
    bind_cols(VStableScenMM %>% dplyr::select(category, animal, UK.Region, B0), .) %>%
    # multiply the total summed by B0
    mutate(EF_right = total * B0)
  head(summedMcfAwms)
  
  ## ---------- notes ---------- # 
  ## now the VS table is used to multiply the values obtained above
  ## first daily CH4 from manure management will be obtained, with the amount  
  ## then be multiplied by 365 to get yearly values
  ## --------------------------- # 
  # ensure that rows completely match first
  ## reorganise so it matches previous animal categories
  x <- match(paste(VStable$category, VStable$animal, VStable$UK.Region)
             , paste(summedMcfAwms$category, summedMcfAwms$animal, summedMcfAwms$UK.Region))
  # reorder by the index
  VStable <- VStable[order(x), ]
  # check
  stopifnot(summedMcfAwms$animal == VStable$animal)
  
  # select just the VS values from the inside proportion, otherwise manure is left on-field
  VSonly <- VStable %>%
    dplyr::select((all_of(insideName))) %>%
    rename("kgVSdayIn" = insideName) 
  head(VSonly)
  
  # multiply the 'right side' calculation by each VS
  dailyMMCH4 <- as.data.frame(summedMcfAwms$EF_right * VSonly)
  # rename
  names(dailyMMCH4) <- gsub("kgVSday", "MMkgCHday", names(dailyMMCH4))
  head(dailyMMCH4)
  
  # get 6-monthly data
  sixMonthMMCH4 <- as.data.frame(dailyMMCH4 * (365 / 2)) %>%
    # add animal info back in
    bind_cols(VStable %>% select(animal, category, UK.Region), .)
  # rename
  names(sixMonthMMCH4) <- gsub("MMkgCHday", "MMkgCH6Mn", names(sixMonthMMCH4))
  head(sixMonthMMCH4)
  
  # get means for all of the same animal, category
  summariseCH <- sixMonthMMCH4 %>%
    filter(grepl("dairy", category)| grepl("dairy", animal) & category != "calf") %>%
    group_by(animal, category) %>%
    dplyr::summarise(across(starts_with("MMkgCH6"), list(mean = ~ mean(., na.rm = TRUE))))
  head(summariseCH)
  
  # get range
  CHcheckMM <- summariseCH %>% 
    rowwise %>% 
    do(as.data.frame(.) %>% { 
      subs <- select(., MMkgCH6MnIn_mean)
      mutate(., Min = subs %>% min,
             Max = subs %>% max) 
    } ) %>%
    ungroup
  head(CHcheckMM)
  
  ###### 5f - convert CH4 to CO2e ######
  # convert to CO2e
  sixMonthMMCO2e <- sixMonthMMCH4 %>%
    mutate(across(where(is.numeric), ~ . * 25))
  # rename
  names(sixMonthMMCO2e) <- gsub("MMkgCH", "MMkgCO2e", names(sixMonthMMCO2e))
  head(sixMonthMMCO2e)
  
  # save tables
  # CH4
  fwrite(sixMonthMMCH4, file.path(savePath, "sixMonthMMCH4.csv"), row.names = F)
  # CO2
  fwrite(sixMonthMMCO2e, file.path(savePath, "sixMonthMMCO2e.csv"), row.names = F)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(savePath, "readme_sixMonthMM.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-07-13
  Last update:",  format(Sys.Date()), "
  
  The files show the calculated amount of emissions from manure management per animals, per six months, to do with the animals' situations i.e. the quality of their diet, and their freedom of moved (e.g. housed, or allowed to roam)  
  
  Files:
    sixMonthMMCH4.csv
    sixMonthMMCO2e.csv
  
  Columns of 'sixMonthMM' files:
      'animal' = the category an animal was assigned to, based on either CTS or AgCensus data
      'category' = the category each animal type was assigned to in this work
      'UK.Region' = different regions of the UK, with 'other' representing non-English regions
      'MMkg[ghg]6MnIn' = emissions produced of [ghg] from manure management, in kg, per-animal for a six-month stay indoors, based on diet
          where [ghg] = either CO2e (carbon dioxide equivalent) or CH4 (methane)
          
  Data: per-animal  production
  units: kg VS per day per animal")
  sink(file = NULL)
  
} # end of 'part5ManureManagement'

#### 6 - part6excretions ####

if(part6excretions){
  
  stop("part 6")
  
  # set input location
  currPath <- file.path("./data_in", "animals")
  # set output location
  savePath <- file.path("./results", "animals")
  
  # import GE inputs per animal
  GEinput <- fread(file.path(savePath, "mjDayGE.csv"))%>%
    as.data.frame()
  head(GEinput)
  
  # import dry matter required per animal
  DMIintake <- fread(file.path(currPath, "animal_coefs_net_energy.csv"))
  head(DMIintake)
  
  ##### 6a - N balance #####
  ## ---------- notes ---------- # 
  ## N balance (specifically N export and intake) is related to diet - the amount
  ## and quality of food. This balance is therefore linked to grazing of cattle.
  ## The excretion rate is the difference between the feed nitrogen intake 
  ## and the nitrogen retention (in body or milk).
  
  ## further emissions (but indirect ones) come from emissions from leaching
  ## and runoff and emissions from N volatilization
  ## --------------------------- # 
  
  ###### 6a1 - calculate N intake ######
  ## ---------- notes ---------- # 
  ## N intake requires percent crude protein of dry matter
  ## for Western Europe generally it is 16.1 (Table 10A.1)
  ## however, this changes based on activity
  ## here, using 'Feed' Table in Nix (2020; p304), weighted averages were calculated for three scens:
  ## always inside, half in/ half out (small), always out (large)
  ## --------------------------- # 
  
  # daily N intake is different for cows and sheep compared to others
  # Eqs. 10.32 and 10.32A IPCC (2019)
  
  ## food for scens are: 
  ## inside: quarter cereals/ compound foods/ pulses/ grass
  cereals <- mean(c(128, 120, 94, 117, 82, 112, 107, 109))
  compund <- mean(c(518, 378, 371, 371, 168))
  pulses <- mean(c(397, 243, 266, 302))
  grass <- mean(c(233, 233, 127, 134, 134, 177, 177, 79, 193
                  , 157, 44, 182, 134, 177, 197, 91, 197))
  insideCP <- mean(c(cereals, compund, pulses, grass))/10
  # half in/ half out (small): third grass/ third inside ones/ third other
  other <- mean(c(192, 237, 207, 237, 105))
  halfInOut <- mean(c(cereals, pulses, grass, other))/10
  # always out (large): 90% grass, with supplementary other
  alwaysOut <- mean(c(grass, grass, grass, grass
                      , grass, grass, grass, grass
                      , grass, other))/10
  # tidy
  rm(cereals, compund, grass, other, pulses)
  
  # pivot GE table for ease of N intake
  GEpivot <- GEinput %>%
    pivot_longer(!c("animal", "category", "UK.Region")) %>%
    rename(MJday = value)
  head(GEinput)
  head(GEpivot) 
  unique(GEpivot$category)
  
  ###### N intake for bovine and ovine ######
  # assign an N intake category - separate cow and sheep from others
  GEpivot$NinCat <- ifelse(grepl("other_cattle|calf|heifer|pregnant|lactating_dairy|sheep|ewe|ram|lambs"
                                 , GEpivot$category), "cowsheep", "other")
  table(GEpivot$NinCat)
  
  GEpivot %>%
    filter(category == "other_sheep") %>%
    distinct(NinCat)
  head(GEpivot)
  
  # using categories, use correct equation
  # use EQUATION 10.32
  # units: kg N animal-1 day-1
  
  # EQUATION 10.32 (IPCC, 2019)
  csgNintakeDay <- function(GE, CP){
    (GE/18.45) * ((CP/100)/6.25)
  }
  
  # assigning different intakes of N based on feeding regimes
  # for the time being, assume all CP is the Western Europe average of 16.1%
  GEpivot$kgNintakeDay <- ifelse(GEpivot$NinCat == "cowsheep" & grepl("house", GEpivot$name)
                                 # , csgNintakeDay(GEpivot$MJday, insideCP)
                                 , csgNintakeDay(GEpivot$MJday, 16.1)
                                 , ifelse(GEpivot$NinCat == "cowsheep" & grepl("sml", GEpivot$name)
                                          # , csgNintakeDay(GEpivot$MJday, halfInOut)
                                          , csgNintakeDay(GEpivot$MJday, 16.1)
                                          , ifelse(GEpivot$NinCat == "cowsheep" & grepl("lrg", GEpivot$name)
                                                   # , csgNintakeDay(GEpivot$MJday, alwaysOut)
                                                   , csgNintakeDay(GEpivot$MJday, 16.1)
                                                   , NA)))
  head(GEpivot)
  
  ###### for porcine, poultry, equines ######
  # use EQUATION 10.32A
  # units: kg N animal-1 day-1
  
  # EQUATION 10.32A (IPCC, 2019)
  otherNintakeDay <- function(DMI, CP){
    (DMI) * ((CP/100)/6.25)
  }
  
  # cycle through columns assigning different CP values based on feeding regimes
  # but DMI was not differentiated indoor/ out etc.
  # so repeat with the different categories
  # combine DMI with long df
  GEpivot <- merge(GEpivot, DMIintake
                   , by = names(DMIintake)[1:3]
                   , all = T)
  head(GEpivot)
  GEpivot$kgNintakeDay <- ifelse(GEpivot$NinCat == "other"
                                 # , otherNintakeDay(GEpivot$DMIkgDayBW, halfInOut/10)
                                 , otherNintakeDay(GEpivot$DMIkgDayBW, 16.1)
                                 , GEpivot$kgNintakeDay)
  head(GEpivot)
  
  ## convert back to wide
  GEpivotWide <- GEpivot %>%
    pivot_wider(id_cols = c("animal", "category", "UK.Region")
                , values_from = kgNintakeDay)
  names(GEpivotWide) <- gsub("MJdayGE", "kgNinDay", names(GEpivotWide))
  head(GEpivotWide)
  
  # combine all animals
  NintakeDayCombined <- GEpivotWide
  
  # tidy
  rm(GEpivot, insideCP, alwaysOut, halfInOut, GEpivotWide)
  
  ###### 6a2 - calculate N retention ######
  ## ---------- notes ---------- # 
  ## there are two ways to calculate N retention. See Options 1 and 2
  ## Eqs. 10.31 & 10.32 in IPCC (2019)
  ## both are below
  
  ## option 1 was used for cows
  ## due to lacking other info, such as eggs produced, or piglets weaned, 
  ## other animals' N retention was calculated using option 2
  
  ### option 1 example 
  ## daily N retention is different for dairy cows compared to all else
  ## Eqs. 10.33 (IPCC, 2019)
  ## requires daily weight gain, NEg, milk amount and fat content of milk, which
  ## have all previously been calculated (see animal_attributes_create.R)
  
  ### option 2 
  # fraction of consumed N that is retained by the animal
  # Table 10.20 (IPCC) has the values used 
  ## --------------------------- # 
  
  ###### for bovine ######
  # option 1
  # from the animal attributes table, get the name, category, region, 
  # and weight of the animals
  ## ------------ Notes --------------  ##
  ## bovine retention values are independent of N intake using option 1, 
  ## instead giving a separate, single value per each animal in a category 
  ## ------------ ----- --------------  ##
  
  ## units: N retained by the animal (in kg animal-1 day-1 or year-1)
  
  # N retention calcs for cattle
  ## EQUATION 10.33 (IPCC, 2019) for dairy cows
  dairyNretent <- function(milk, fat, wg, NEg){
    
    milkPR <- (1.9 + 0.4 * fat) / 100
    left <- (milk * milkPR) / 6.38
    
    NEG <- 268 - ((7.03 * NEg) / wg)
    NEG1000 <- NEG / 1000
    right <- (wg * NEG1000)/6.25
    
    rl <- right + left
    return(rl)
  }
  
  ## for non-dairy bovine
  nondairyNretent <- function(wg, NEg){
    NEG <- 268 - ((7.03 * NEg) / wg)
    NEG1000 <- NEG / 1000
    right <- (wg * NEG1000)/6.25
    return(right)
  }
  
  ## milk - volume, fat, protein ##
  # for cows, 1 litre of milk is 1.033 kg - the litre values are found in John Nix,
  # but kg values need to be input for calculation
  # Friesians produced an average of 8,000 litres a year (Nix; p47)
  # which is equiv to 22 l a day, so 22.6411 kg day-1
  # (similar average to all British cows: https://ahdb.org.uk/dairy/uk-milk-yield)
  kgDayMilk <- 22.6411
  
  # values for fat and true protein content for the 'average' breed
  fat <- (3.9+3.9+4.2+8+4.8+3.6+3.6+5+4.5+3.8+4.2+4.5+3.7)/14
  prot <- (3.8+3.3+3.7+3.5+3.5+3.6+3.2+3.2+4.2+4+3.3+3.5+4+3.2)/14
  
  ## ------------ Notes --------------  ##
  ## for these equations, dairy and non-dairy bovines are calculated differently
  ## in terms of their N retention. The difference is due to milk production.
  ## ------------ ----- --------------  ##
  
  NretainCows <- DMIintake %>%
    dplyr::select(c("animal", "category", "UK.Region", "Weight_Kg", "wGain_Kg", "NEg")) %>%
    # split if bovine
    filter(category %in% c("other_cattle",  "calf", "heifer", "pregnant", "lactating_dairy")) %>%
    # separate into dairy and non-dairy
    mutate(kgDayNreten = if_else(grepl("lactating", category)
                                 , dairyNretent(kgDayMilk, fat, .$wGain_Kg, .$NEg)
                                 , nondairyNretent(.$wGain_Kg, .$NEg))) 
  head(NretainCows)
  
  ## combine with N intake, so that retention can be extracted from intake
  Nintake.retainCows <- NintakeDayCombined %>%
    merge(., NretainCows
          , by = c("animal", "category", "UK.Region")
          , all = F)
  head(Nintake.retainCows)
  
  ###### for non-bovine ######
  # option 2 (retained by fractions of weight)
  # see Table 10.20 (IPCC)
  ## ------------ Notes --------------  ##
  ## non-bovine retention values are dependent of N intake when using option 2, 
  ## whereby a proportion of N intake is consider excretion 
  ## ------------ ----- --------------  ##
  Nintake.retainNoBovine <- NintakeDayCombined %>%
    # split if non-bovine
    filter(!category %in% c("other_cattle",  "calf", "heifer", "pregnant", "lactating_dairy")) %>%
    # separate into different animal groups; sums indicate proportion of daily N (in kg) intake that is retained 
    mutate(kgDayNPropRetain = if_else(grepl("ewe|sheep|ram|lamb", category)
                                      , 0.1
                                      , if_else(grepl("pig|poultry", category)
                                                , 0.3
                                                , if_else(grepl("horse", category)
                                                          , 0.07, 0))))
  # times the Nintake by the fraction retained
  Nintake.retainNoBovine <- Nintake.retainNoBovine %>%
    mutate(across(contains("kgNinDay"), ~ . * kgDayNPropRetain, .names = "kgNRetainDay_{.col}"))
  # correct the names
  names(Nintake.retainNoBovine) <- sub("kgNRetainDay_kgNinDay", "kgNRetainDay_", names(Nintake.retainNoBovine))
  head(Nintake.retainNoBovine)
  
  ###### 6a3 - calculate N excretion ######
  ## ------------ Notes --------------  ##
  ## using retention and intake values for bovine, EQUATION 10.31A (IPCC, 2019)
  ## will be used to determine the difference - i.e. the excretion. For
  ## milk-producing cows, retention is always the same, regardless of freedom 
  ## ------------ ----- --------------  ##
  
  ###### for bovine ######
  Nintake.retain.excreteCows <- Nintake.retainCows %>%
    # determine difference between intake and excretion per day
    # multiply all intake columns by the single retention column
    mutate(across(contains("kgNinDay"), ~ . - kgDayNreten, .names = "kgNExcreteDay_{.col}"))
  head(Nintake.retain.excreteCows)
  
  # correct the names
  names(Nintake.retain.excreteCows) <- sub("kgNExcreteDay_kgNinDay", "kgNExcreteDay_", names(Nintake.retain.excreteCows))
  head(Nintake.retain.excreteCows)
  # check minimum of each column
  apply(Nintake.retain.excreteCows,2,min)
  
  # get yearly values based on scenario
  ## extract the two bits required
  NexcreteCows <- Nintake.retain.excreteCows %>%
    dplyr::select(c(animal, category, UK.Region
                    , kgDayNreten
                    , contains(c(
                      paste0("housed_", as.character(DEinside))
                      , paste0("sml_", as.character(DEoutside)))))) %>%
    # multiply each of the excretions by 6 months (i.e. half a year)
    ## half a year inside
    mutate(across(contains("kgNExcreteDay_housed"), ~ . * (365/2), .names = "kgNExcreteYr_housed")) %>%
    ## half a year outside
    mutate(across(contains("kgNExcreteDay_sml"), ~ . * (365/2), .names = "kgNExcreteYr_sml")) %>%
    ## sum together, to get the year total
    mutate(kgNExcreteYr = kgNExcreteYr_housed + kgNExcreteYr_sml)
  head(NexcreteCows)
  
  # tidy
  rm(Nintake.retainCows, Nintake.retain.excreteCows, NretainCows)
  
  ###### for non-bovine ######
  # times the Nintake by the fraction *NOT* retained (i.e. 1 - fraction)
  Nintake.retain.excreteNoBovine <- Nintake.retainNoBovine %>%
    mutate(across(contains("kgNinDay"), ~ . * (1-kgDayNPropRetain), .names = "kgNExcreteDay_{.col}"))
  head(Nintake.retain.excreteNoBovine)
  # correct the names
  names(Nintake.retain.excreteNoBovine) <- sub("kgNExcreteDay_kgNinDay", "kgNExcreteDay_", names(Nintake.retain.excreteNoBovine))
  head(Nintake.retain.excreteNoBovine)
  names(Nintake.retain.excreteNoBovine)
  
  # get yearly values based on scenario
  ## extract the two bits required
  NexcreteNoBovine <- Nintake.retain.excreteNoBovine %>%
    dplyr::select(c(animal, category, UK.Region
                    , contains(c(
                      paste0("housed_", as.character(DEinside))
                      , paste0("sml_", as.character(DEoutside)))))) %>%
    # multiply each of the excretions by 6 months (i.e. half a year)
    ## half a year inside
    mutate(across(contains("kgNExcreteDay_housed"), ~ . * (365/2), .names = "kgNExcreteYr_housed")) %>%
    ## half a year outside
    mutate(across(contains("kgNExcreteDay_sml"), ~ . * (365/2), .names = "kgNExcreteYr_sml")) %>%
    ## sum together, to get the year total
    mutate(kgNExcreteYr = kgNExcreteYr_housed + kgNExcreteYr_sml)
  head(NexcreteNoBovine)
  
  ## ------------ Notes --------------  ##
  ## now, both excretions for bovines and non-bovines have been calculated. They
  ## can be combined, as excretion values are the important ones for GHG
  ## ------------ ----- --------------  ##
  
  ###### bovine and non-bovine combined ######
  # combine excretions of bovines and non-bovines
  excreteNcombined <- bind_rows(NexcreteCows
                                , NexcreteNoBovine) 
  head(excreteNcombined)
  
  # get means for all of the same animal, category
  summariseN <- excreteNcombined %>%
    group_by(animal, category) %>%
    dplyr::summarise(across(where(is.numeric)
                            , list(mean = ~ mean(., na.rm = TRUE)))) %>%
    filter(grepl("dairy", category)) %>%
    dplyr::select(animal, category, contains("ExcreteYr"))
  head(summariseN)
  # should be near: 89 kg N/cow
  
  # save
  fwrite(excreteNcombined
         , file.path(savePath, "annual_excretion_rates.csv"), row.names = F)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(savePath, "readme_annual_excretion_rates.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-07-16
  Last update:",  format(Sys.Date()), "
  
  Description of annual_excretion_rates:
  The results in 'annual_excretion_rates.csv' show the annual N excretion rates in units: kg N animal-1 for each category of animal, based on different conditions.

  Columns of annual_excretion_rates.csv:
    'animal' = the category an animal was assigned to, based on either CTS or AgCensus data
    'category' = the category each animal type was assigned to in this work
    'UK.Region' = different regions of the UK, with 'other' representing non-English regions
    other column names combine different aspects:
      'kgN' = the amount of N, in kg, either retained, taken in, or excreted per animal, per time period, and per living situation
        'kgDayNreten' or 'RetainDay' = the amount of N, in kg, retained per day
        'inDay' = the amount of N, in kg, taken in per day
        'kgNExcreteDay' = the amount of N, in kg, excreted per day
        'housed' = value for animals assumed to be housed
        'sml', and 'lrg' = values for all animals with freedom to roam small or large fields, respectively
        '_[x]' = the digestible energy of an animal's diet, with 'x' being the amount in terms of percentage
    'kgNExcreteYr_housed' = the amount of N, in kg, excreted per 6-months of living in a housed situation
    'kgNExcreteYr_sml' = the amount of N, in kg, excreted per 6-months of living in a small field situation
    'kgNExcreteYr' = the final amount of N excreted, in kg, per year per animal. It is the summed amount of the 6-months housed and 6 months small field
          
  Data: per-animal N excretion
  units: kg N excreted per day or year per animal")
  
  sink(file = NULL)
  
  # tidy
  rm(summariseN, csgNintakeDay, NexcreteCows, NexcreteNoBovine, Nintake.retainNoBovine, Nintake.retain.excreteNoBovine)
  
  ## important note: N export values generally start much higher (i.e. when DE is lower), but
  # then get closer/ surpass the test values. The N intake values are also possibly high based on 
  # Table 10.19 values. The values from literature are often lower (see GDoc)
  
  ##### 6b - calculate emissions from grazing #####
  ###### 6b1 - calculate direct N emissions ######
  # only calculated for cows that are grazing - other become manure and
  # other animals' faeces does not contribute much
  # this is determined by the scenario at the top of this code
  # Eq. 3.2.16 in CFT
  # EQUATION 10.25 (UPDATED) (IPCC, 2019)
  # units: Direct N2O emissions [kg N2O] per year
  
  # exclude 'housed' as they never graze
  # therefore, due to the scenario, we are only interested in small fields
  dirNemissionsSml <- excreteNcombined %>%
    select(contains(c("animal", "UK.", "categ", "sml", "ExcreteYr"))) %>%
    select(-contains(c("kgNRetainDay", "kgNinDay"))) %>%
    filter(grepl("cattle|calf|lacta|pregnant|heifer", category))
  head(dirNemissionsSml)
  
  # we assume the turnout is half a year
  turnout <- 365/2
  # and that it is 24 hours a day
  hours <- 24
  # multiply
  timOut <- turnout * hours
  # totalTimePoss
  totalTimePoss <- 365 * 24
  timePar <- timOut / totalTimePoss
  
  # calculate Direct N2O emissions [kg N2O]
  ## multiply the average excretion by time outdoors and emission factors
  NemissionsGraz <- dirNemissionsSml %>%
    mutate(dirNemiskgN2O = kgNExcreteYr_sml * timePar * 0.02 * (44/28))
  head(NemissionsGraz)
  
  ###### 6b2 - calculate indirect N emissions ######
  ## ------------ Notes --------------  ##
  ## indirect emissions come from leaching and volatisation of excreta that
  ## is left on-field
  ## leach and volatisation calculations are similar to direct emissions
  ## Eq. 3.2.17 CFT
  ## ------------ ----- --------------  ##
  
  ###### leaching ######
  # these are emissions from leaching and runoff - the loss of water-soluble
  # plant nutrients from the soil
  
  # calculate leaching N2O emissions [kg N2O]
  ## multiply the average excretion by time outdoors and emission factors
  ## and leaching factor of 0.3
  NemissionsGraz <- NemissionsGraz %>%
    mutate(leachN_kgN2O = kgNExcreteYr_sml * timePar * 0.3 * 0.0075 * (44/28))
  head(NemissionsGraz)
  
  ###### volatisation ######
  # Volatilization is the conversion of a liquid chemical into a vapour, which escapes into the atmosphere
  # Eq. 3.2.19 CFT
  
  # calculate volatisation N2O emissions [kg N2O]
  ## multiply the average excretion by time outdoors and emission factors
  ## and Volatilisation factor (FracGASM) of 0.2
  NemissionsGraz <- NemissionsGraz %>%
    mutate(volaN_kgN2O = kgNExcreteYr_sml * timePar * 0.2 * 0.01 * (44/28))
  head(NemissionsGraz)
  
  ###### 6b3 - calculate total N emissions from grazing ######
  NemissionsGrazFinal <- NemissionsGraz %>%
    # add direct and indirect together
    mutate(grazEmisYr_kgN2O = dirNemiskgN2O + leachN_kgN2O + volaN_kgN2O) %>%
    # times all by 298 to convert to CO2e
    mutate(grazEmisYr_kgCO2e = grazEmisYr_kgN2O * 298)
  head(NemissionsGrazFinal)
  
  # select certain columns
  NemissionsGrazFinal <- NemissionsGrazFinal %>%
    dplyr::select(animal:category, dirNemiskgN2O:grazEmisYr_kgCO2e)
  head(NemissionsGrazFinal)
  
  # add back in all animals, just at 0
  NemissionsGrazFinal <- NemissionsGrazFinal %>%
    merge(., excreteNcombined
          , by = c("animal", "category", "UK.Region")
          , all = T)
  head(NemissionsGrazFinal)
  
  # save
  fwrite(NemissionsGrazFinal
         , file.path(savePath, "annual_grazing_emissions.csv"), row.names = F)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(savePath, "readme_annual_grazing_emissions.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-06-24
  Last update:",  format(Sys.Date()), "
  
  The results in 'annual_grazing_emissions.csv' show the annual CO2e emissions, in kg, from direct and indirect (leaching and volatisation) 
  The results were derived from the time the animal spent outside. They have been converted from N2O values into CO2e.
  
  Columns of annual_grazing_emissions.csv:
    'animal' = the category an animal was assigned to, based on either CTS or AgCensus data
    'category' = the category each animal type was assigned to in this work
    'UK.Region' = different regions of the UK, with 'other' representing non-English regions
    'dirNemiskgN2O' = annual direct N emissions, in kg N2O, from the excretion process per animal
    'leachN_kgN2O' = annual indirect N emissions, in kg N2O, from leaching of the excreta
    'volaN_kgN2O' = annual indirect N emissions, in kg N2O, from volatisation of the excreta
    'grazEmisYr_kgN2O' = total (direct + indirect) N emissions, in kg N2O, from grazing 
    'grazEmisYr_kgCO2e' = total (direct + indirect) N emissions, in kg CO2e, from grazing. A GWP of 298 was used to convert from N2O to CO2e 

  Data: per-animal N grazing emissions
  units: CO2e emissions per year per animal")
  
  sink(file = NULL)
  
  ##### 6c - calculate the amount of N in manure #####
  ## ------------ Notes --------------  ##
  ## The amount of N in manure is an important factor when considering fertiliser
  ## amounts, so this should be calculated per km cell. Manure is created
  ## from the excreta left when animals are in housing, so only the 'housed' 
  ## amounts of N are important
  ## ------------ ----- --------------  ##
  
  head(excreteNcombined)
  
  # the main calculation that is required is getting how many animals of each 
  # type are in a 1 km2 cell, and then multiplying by the amount of housed 
  # excreta N each animal produces
  
  # load the animal grid, which is at the 1 km2
  animalGrid <- st_read(file.path(currPath, "all_animals_1km.gpkg"))
  head(animalGrid)
  
  ## ------------ Notes --------------  ##
  ## The below section tests two methods, the first just confirming that the
  ## second (which is faster) is accurate
  ## ------------ ----- --------------  ##
  
  ###### first (test) long multiplication method ######
  # make long
  animalGridLong <- animalGrid %>%
    st_drop_geometry() %>%
    pivot_longer(!c(rcFid_1km, region))
  head(animalGridLong)
  
  # merge the two dataframes, ensuring regions match
  animalGridExcreteN <- animalGridLong %>% as.data.frame() %>%
    filter(region == "East of England") %>%
    merge(., excreteNcombined %>% 
            dplyr::select(animal:UK.Region, kgNExcreteYr_housed) %>%
            filter(UK.Region %in% c("East_of_England", "All"))
          , by.x = c("name"
                     # , "region"
          )
          , by.y = c("animal"
                     # , "UK.Region"
          )
          , all = T)
  
  head(animalGridExcreteN)
  table(animalGridExcreteN$region)
  table(animalGridExcreteN$UK.Region)
  
  # multiply N excreted by value (number of animals in that km2)
  animalGridExcreteN2 <- animalGridExcreteN %>%
    mutate(totalN = value * kgNExcreteYr_housed)
  head(animalGridExcreteN2)
  # sum by pixel
  sumExcreteN <- animalGridExcreteN2 %>%
    group_by(rcFid_1km) %>%
    summarise(total1km_kgN = sum(totalN, na.rm = T))
  head(sumExcreteN)
  
  ###### second (final) matrix multiplication method ######
  
  # add underscores to regions, to make them match output df
  animalGrid$region <- gsub(" ", "_", animalGrid$region)
  # and make London the same as South East
  animalGrid$region <- ifelse(animalGrid$region == "London", "South_East"
                              , animalGrid$region)
  # convert the humber, to match excrete N df
  animalGrid$region <- ifelse(animalGrid$region == "Yorkshire_and_The_Humber"
                              , "Yorkshire_and_the_Humber"
                              , animalGrid$region)
  unique(excreteNcombined$UK.Region)
  unique(animalGrid$region)
  
  ## get all different possible regions
  ukRegionsPossible <- unique(animalGrid$region)
  
  x <- "Highlands_and_Islands"
  ### run it through a function, extracting one region at a time
  regionsNitrogen <- pblapply(ukRegionsPossible, function(x){
    
    # get current region
    cat(x)
    currRegion <- ifelse(x %in% c("Highlands_and_Islands", "Mid_and_West_Wales"
                                  , "West_Scotland", "South_Scotland", "North_Wales"
                                  , "Mid_Scotland_and_Fife", "South_Wales_West"
                                  , "Glasgow", "Central_Scotland", "South_Wales_Central", "Lothian"
                                  , "North_East_Scotland", "South_Wales_East")
                         , "other"
                         , x)
    cat(" |", currRegion, "\n")
    
    # use region to extract from over all grid
    regionGrid <- animalGrid %>% st_drop_geometry() %>%
      filter(region == currRegion)
    
    # use the same information for the N excrete grid, but add 'All' as well
    regionN <- excreteNcombined %>% 
      dplyr::select(animal:UK.Region, kgNExcreteYr_housed) %>%
      filter(UK.Region %in% c(currRegion, "All"))
    head(regionGrid)
    head(regionN)
    
    # ensure lengths match
    stopifnot((ncol(regionGrid)-2) == nrow(regionN))
    
    # ensure colnames of the grid and row names of the N amounts match
    ## ensure that across and down are in the same positions,
    ## otherwise change their positions
    if(identical(c(names(regionGrid)[2:(ncol(regionGrid)-1)]) , c(regionN$animal))){
      regionalMix <- regionN
      cat("all names match\n")
      
    } else {
      
      # see if any names are not present
      wx <- colnames(regionGrid)[!which(colnames(regionGrid) %in% regionN$animal)]
      stopifnot(length(wx) == 0)
      
      regionalMix <- regionN[match(names(regionGrid)[2:(ncol(regionGrid)-1)]
                                   , regionN$animal), ]   
      
      # recheck - identical?
      stopifnot(identical(c(names(regionGrid)[2:(ncol(regionGrid)-1)]), c(regionalMix$animal)))
      cat("names match after reshuffle\n")
    }
    
    # multiply the animals in that region by their emissions
    ## shorten regional animals emissions
    regionalMix0 <- regionalMix
    regionalMix0[is.na(regionalMix0)] <- 0
    
    # shorten the animal numbers
    anr <- regionGrid %>%
      dplyr::select(-c(rcFid_1km, region)) %>%
      st_drop_geometry()
    names(anr)
    
    stopifnot(length(regionalMix0$animal[which(!regionalMix0$animal %in% names(anr))]) == 0)
    
    # multiply for result
    anrTimes <- as.matrix(anr) %*% diag(regionalMix0$kgNExcreteYr_housed) %>%
      as.data.frame() %>%
      mutate(totalN = rowSums(.[], na.rm = T)) 
    head(anrTimes)
    
    # run some checks
    stopifnot(anr[20, 20] * regionalMix0$kgNExcreteYr_housed[20] == anrTimes[20, 20])
    stopifnot(anr[100, 100] * regionalMix0$kgNExcreteYr_housed[100] == anrTimes[100, 100])
    stopifnot(anr[200, 200] * regionalMix0$kgNExcreteYr_housed[200] == anrTimes[200, 200])
    
    anrTimes <- anrTimes %>%
      # include id col back in
      bind_cols(regionGrid %>% dplyr::select(rcFid_1km), .)
    head(anrTimes)
    
    # compare with the test one above - for east of England
    if(currRegion == "East_of_England"){
      resultsCompare <- anrTimes %>% dplyr::select(rcFid_1km, totalN) %>%
        merge(., sumExcreteN)
      head(resultsCompare)
    }
    return(anrTimes)
  }) # end of function
  
  # bind all the individual regions, and make them spatial using rcFid1km
  regionsGridNitrogen <- bind_rows(regionsNitrogen) %>%
    dplyr::select(rcFid_1km, totalN)
  head(regionsGridNitrogen)
  ## make spatial
  spatialGridN <- ghgGridGB %>%
    merge(., regionsGridNitrogen, by = "rcFid_1km")
  head(spatialGridN)
  ### save
  st_write(spatialGridN, file.path(savePath, "spatialN1km_kgN.gpkg"), append = F)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(savePath, "readme_spatialN1km_kgN.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-06-06
  Last update:",  format(Sys.Date()), "
  
  Description of 'spatialN1km_kgN.gpkg':
  This file contains the amount of nitrogen, in kg, created in a 1 km2 by using the value of N excreted in manure (i.e. excreta collected when animals are housed).
  It is assumed that this N will be kept within a 1 km cell, and used before any additional fertiliser is used. 
  The amount of N is dependent on the diet of all the animals.
  
  Columns of 'spatialN1km_kgN.gpkg':
    'rcFid_1km' = grid reference of that cell (x and y coordinates of the centroid point separated by a'-')
    'totalN' = amount of nitrogen, in kg, excreted in a 1 km2 by all housed animals in that cell
    
  Data: total N excreted
  units: kg N per year for a 1 km2
  Spatial resolution: 1 km2
  Spatial extent: UK
  Temporal coverage: 2020
  Native projection: 27700")
  
  sink(file = NULL)
  
} # end 'part6excretions'

# tidy
rm(DEinside, DEoutside)