## ---------------------------
##
## Script name: total_ghg_script.R
##
## Purpose of script: to combine all of the different elements of greenhouse gas
##                    emissions calculation into one script                   
##                    
## Run after: ghgtool_setup.R
##            nonBovines.R
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
## 1 - create a 1 km2 GB-wide grid that will be used to record the results of all emissions / sequestrations
## 2 - determine the number of 'livestock' animals per km - this is composed of bovine and non-bovive (sheep, pigs, poultry) animals
## 3 - create an attribute table that contains all information (i.e. life history properties) about all animal types
##        This includes their roles, and the regions of GB they are found in.
## 4 - Calculate the emissions due to enteric fermentation - the process by which animals create methane
## 5 - Calculate emissions from manure management
## 6 - Calculate emissions from grazing emissions
## 7 - Calculate emissions from fertiliser use
## 8 - Calculate CO2e balance from land use
## 9 - Calculate emissions from residue management
## 10 - Calculate emissions from on-farm energy use
## 11 - combine all emissions into one final CO2e balance layer
## 13 - calculate emissions for different scenarios
##
## list of final outputs:
##    data_in/ghgGridGB.gpkg -> 1 km2 grid of inland GB that contains 'rcFid_1km', which is a unique cell identifier
##    data_in/land_cover_table.gpkg -> 
##    data_in/animals/CTSCowsGrid.gpkg ->
##    data_in/land_cover/land_cover_table_19.gpkg ->
##    data_in/animals/cow_grid_2020_1km.gpkg ->
##    data_in/animals/nonBovine_grid_2015_1km.gpkg ->
##    data_in/animals/CTSCowsGrid.gpkg ->
##    data_in/animals/all_animals_1km.gpkg ->
##    data_in/animals/animal_coefs_net_energy.csv ->
##    results/animals/mjDayGE.csv
##    results/animals/entericFermCH4CO2e.csv
##    results/animals/entericFerm_uk_CH4CO2e.csv
##    results/animals/dailyVS.csv
##    results/animals/sixMonthMMCO2e.csv / sixMonthMMCH4.csv
##    results/animals/annual_grazing_emissions.csv
##    results/animals/annual_excretion_rates.csv
##    results/animals/spatialN1km_kgN.gpkg
##    data_in/crops/KgNapplied_fert.csv/.gpkg
##    data_in/crops/sns_and_crops.csv/.gpkg
##    results/arable/liming_eval.csv / .gpkg
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
gc()

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

#### 0 - set initial paths ####
englishRegions <- "N:/Data/UK/Boundary/region/RGN_DEC_2021_EN_BUC.shp"
soilsData <- "P:/NEC07065_AGLAND/WP1/SpatialData/Soils_1km/"
scenarioResultsPath <- file.path("scenario", "results") # create scenario save path

#### 0 - run which parts? ####
# select which bits need running
## ------------ Notes --------------  ##
## This is the section where the user chooses which section they want to run by
## putting a 'T' (or 'TRUE') to the right of the specific arrows
## ------------ ----- --------------  ##
# Do you want to run all of the script?
# if you do, you just need to convert the below line 'runAll' to 'TRUE' (or 'T') - you do not need to adjust any of the other sections
runAll <- F

# part 0 - create the grid
createGrid <- F

# part 1 creates a grid that contains the composition of 1 km2 in terms of land 
# cover. The land cover information comes from CEH's LCM2015PlusCrops (LCM2015) map
part1LandCoverGrid <- F

# part 2 creates a grid that contains the number of animals present in a 1 km2
# cell. It uses Defra CTS data and AgCensus data, and Scotland and Wales maps to derive the animal numbers
part2AnimalGrid <- F

# part 3 creates attribute tables for each type of animal
part3AnimalAttributes <- F

# part 4 calculates the CO2e values from enteric fermentation for individual animals
part4Enteric <- F

# part 5 calculates the amount of manure per animal and the emissions from its management 
part5ManureManagement <- T

# part 6 calculates the amount of excretions produced per animal
part6excretions <- F
## scenarios as part of that?
nManureScens <- F

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

# combine all the final emissions together
part11combineemissions <- T

# run all the calculations for the scenarios as well
runAllScenarios <- F

if(runAll){
  createGrid <- part1LandCoverGrid <- part2AnimalGrid <- T
  part3AnimalAttributes <- part4Enteric <- part5ManureManagement <- T
  part6excretions <- part7fertiliserUse <- part8landuse <- T
  part9residue <- part10onfarmenergy <- part11combineemissions <- T
  runAllScenarios <- T
}

# Filter the objects to keep
objKeep <- ls()[grep("part", ls())]

#### 0 - select animal scenario ####
## ------------ Notes --------------  ##
## This sections allows you to choose what the set-up of animals will be. This
## relates to the quality of food, for both the inside portion (6 months) and
## the outside portion (6 months). The range should be between 60 and 85% 
## digestible energy (DE). The defaults below, x and y, are recommended as that
## produces enteric fermentation and volatile solids (in manure management)
## similar to what literature suggests (for the UK)
## ------------ ----- --------------  ##
DEinside <- 85 # 'housed' value
DEoutside <- 79 # outside value
inDE <- paste0("housed_", DEinside)
outDE <- paste0("sml_", DEoutside)

#### 0 - select global warming potentials (GWP) ####
## ------------ Notes --------------  ##
## GWPs have changed over the years, and can converted including feedback or not, 
## and across different time periods.
## Below, we have used, following AR6 (2014), GWPs across a 100-year time period,
## and with no feedback included. Therefore, GWPs used were:
## CO2 = 1
## CH4 = 28
## N2O = 265
## ------------ ----- --------------  ##
gwpCH4 <- 28 
gwpN2O <- 265

#### 0 - run which scenarios? ####
## ------------ Notes --------------  ##
## In this section, the user should select which scenarios to run by providing
## paths to the polygons that contain the relevant information on animal
## numbers and land cover composition. 
## Note: for this example, scenarios were produced as part of the 'agland_ghg_scenarios.R'
## script, focusing on expansion of arable or grassland
## Note also: the first three sections of this script focus on creating the initial
## 2015 data map. Any further scenarios introduced in this section will be considered 
## from part 4 onwards
## ------------ ----- --------------  ##
scenPath <- file.path("scenario", "scen_maps")
scenarios <- c(
  # write paths here (with a comma at the start, with the exception of path 1)
  file.path(scenPath, "Ag_15.gpkg") # path 1
  , file.path(scenPath, "Ag_30.gpkg") # path 2
  , file.path(scenPath, "baseline2007.gpkg") # path 3
  , file.path(scenPath, "Gr_15.gpkg") # path 4
  , file.path(scenPath, "Gr_30.gpkg") # path 5
)

# if you need to separate them into animals and land cover, put 'T' below
sepAnimalsLand <- F
# if the above is 'T', this is carried out between parts 2 and 3

# load regions of England
st_read(englishRegions) %>%
  # save locally as gpkg
  st_write("data_in/RGN_DEC_2021_EN_BUC.gpkg", append=F)
regions <- st_read("data_in/RGN_DEC_2021_EN_BUC.gpkg")

#### 0 - load the GB grid ####
## ------------ Notes --------------  ##
## The GB grid should have been created in the 'ghgtool_setup.R' script, and
## saved as 'agland_grid_ESW.gpkg'. That grid needs to be loaded. The script 
## below modifies it as to only keep the 1 km2 of inland GB.
## ------------ ----- --------------  ##
tic("Time taken to create a grid")
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
            , 454677, 454052, 455301, 265407, 265419, 296532, 295907, 294030, 294656
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
toc(log = T)

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
  stopifnot(which(!lc1km$rcFid_1km %in% ghgGridGB$rcFid_1km) == 0)
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
  
  # stop("part 2 - top")
  
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
    
    notess
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
    
    ## ------------ Notes --------------  ##
    ## on 6 cores the below took approximately
    ## 7 minutes
    ## ------------ ----- --------------  ##
    tic("total extraction time for grass")
    # Set up parallel processing (adjust the number of cores as needed)
    registerDoParallel(cores = 6)
    
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
                         cowsAndGrass <- bind_cols(grassRow %>% st_drop_geometry() %>%
                                                     dplyr::select(rcFid_1km), cowsAndGrass)
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
    apply(cAndG %>% st_drop_geometry(), 2, function(x) max(x, na.rm = TRUE))
    
    # merge with GB grid
    cows1kmGB <- merge(ghgGridGB, cAndG, by = "rcFid_1km", all = T)
    cows1kmGB[is.na(cows1kmGB)] <- 0
    head(cows1kmGB[, 1:20])
    plot(cows1kmGB[which(names(cows1kmGB) == "Continental_0 - 3 months_Females for slaughter")])
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
  ## We could only get bovine data from Defra, which covered England, Scotland,
  ## and Wales. Data on other animals (for England only) came from
  ## AgCensus, of which data is available at 5 km2 resolution. These data included
  ## numbers of sheep, pigs, and poultry.
  ## Assuming these animals are always on grasslands, a similar method will be
  ## used for these animals as it was for the bovines: they were calculated per
  ## grassland proportion. As the AgCensus data were from 2010, the numbers 
  ## needed to be updated to 2020 i.e. to match the cattle data
  ## This was done using the total number of pigs, sheep, and poultry between
  ## 2010 and 2020, and adjusting the proportion appropriately. In practice, this
  ## meant that each animal had more space / less space per individual depending on whether
  ## that animal type increased or decreased.
  
  ## In England, pigs decreased by 26% between 2010 and 2020, while sheep decreased 21%
  ## and poultry increased by 7%
  
  ## Non-bovine data for Scotland and Wales were not 
  ## available from the AgCensus dataset. These data were therefore derived 
  ## from government datasets that were released in the form of spatial images. 
  ## The non-bovine animal groups obtained in this way were sheep, pigs, and 
  ## poultry. Each spatial image contained group data on the number of maximum 
  ## capacities of the registered holdings of that type per km2. The maps were 
  ## produced using a smoothing kernel density function based on the holdings. 
  ## Each of the maps were created as raster after the were georeferenced 
  ## based on the nation outlines. The '2b5b' section below deals with extracting
  ## the numbers of these animals to Scotland and Wales. This is based on the 
  ## maps and government data collected in 2020. 
  ## ------------ ----- --------------  ##
  
  # create if it does not already exist - the non-bovine grid for England
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
      ## should take between 1h 20 mins - 1h 40 mins
      ## ------------ ----- --------------  ##
      tic("total time for non-bovine grid")
      registerDoParallel(cl <- makeCluster(6))
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
      apply(nonBov1km %>% st_drop_geometry(), 2, function(x) max(x, na.rm = TRUE))
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
    plot(regions.WS[1])
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
    # plot(cellRegions[2])
    
    # save
    fwrite(cellRegions %>% st_drop_geometry()
           , file.path("data_in", "regions_assigned.csv"), row.names = F)
    
    ###### 2b5 - proportional relationship of non-bovine animals between 2010 and 2020 ######
    ## ------------ Notes --------------  ##
    ## in this section, proportional relationships between non-bovine animals recorded in
    ## 2010 (i.e. same time as AgCensus data) and 2020 will be used to adjust the 
    ## original 2010 AgCensus numbers. This is only for England. For Scotland and Wales, 
    ## see the next sections (2b5b, 2b5c)
    
    ## Pigs increased by 12% between 2010 and 2020, while sheep increased 5%
    ## and poultry increased by 11%
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
      mutate(across(contains("pigs"), ~.*1.12)
             ## sheep
             , across(contains(c("ewe", "sheep", "rams")), ~.*1.05)
             ## sheep
             , across(contains("poultry"), ~.*1.11))
    head(nonBov1kmGB.prop2)
    
    st_write(nonBov1kmGB.prop2, file.path(currPath, "nonBovine_grid_2020_1km.gpkg"), append = F)
    
    # write readme
    sink(file = NULL)
    sink(file = file.path(currPath, "readme_nonBovine_grid_2015_1km.md"))
    cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-08-10
  Last update:",  format(Sys.Date()), "
  
  Description of 'nonBovine_grid_2020_1km.gpkg':
  This is a gridded dataset that contains the quantity of different farmed non-bovine animals, based on the underlying grass area. 
  The animals are split up into different species and life stages e.g. breeding ewes intended for further breeding
  These data were obtained from 2015 AgCensus, which was at the 5 km2 resolution and featured data collected in 2010
  Using the amount of grass in each 1 km2 underlying the each 5 km2, the non-bovine numbers in a 5 km2 cell where proportied to each 1 km2 based on proportion of grass
  These numbers were then further adjusted based on changes in animals numbers between 2010 and 2020: Pigs decreased by 26% between 2010 and 2020, while sheep decreased 21% and poultry increased by 7%
  
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
  ## the lowest number of each bracket was used 
  ## These data were geoprocessed from the original, polygonised.
  ## This was calculated in 'nonBovines.R'
  ## Below, they were loaded back in
  ## ------------ ----- --------------  ##
  
  # stop("grass roots")
  # rm(list=setdiff(ls(), c("currPath")))
  # create if it does not already exist - the non-bovine grid for Scotland and Wales
  if(file.exists(file.path(currPath, "nonBovine_SW_2015_1km.gpkg"))){
    nbScotWal <- st_read(file.path(currPath, "nonBovine_SW_2015_1km.gpk"))
  } else {
    
    # read in polygon
    nonBov1kmSW <- st_read(file.path("results", "animals", "nonBov1kmSW.gpkg"))
    nonBov1kmSW
    
    #### pig proportions
    # extract to polygon
    nonBov1kmGBprop2 <- st_read("nonBov1kmGBprop2.gpkg")
    head(nonBov1kmGBprop2)
    ## use the proportions of different roles of different pigs and sheep to determine those in scotland
    nonBov1km.spl <- nonBov1kmGBprop2 %>% 
      st_drop_geometry() %>%
      filter(total.pigs..c47. > 0) %>%
      dplyr::select(contains(c("breeding.pig", "fattening.pig"))) %>%
      ### get totals, and then proportion of each categories
      mutate(total = rowSums(., na.rm = T))
    # Calculate percentages
    percentages <- nonBov1km.spl %>%
      mutate_at(vars(-total), ~ . / total * 100)
    # get the final sums
    percentages2.pigs <- colMeans(percentages[, 1:(ncol(percentages)-1)])
    ### breeding pigs make up 0.9133 of all pigs
    
    #### sheep proportions
    ## use the proportions of different roles of different pigs and sheep to determine those in scotland
    nonBov1km.spl <- nonBov1kmGBprop2 %>% 
      st_drop_geometry() %>%
      filter(total.pigs..c47. > 0) %>%
      dplyr::select(contains(c("ewe", "sheep", "lamb", "ram"))) %>%
      dplyr::select(-total.sheep.and.lambs..c44.) %>%
      ### get totals, and then proportion of each categories
      mutate(total = rowSums(., na.rm = T))
    # Calculate percentages
    percentages <- nonBov1km.spl %>%
      mutate_at(vars(-total), ~ . / total * 100)
    # get the final sums
    percentages2.sheep <- colMeans(percentages[, 1:(ncol(percentages)-1)])
    percentages2.sheep
    ### breeding ewes make up 0.54 and 0.36, first time breeders are 0.03, non-breeding sheep are 0.01
    ### and rams for service are 0.06 
    
    head(nonBov1kmSW)
    # multiply, to get final values for different groups and roles
    ## for pigs
    nonBov1kmSW.pigs <- bind_cols(lapply(1:length(percentages2.pigs), function(j){
      nonBov1kmSW2 <- percentages2.pigs[j] / 100 * nonBov1kmSW$pig_per_km
    })) %>% as.data.frame() %>%
      setNames(names(percentages2.pigs))
    
    ## for sheep
    nonBov1kmSW.sheep <- bind_cols(lapply(1:length(percentages2.sheep), function(j){
      percentages2.sheep[j] / 100 * nonBov1kmSW$sheep_per_km
    })) %>% as.data.frame() %>%
      setNames(names(percentages2.sheep))
    
    # st_write(nonBov1km.part2, "nonBov1km2.gpkg", append=F)
    
    ###### 2b5c - proportional relationship of non-bovine animals - Poultry
    ## ------------ Notes --------------  ##
    ## according to Scotland Gov agricultural survey results [https://www.gov.scot/publications/results-december-2020-agricultural-survey/]
    ## there were 6.68 million poultry for meat production in 2020, and 6.6 million for egg production,
    ## which equals 13.28 together.
    
    ## for Wales, there were 2,438,761, 1,787,890, and 5,326,057 poultry for 
    ## North Wales, South Wales, Mid and West Wales, respectively (9,552,708 combined). 
    ## Scotland and Wales together make 22.55 million poultry
    
    ## like sheep and pigs, poultry was calculated using georeferenced maps, with 
    ## the overall amounts limited to the values above. In the df below, highest value of
    ## a range were used. The ranges were:
    ## 0-50, 51-250, 251 - 600, 601-1300, 1301-2500, 2501-39474 for poultry
    ## http://apha.defra.gov.uk/documents/surveillance/diseases/lddg-pop-report-avian2020.pdf
    ## ------------ ----- --------------  ##
    
    ### combine the values of all animals
    nonBov1kmSW.end <- bind_cols(nonBov1kmSW %>% dplyr::select(rcFid_1km, poul_per_km) %>%
                                   rename("total.poultry..c48." = poul_per_km)
                                 , nonBov1kmSW.pigs
                                 , nonBov1kmSW.sheep)
    head(nonBov1kmSW.end)
    
  } # end of non-bovine grid
  
  ##### 2b6 - add non-bovine animals back in #####
  ## ------------ Notes --------------  ##
  ## There are four important elements for below:
  ## non-bovine from Scotland and Wales
  ## non-bovine from England
  ## bovine from Scotland and Wales
  ## bovine from England
  
  ## get all those dfs separately
  ## ------------ ----- --------------  ##
  
  ## non-bovine from Scotland and Wales
  nbSW <- nonBov1kmSW.end
  plot(nbSW[2])
  ## non-bovine from England
  nbE <- st_read(file.path(currPath, "nonBovine_grid_2020_1km.gpkg")) %>%
    # restrict to England, and non-bovine
    filter(!rcFid_1km %in% nbSW$rcFid_1km) %>%
    dplyr::select(any_of(c(names(nbSW))))
  plot(nbE[2])
  ## bovine from Scotland and Wales
  bSW <- st_read(file.path(currPath, "cow_grid_2020_1km.gpkg"))  %>%
    # restrict to non-England, and bovine
    filter(rcFid_1km %in% nbSW$rcFid_1km)
  names(bSW)
  ## bovine from England
  bE <- st_read(file.path(currPath, "cow_grid_2020_1km.gpkg"))  %>%
    # restrict to England
    filter(!rcFid_1km %in% nbSW$rcFid_1km)
  names(bE)
  
  ## bind cows and non-bovine
  nbAll <- bind_rows(nbSW, nbE)
  bAll <- bind_rows(bSW, bE)
  stopifnot(length(which(names(nbAll) %in% names(bAll))) == 2)
  
  ### merge both
  animalsFinal <- merge(bAll
                        , nbAll %>% st_drop_geometry()
                        , by = "rcFid_1km"
                        , all = T) %>% filter(!st_is_empty(.))
  head(animalsFinal)
  sort(unique(animalsFinal$region))
  class(animalsFinal)
  apply(animalsFinal %>% st_drop_geometry(), 2, function(x) max(x, na.rm = TRUE))
  
  # save 
  st_write(animalsFinal, file.path(currPath, "all_animals_1km.gpkg"), append = F)
  
  # tidy
  rm(bSW, bE, nbSW, nbE,nonBov1kmGB.prop2, nonBov1kmGB2, nbAll, bAll
     , nonBov1kmSW
     , percentages, percentages2.pigs, percentages2.sheep)
  gc()
  
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

if(sepAnimalsLand){
  
  # load in created animal and land cover grids
  # set location of files that relate to animal data
  currPath <- file.path("./data_in", "animals")
  # load in animal numbers
  animalNumbers <- st_read(file.path(currPath, "all_animals_1km.gpkg"))
  names(animalNumbers)
  
  # set location of files that relate to animal data
  currPath <- file.path("./data_in", "land_cover")
  # load in animal numbers
  landArea <- st_read(file.path(currPath, "land_cover_table.gpkg"))
  names(landArea)
  
  # combine the names, for the purposes of selecting those from the scenario data
  keepNames <- c(names(animalNumbers), names(landArea))
  
  # create saving directory for separate land and animal aspects
  dir.create(file.path(scenPath, "finer_detail"), showWarnings = F)
  
  for(i in 1:length(scenarios)){
    
    cat("i =", i, scenarios[[i]], "\n")
    
    # load in full dataset
    fulldf.in <- st_read(scenarios[[i]], quiet = T) %>%
      # keep only the columns with the names from the original data
      dplyr::select(all_of(keepNames))
    
    names(fulldf.in)[which(!names(fulldf.in) %in% keepNames)]
    names(fulldf.in)
    
    # split into land area and animals separately keeping ID
    ## land
    fulldf.land <- fulldf.in %>%
      dplyr::select(rcFid_1km, contains("_ha")) %>%
      replace(is.na(.), 0)
    ## animals
    fulldf.animals <- fulldf.in %>%
      dplyr::select(-c(contains("_ha"))) %>%
      replace(is.na(.), 0)
    
    ### save those as separate
    cat("saving ...", basename(gsub(".gpkg", paste0("_land.gpkg"),  scenarios[[i]])), "...\n")
    st_write(fulldf.land, file.path(scenPath, "finer_detail"
                                    , basename(gsub(".gpkg", paste0("_land.gpkg"),  scenarios[[i]])))
             , append = F)
    st_write(fulldf.animals, file.path(scenPath, "finer_detail"
                                       , basename(gsub(".gpkg", paste0("_animal.gpkg"),  scenarios[[i]])))
             , append = F)
  }
}

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
                                                      # cows are impregnated around 15 months
                                                      , ifelse(grepl("15...18|18...21|21...24|24...27|27...30|30...33|airy", cowNm$animal) & !grepl("_Beef", cowNm$animal), "lactating_dairy"
                                                               , "other_cattle")))))
  ## see all combinations
  uc <- cowNm %>%
    group_by(animal, category) %>%
    reframe(n = n())
  
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
  ## growing cattle, respectively, the values of DMI per day were slightly too low
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
  
  # stop("dmi calcs")
  
  ## ------------ Notes --------------  ##
  ## this section shows previous examples for growing cattle and calves, based on IPCC equations
  ## This was not used in the final methodology
  ## It assumes a NEmf of 7.5, based on high quality diet from Table 10.8a
  
  NEmf <- 7.5
  # calculate DMI for calves
  DMIkgDay.calves <- DMIinput %>%
    filter(category == "calf") %>%
    dplyr::select(animal, category, Weight_Kg, UK.Region) %>%
    mutate(DMI = Weight_Kg^0.75 * 
             (
               (0.0582 * NEmf - 0.00266 * NEmf^2 - 0.0869) # top line (numerator) 
               /
                 (0.239 * NEmf) # denominator
             )
           # determine DMI as percentage of bodyweight
           , BWDMI = DMI / Weight_Kg
    )
  head(DMIkgDay.calves)
  # determine range per cent of BW
  range(DMIkgDay.calves$BWDMI)
  print(range(DMIkgDay.calves$BWDMI))
  
  # calculate DMI for growing cattle
  DMIkgDay.grow <- DMIinput %>%
    dplyr::select(animal, category, Weight_Kg, UK.Region) %>%
    filter(category == "other_cattle") %>%
    mutate(DMI = Weight_Kg^0.75 * 
             (
               (0.0582 * 7.5 - 0.00266 * 7.5^2 - 0.1128) # top line (numerator) 
               /
                 (0.239 * 7.5) # denominator
             )
           # determine DMI as percentage of bodyweight
           , BWDMI = DMI / Weight_Kg
    )
  print(range(DMIkgDay.grow$BWDMI))
  ## ------------ ----- --------------  ##
  
  ### ----- final used methodology ---- ###
  # calculate DMI - based on specific equations, or bodyweight for non-bovines or growing cattle
  DMIinput$DMIkgDay <- ifelse(grepl("calf", DMIinput$category), DMIinput$Weight_Kg * 0.021
                              # lambs (4%)
                              , ifelse(grepl("lamb", DMIinput$animal), DMIinput$Weight_Kg * 0.04
                                       # rams (1.5%)
                                       , ifelse(grepl("ram", DMIinput$animal), DMIinput$Weight_Kg * 0.015
                                                # ewes ('dry')
                                                , ifelse(DMIinput$animal == "breeding.ewes.intended.for.slaughter..c39.", DMIinput$Weight_Kg * 0.015
                                                         # ewes (young)
                                                         , ifelse(DMIinput$animal == "ewes.intended.for.first.time.breeding..c40.", DMIinput$Weight_Kg * 0.03
                                                                  # ewes (mid)
                                                                  , ifelse(grepl("ewe", DMIinput$animal), DMIinput$Weight_Kg * 0.0275
                                                                           # heifers (EQUATION 10.18A)
                                                                           , ifelse(DMIinput$category == "heifer", 3.184 + 0.01536 * DMIinput$Weight_Kg * 0.96
                                                                                    # bulls/steers (EQUATION 10.18A)
                                                                                    , ifelse(grepl("teer|ull", DMIinput$animal), 3.83 + 0.0143 * DMIinput$Weight_Kg * 0.96
                                                                                             # lactating dairy cows (EQUATION 10.18B)
                                                                                             , ifelse(grepl("lactating|pregn", DMIinput$category), 0.0185 * DMIinput$Weight_Kg + 0.305  * FCM
                                                                                                      # EQUATION 10.18
                                                                                                      , ifelse(grepl("other_cat", DMIinput$category), DMIinput$Weight_Kg * 0.021
                                                                                                               # pigs (2%)
                                                                                                               , ifelse(grepl("pig", DMIinput$category), DMIinput$Weight_Kg * 0.02
                                                                                                                        # 2.5% bodyweight assumed for poultry DMI
                                                                                                                        , 0.07 * DMIinput$Weight_Kg)))))))))))
  # check the %s of BW
  pcWeight <- (DMIinput$DMIkgDay/DMIinput$Weight_Kg)*100
  # tidy
  rm(pcWeight)
  head(DMIinput)
  
  # run checks to ensure DMI is sensible
  ## for dairy cows, it should be between 16 - 20 kg (based on Feed into Milk, 2004)
  DMIinput %>%
    filter(category %in% c("lactating_dairy", "pregnant")) %>%
    dplyr::select(DMIkgDay) %>%
    range()
  
  ## for calves
  DMIinput %>%
    filter(category %in% c("calf")) %>%
    dplyr::select(DMIkgDay) %>%
    range()
  
  DMIinput %>%
    filter(category %in% c("lactating_dairy", "pregnant")) %>%
    dplyr::select(DMIkgDay) %>% unlist %>% as.numeric() %>% mean(na.rm = T)
  ## for heifer
  DMIinput %>%
    filter(category == "heifer") %>%
    dplyr::select(DMIkgDay) %>%
    range()
  ## for beef cows, 11.5 kg should be max (based on AFRC, 1993)
  DMIinput %>%
    filter(grepl("beef|Beef", animal)) %>%
    dplyr::select(DMIkgDay) %>%
    range()
  ### cap beef DMI at 11.5
  DMIinput$DMIkgDay <- ifelse(grepl("beef|Beef", DMIinput$animal) & DMIinput$DMIkgDay > 11.5
                              , 11.5, DMIinput$DMIkgDay)
  DMIinput %>%
    filter(grepl("beef|Beef", animal)) %>%
    dplyr::select(DMIkgDay) %>%
    range()
  ## get median for all categories
  medmed <- DMIinput %>%
    group_by(animal, category) %>%
    reframe(median = median(DMIkgDay)
            , max = max(DMIkgDay))
  
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
Produced by 'total_ghg_script.R'

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
          'NEm' = Net energy required for maintenance, in MJ day-1 
          'NEa_[x]' = net energy for activity based on an animal’s feeding situation (x) for a typical animal in this category and region combination, with 'x' being 'housed' (housed all the time), 'flat_sml_past' (able to move across a small, flat pasture), and 'hill_lrg_past' (able to roam across a large, hilly, pasture)
          'NEg' = net energy required for growth, in MJ day-1 
          'NEl' = net energy required for lacatation, in MJ day-1 
          'NEw' = net energy required for wool production, in MJ day-1 
          'NEp' = net energy required for pregnancy, in MJ day-1
")
  
  sink(file = NULL)
  
  # tify
  rm(part3AnimalAttributes)
  
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
  units: MJ per day")
  
  sink(file = NULL)
  # tidy
  rm(GEremHoused, GEreg, GEremLrg, GEremregHoused, GEremregSml, GEremSml
     , GEremregLrg, de)
  
  ## ------------ Notes --------------  ##
  ## check to ensure that it is about right i.e. 1.5 percent and 3.0 percent of the animal’s weight.
  ## to check  the estimate can be converted in daily intake in kilograms by dividing by 18.45 MJ/kg.
  ## check middling column.
  ##
  ## Gibb (2002): 'To check the estimate of daily gross energy intake from Equation 14, the estimate can be converted in daily
  ## intake in kilograms by dividing by 18.45 MJ/kg. This estimate of intake in kilograms should generally be
  ## between 1.5 percent and 3.0 percent of the animal’s weight'
  ## ------------ ----- --------------  ##
  
  GEcheck <- as.data.frame(cbind(animal = netEnergyTable$animal
                                 , category = netEnergyTable$category
                                 , UK.Region = netEnergyTable$UK.Region
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
                                             , if_else(GE85 > weightTo, "NO - too high", "NO")))) %>%
    mutate(difference65 = ifelse(fitCriteria65 != "Yes", (abs(GE65-weightFrom)/actualWeight)* 100, 0)) %>%
    mutate(difference85 = if_else(fitCriteria85 != "Yes", (abs(GE85-weightFrom)/actualWeight)* 100, 0)) %>%
    group_by(animal, UK.Region) %>%
    summarise(cc = mean(difference85))
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
      # select just housed 85 and small field 79
      subs <- dplyr::select(., MJdayGEhoused_85, MJdayGEsml_79)
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
  
  # stop("testing limits")
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
  table(colSums(dairyGEcheck2[, 3:ncol(dairyGEcheck2)]))
  table(as.matrix(dairyGEcheck2[, 3:ncol(dairyGEcheck2)]))
  which(colSums(dairyGEcheck2[, 3:ncol(dairyGEcheck2)]) > 23)
  cn <- names(dairyGEcheck2 %>% ungroup() %>%
                dplyr::select(MJdayGEhoused_65:MJdayGElrg_85))[which(colSums(dairyGEcheck2[, 3:ncol(dairyGEcheck2)]) > 23)]
  # save the names that meet the criteria
  criteria.match <- cn
  save(criteria.match, file = "critMatch.RData")
  
  # it appears that the general trend is that young animals e.g. cows < 6 months (i.e., heifers)
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
  
  for(de in seq(65, 85, 1)){
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
    # convert to CO2e (using non-fossil origin AR6 2014 value: 28)
    kgCO2GE <- kgCH4GE * gwpCH4
    
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
    mutate(fitCriteria = if_else(Min <= 125.44 & Max >= 125.44
                                 , "Yes"
                                 , "NO"))
  head(dairyCH4)
  
  # refine to those that met the criteria earlier
  dairyCH4.crit <- dairyCH4 %>%
    # remove heifers
    filter(category != "heifer") %>%
    dplyr::select(contains(criteria.match))
  ## max min
  dairyCH4.crit$min <- do.call(pmin, dairyCH4.crit)
  dairyCH4.crit$max <- do.call(pmax, dairyCH4.crit)
  dairyCH4.crit <- dairyCH4.crit %>%
    mutate(fitCriteria = if_else(min <= 125.44 & max >= 125.44 #covers the range
                                 , "Yes"
                                 , "NO"))
  
  # repeat based on current DE used
  dairyCH4.critCurr <- dairyCH4 %>%
    # remove heifers
    filter(category != "heifer") %>%
    dplyr::select(contains(c("housed_85", "sml_79")))
  ## max min
  dairyCH4.critCurr$min <- do.call(pmin, dairyCH4.critCurr)
  dairyCH4.critCurr$max <- do.call(pmax, dairyCH4.critCurr)
  dairyCH4.critCurr <- dairyCH4.critCurr %>%
    mutate(fitCriteria = if_else(min <= 125.44 & max >= 125.44 #covers the range
                                 , "Yes"
                                 , "NO"))
  
  # for dairy cows, 125.44 kg CH4/head/year (GHGi, 2019) needs to be covered, which it is
  # for dairy heifers, 53.20 kg CH4/head/year (GHGi, 2019) needs to be covered, which it is
  stopifnot(length(unique(dairyCH4.crit$fitCriteria)) == 1)
  
  # use the scenario described in the '0 - select animal scenario' section 
  # to give a per-animal enteric fermentation emission
  
  # as we are going to use the scenario of 50% time housed with 83% DE and 50% grazing with 79% DE
  # see whether the half and half is ok
  
  # stop("before enteric trial")
  ## ------------ Notes --------------  ##
  ## you can comment out this section to test other values
  # DEinside <- 81 # 'housed' value
  # DEoutside <- 79 # outside value
  # inDE <- paste0("housed_", DEinside)
  # outDE <- paste0("sml_", DEoutside)
  ## ------------ ----- --------------  ##
  
  dairyCH4v2 <- dairyCH4 %>%
    dplyr::select(animal, category
                  , contains(c(inDE, outDE)))
  ## should be four columns. If so, rename
  if(ncol(dairyCH4v2) == 4){
    names(dairyCH4v2)[3:4] <- c(paste0("kgCH4GE.MJday", c("in", "out"))) 
    names(dairyCH4v2)
  } else {
    stop("wrong columns in 4a2")
  }
  dairyCH4v2 <- dairyCH4v2 %>%
    mutate(kgCH4GE.50pcScen = (kgCH4GE.MJdayout/2) + (kgCH4GE.MJdayin/2)) %>%
    rowwise() %>%
    mutate(litDifference = kgCH4GE.50pcScen - 125.44)
  head(dairyCH4v2)
  dairyCH4v3 <- dairyCH4v2 %>%
    filter(category != "heifer")
  mean(dairyCH4v3$litDifference)
  
  # see with all possibilities
  dairyCH4v3 <- dairyCH4 %>%
    dplyr::select(-c(Min, Max, fitCriteria)) %>%
    mutate(across(3:last_col(), ~ ./2)) %>%
    # refine to just the criterai matched ones
    dplyr::select(contains(c(criteria.match)))
  
  # stop(dairy")
  
  # although lambs always housed
  CH4EntFerm <- DEtable2 %>%
    dplyr::select(animal, category, UK.Region
                  , contains(c(inDE, outDE)))
  ## should be five columns. If so, rename
  if(ncol(CH4EntFerm) == 5){
    names(CH4EntFerm)[4:5] <- c(paste0("kgCH4GE.MJday", c("in", "out"))) 
    names(CH4EntFerm)
  } else {
    stop("wrong columns in 4a2")
  }
  CH4EntFerm <- CH4EntFerm %>%
    rowwise() %>%
    mutate(kgCH4yr.50pcScen = ifelse(grepl("lamb", category), kgCH4GE.MJdayin
                                     , (kgCH4GE.MJdayout/2) + (kgCH4GE.MJdayin/2))) %>%
    # convert to CO2e (non-fossil fuel)
    mutate(kgCO2.50pcScen = kgCH4yr.50pcScen * gwpCH4)
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
  
  ##### 4b - CH4 based on UK-specific equation #####
  # stop("uk time")
  ## ------------ Notes --------------  ##
  ## the calculations below follow the methodology of Brown et al. (2021), where
  ## bovine, and ovine CH4 were calculated using equations based on dry matter intake,
  ## whereas other animals (i.e., pigs) used a Tier 1 approach from IPCC
  ## ------------ ----- --------------  ##
  
  # load in table with Net Energy amounts per animal category
  animAttrib <- fread(file.path(currPath, "animal_coefs_net_energy.csv")) %>%
    as.data.frame()
  unique(animAttrib$category)
  
  # units = gCH4 day-1
  animAttribCH4 <- animAttrib %>%
    dplyr::select(1:3, DMIkgDay) %>%
    # if dairy cattle
    mutate(ch4_g_day = ifelse(category %in% c("lactating_dairy", "pregnant", "heifer")
                              , (15.8185 * DMIkgDay) + 88.6002
                              # other cattle
                              , ifelse(category %in% c("other_cattle", "calf")
                                       , (17.5653 * DMIkgDay) + 45.8688
                                       # sheep
                                       , ifelse(category %in% c("ewes", "lambs", "other_sheep", "rams")
                                                , (12.3894 * DMIkgDay) + 5.1595
                                                # pigs (Tier 1)
                                                , ifelse(category %in% c("pigs")
                                                         , 1.5 / 365  
                                                         # poultry (Tier 1)
                                                         , 0)))))
  head(animAttribCH4)
  
  # get annual values, in kg
  animAttribCH4Year <- animAttribCH4 %>%
    mutate(ch4_kg_year = (ch4_g_day * 365)/1000) %>%
    # convert to co2e
    mutate(co2e_kg_year = ch4_kg_year * gwpCH4)
  head(animAttribCH4Year)
  
  # compare with the longer method
  ch4Comp <- merge(animAttribCH4Year
                   , CH4EntFerm %>% dplyr::select(1:3, kgCO2.50pcScen))
  head(ch4Comp)
  
  # save final per-animal amount of CH4 and CO2e enteric fermentation
  fwrite(animAttribCH4Year
         , file.path(savePath, "entericFerm_uk_CH4CO2e.csv")
         , row.names = F)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(savePath, "readme_entericFerm_uk_CH4CO2e.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-02-13
  Last update:",  format(Sys.Date()), "

  Files:
    entericFerm_uk_CH4CO2e.csv
  
  Columns of entericFerm_uk_CH4CO2e.csv:
      'animal' = the category an animal was assigned to, based on either CTS or AgCensus data
      'category' = the category each animal type was assigned to in this work
      'UK.Region' = different regions of the UK, with 'other' representing non-English regions
      'DMIkgDay' = the amount of dry matter intake (DMI) an animal requires in kg per day, based on the weight of a typical animal in this category and region combination
      'ch4_g_day' = CH4 emissions per day, in g, per animal type, based on DMI, calculated using Brown's equations, and multiplied from ch4_g_day
      'ch4_kg_year' = CH4 emissions per year, in kg, per animal type, based on DMI, and calculated using Brown's (2021) equations
      'co2e_kg_year' = the final amount of per-animal enteric fermentation emissions, in kg CO2e, based on Brown (2021) equations
    
  Data: per-animal enteric fermentation emissions, in kg CO2e or kg cH4
  units: kg CO2e or kg cH4 per year")
  sink(file = NULL)
  
  # stop("before end part 4")
  
  # tidy
  rm(part4Enteric)
  
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
  
  ##### 5.0 - table for manure storage fraction / mcf #####
  ## ------------ Notes --------------  ##
  ## for MM_table:
  ##      MCF values are from Brown et al. (2021)
  ##      EF [emission factors, in g CH4 kg VS-1] are from IPCC (2019), TABLE 10.14
  ## Stor_frac [storage fraction, based on storage type / Animal Waste Management Systems (%)]
  ## are from Brown et al. (2021)
  ## ------------ ----- --------------  ##
  
  MM_data <- tibble(
    factor = c("MCF", "EF"),
    awms_slurry = c(17, 33.8),            # slurry is liquid in Brown et al. (2021)
    awms_daily = c(0.1, 0.2),
    awms_FYM_cows_pigs = c(17, 3.2),      # deep bedding / farmyard manure
    awms_FYM_sheep = c(2, 3.2),          # deep bedding / farmyard manure
    awms_past = c(0.08, 1.6),            # pasture / dry lot
    awms_FYM_poultry = c(0.2, 3.2),      # assumed solid storage
    awms_anaerobic = c(0.0, 3.2)         # from IPCC (2019)
  ) %>%
    # Transpose the table 
    pivot_longer(-factor, names_to = "management_system", values_to = "value") %>%
    pivot_wider(names_from = factor, values_from = value)
  head(MM_data)
  
  # load in the GE required per day
  dailyGE <- fread(file.path(savePath, "mjDayGE.csv")) %>%
    as.data.frame()
  head(dailyGE)
  # info from Brown et al. (2021) Table A 3.3.5
  ## get all unique categories first
  cats <- unique(dailyGE$category)
  Stor_frac <- data.frame(category = cats
                          , awms_slurry = NA
                          , awms_daily = NA
                          , awms_past = NA 
                          , awms_anaerobic = NA
                          , awms_FYM_poultry = NA
                          , awms_FYM_cows_pigs = NA
                          , awms_FYM_sheep = NA) %>%
    t() %>% as.data.frame()
  names(Stor_frac) <- Stor_frac[1,]
  Stor_frac <- Stor_frac[-1,]
  head(Stor_frac)
  Stor_frac2 <- Stor_frac %>%
    mutate(ewes = c(0, 0, 92, 0, 0, 0, 8)
           , pigs = c(36, 14, 10, 4, 0, 36, 0)
           , other_cattle = c(18, 12, 48, 1, 0, 21, 0)
           , calf = c(18, 12, 48, 1, 0, 21, 0)
           , heifer = c(18, 12, 48, 1, 0, 21, 0)
           , lactating_dairy = c(60, 8, 22, 2, 0, 9, 0)
           , other_sheep = c(0, 0, 99, 0, 0, 0, 1)
           , rams = c(0, 0, 99, 0, 0, 0, 1)
           , poultry = c(0, 39, 3, 8, 50, 0, 0)
    ) %>%
    mutate(management_system = rownames(.))
  head(Stor_frac2)

  # merge, to get all info
  MMStorage <- Stor_frac2 %>%
    pivot_longer(!c(management_system)
                 , names_to = "category", values_to = "awms_pc") %>%
    merge(., MM_data)
  head(MMStorage)
  lStor_frac2 <- MMStorage %>%
    filter(category == "lactating_dairy")
  
  ###### 5a - volatile solids (VS) ######
  ## ---------- notes ---------- # 
  ## Default VS value from IPCC (2019) is 8.4 (1000 kg animal mass-1) day-1 
  
  ## If average daily VS excretion rates are not available, country-specific VS
  ## excretion rates can be estimated from feed intake levels (IPCC, 2019). These
  ## were calculated in the ghg_excretions.
  
  ## Here, values from the most recent UK inventory is used (Ward et al., 2024) 
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
  
  # load in criteria match
  load("critMatch.RData")
  
  # get means for all of the same animal, category
  summariseVS <- VStable %>%
    group_by(animal, category) %>%
    dplyr::summarise(across(starts_with("kgVS"), list(mean = ~ mean(., na.rm = TRUE)))) %>%
    filter(grepl("lactat", category)) %>%
    # refine to just the criteria matched ones
    dplyr::select(1:2, contains(c(gsub("MJdayGE", "", criteria.match))))
  head(summariseVS)
  
  vss <- VStable %>%
    filter(grepl("lactat", category)) %>%
    # refine to just the criteria matched ones
    dplyr::select(1:2
                  , contains(c(gsub("MJdayGE", "", criteria.match)))
    )
  
  ## ------------ Notes --------------  ##
  ## from Appuhamy et al. (2018), the mean VS (of lactating dairy cows) was 5.73 kg/d,
  ## with a SD of 1.79, and a min and max of 1.67 and 10.86, respectively. All of the
  ## the values selected fit in between the 3.94 - 7.52 (mean +/- SD) range. Those that did not were removed
  ## (see vss.criteria)
  ## ------------ ----- --------------  ##
  
  UpperLim <- 7.52
  lowerLim <- 3.94
  
  # stop("before VS trial")
  ## ------------ Notes --------------  ##
  ## get the average of the selected DE compared published values
  ## you can comment out this section to test other values
  # DEinside <- 85 # 'housed' value
  # DEoutside <- 79 # outside value
  # inDE <- paste0("housed_", DEinside)
  # outDE <- paste0("sml_", DEoutside)
  # 
  # limitTest <- vss %>%
  #   dplyr::select(animal, category
  #                 , contains(c(inDE, outDE)))
  # range(limitTest[, 3]); mean(limitTest[, 3])
  # range(limitTest[, 4]); mean(limitTest[, 4])
  # mean(c(mean(limitTest[, 3]), mean(limitTest[, 4])))
  # mean(c(mean(limitTest[, 3]), mean(limitTest[, 4]))) - 5.73
  ## ------------ ----- --------------  ##
  
  # see which are too high or low
  VScheck <- VStable %>%
    filter(grepl("lactat", category)) %>%
    # group by animal and category, to get average regardless of region
    group_by(animal, category) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    mutate(across(where(is.numeric), ~ ifelse(. > lowerLim & . < UpperLim, 1, 0)))
  
  # check the final amount
  colSums(VScheck[, 3:ncol(VScheck)])
  table(colSums(VScheck[, 3:ncol(VScheck)]))
  table(as.matrix(VScheck[, 3:ncol(VScheck)]))
  which(colSums(VScheck[, 3:ncol(VScheck)]) == 24)
  cn <- names(VScheck %>% ungroup() %>%
                dplyr::select(kgVSdayhoused_65:kgVSdaylrg_85))[which(colSums(VScheck[, 3:ncol(VScheck)]) == 24)]
  cn
  # reduce criteria match to ones that match both the literature stated EF and VS quantities
  criteria.match2 <- criteria.match[which(grepl(paste(gsub("kgVSday", "", as.character(cn)), collapse =  "|")
                                                , gsub("MJdayGE", "", criteria.match)))]
  criteria.match2 <- gsub("MJdayGE", "", criteria.match2)
  criteria.match2
  
  UpperMax <- 10.86
  lowerMin <- 1.67
  
  # see which are too high or low - based on max and min
  VScheckMM <- VStable %>%
    filter(grepl("lactat", category)) %>%
    # group by animal and category, to get average regardless of region
    group_by(animal, category) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    mutate(across(where(is.numeric), ~ ifelse(. > lowerMin & . < UpperMax, 1, 0)))
  
  ###### 5b - per-animal volatile solids (VS) ######
  ## ------------ Notes --------------  ##
  ## this final output will be based upon the scenario of DE written in the 
  ## '0 - select animal scenario' section above
  ## ------------ ----- --------------  ##
  
  # stop("before VS summary")
  # set names to keep
  insideName <- paste0("kgVSdayhoused_", DEinside)
  outsideName <- paste0("kgVSdaysml_", DEoutside)
  cat("inside value =", insideName, "|", "outside value =", outsideName, "\n")
  
  VStableScen <- VStable %>%
    dplyr::select(animal, category, UK.Region
                  , (all_of(c(insideName)))
                  , (all_of(c(outsideName)))) %>%
    rename("kgVSdayIn" = insideName
           , "kgVSdayOut" = outsideName) %>%
    # get the 12-month value of least, most realistic, DE
    mutate(kgVSyr = kgVSdayOut * (365)) %>%
    dplyr::select(-c(kgVSdayIn, kgVSdayOut))
  head(VStableScen)
  
  ## for dairy cattle, change VS to 5.3 VS day-1 animal-1, to match most recent inventory
  VStableScen <- VStableScen %>%
    mutate(kgVSyr = ifelse(category == "lactating_dairy", 5.3 * 365, kgVSyr))
  head(VStableScen)
  
  ##### Tier 1 calculation #####
  ## ---------- notes ---------- # 
  ## The IPCC Tier 1 method relies on default emission values based on 
  ## average temperature of a region, and the livestock 'species'
  ## it also  requires VS, and fraction of manure management type
  ## --------------------------- # 
  Tier1Manure <- VStableScen %>%
    # get yearly VS, in kg
    dplyr::select(animal, category, UK.Region, kgVSyr) %>%
    merge(MMStorage)
  head(Tier1Manure)
  
  # VS * AMWS (fraction) * EF
  Tier1Manure_perAnimal <- Tier1Manure %>%
    mutate(ch4kg_yr = (kgVSyr/365) * (awms_pc/100) * EF) %>%
    # sum all types of management
    group_by(category, animal, UK.Region) %>%
    summarise(ch4kg_yr = sum(ch4kg_yr))
  head(Tier1Manure_perAnimal)
  
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
  head(VStableScen)
  
  ##### Tier 2 calculation #####
  # stop("before tier 2")
  ## ---------- notes ---------- # 
  ## The IPCC Tier 2 [using VS of 5.3 for lactating cows]
  ## --------------------------- # 
  head(Tier1Manure)
  
  # AMWS (fraction) * MCF (fraction)
  Tier2Manure_perAnimal <- Tier1Manure %>%
    mutate(mcfXawms = (awms_pc/100) * (MCF/100)) %>%
  # sum all types of management
    group_by(category, animal, UK.Region) %>%
    summarise(sumAwmsMcf = sum(mcfXawms)) %>%

    ## ---------- notes ---------- # 
    ## This next part uses the equation 10.23 in IPCC (2019)
    ## It assumes that the different management systems correspond to Table A 3.3.5 in Brown et al. (2021)
    ## --------------------------- # 
    
    # multiply by b0 and conversion (0.67) %>%
    merge(., VStableScen) %>%
    mutate(b0conv = sumAwmsMcf * 0.67 * B0) %>%
    # multiply by VS (yearly)
    mutate(ch4kg_yr = b0conv * kgVSyr)
  head(Tier2Manure_perAnimal)
  
  ## ------------ Notes --------------  ##
  ## testing bit
  ccc <- Tier1Manure %>%
    mutate(mcfXawms = (awms_pc/100) * (MCF/100)) %>%
    # sum all types of management
    group_by(category, animal, UK.Region) %>%
    filter(category == "lactating_dairy") %>%
    filter(management_system %in% c("awms_anaerobic", "awms_past", "awms_FYM_cows_pigs"
                                    , "awms_daily", "awms_slurry"))
  ## ------------ ----- --------------  ##
  
  ## ---------- notes ---------- # 
  ## now the VS table is used to multiply the values obtained above
  ## first daily CH4 from manure management will be obtained, with the amount  
  ## then be multiplied by 365 to get yearly values
  ## --------------------------- # 
  
  # get means for all of the same animal, category
  summariseCH <- Tier2Manure_perAnimal %>%
    group_by(animal, category) %>%
    dplyr::summarise(across(starts_with("ch4kg_yr"), list(mean = ~ mean(., na.rm = TRUE)
                                                         , min = ~ min(., na.rm = TRUE)
                                                         , max = ~ max(., na.rm = TRUE))))
  head(summariseCH)
  
  ###### 5f - convert CH4 to CO2e ######
  # convert to CO2e (non-fossil fuel: 28)
  MMCO2e <- Tier2Manure_perAnimal %>%
    mutate(CO2e_kg_yr = ch4kg_yr * gwpCH4)
  head(MMCO2e)
  
  # save tables
  # CH4 and CO2
  fwrite(MMCO2e, file.path(savePath, "MMCO2e.csv"), row.names = F)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(savePath, "readme_MMCO2e.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-07-13
  Last update:",  format(Sys.Date()), "
  
  The files show the calculated amount of emissions from manure management per animals, per six months, to do with the animals' situations i.e. the quality of their diet, and their freedom of moved (e.g. housed, or allowed to roam)  
  
  Files:
    MMCO2e.csv
  
  Columns of 'MMCO2e' file:
      'animal' = the category an animal was assigned to, based on either CTS or AgCensus data
      'category' = the category each animal type was assigned to in this work
      'UK.Region' = different regions of the UK, with 'other' representing non-English regions
      'CO2e_kg_yr' = annual emissions produced of CO2e (carbon dioxide equivalent) from manure management, in kg, per-animal 
      'ch4kg_yr' = annual emissions produced of CH4 (methane) from manure management, in kg, per-animal 
          
  Data: per-animal  production
  units: kg manure management emissions per animal")
  sink(file = NULL)
  
  # tidy
  rm(part5ManureManagement)
  
} # end of 'part5ManureManagement'

#### 6 - part6excretions ####

if(part6excretions){
  
  # stop("part 6")
  
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
                                 # , otherNintakeDay(GEpivot$DMIkgDay, halfInOut/10)
                                 , otherNintakeDay(GEpivot$DMIkgDay, 16.1)
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
  apply(Nintake.retain.excreteCows, 2, min)
  
  # stop("trial N excretion")
  #### trial DE for N excretion ####
  ## ------------ Notes --------------  ##
  ## you can comment out this section to test other values
  # DEinside <- 83 # 'housed' value
  # DEoutside <- 79 # outside value
  ## ------------ ----- --------------  ##
  
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
  
  # should be near: 130.8 kg N/ dairy cow / yr [based on 358.4 g/d (Bougouin et al., 2022)]
  # which is equivalent to 65.4 for the 6-month portion outside
  ## get means for all of the same animal, category
  summariseN0 <- NexcreteCows %>%
    group_by(animal, category) %>%
    dplyr::summarise(across(where(is.numeric)
                            , list(mean = ~ mean(., na.rm = TRUE)))) %>%
    filter(grepl("dairy", category)) %>%
    dplyr::select(animal, category, contains("ExcreteYr"))
  head(summariseN0)
  
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
  table(excreteNcombined$UK.Region)
  
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
  rm(summariseN, csgNintakeDay, NexcreteCows
     , NexcreteNoBovine, Nintake.retainNoBovine, Nintake.retain.excreteNoBovine)
  
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
    # times all by 265 to convert to CO2e
    mutate(grazEmisYr_kgCO2e = grazEmisYr_kgN2O * gwpN2O)
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
    'grazEmisYr_kgCO2e' = total (direct + indirect) N emissions, in kg CO2e, from grazing. A GWP of 265 was used to convert from N2O to CO2e 

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
  
  # stop("before animal grid loading in part 6")
  
  # load the animal grid, which is at the 1 km2
  animalGrid <- st_read(file.path(currPath, "all_animals_1km.gpkg"))
  head(animalGrid)
  
  ## ------------ Notes --------------  ##
  ## The below section tests two methods, the first just confirming that the
  ## second (which is faster) is accurate
  ## ------------ ----- --------------  ##
  
  ###### first (test) long multiplication method ######
  # make long, after combing to get region
  reg <- st_read("data_in/animals/nonBovine_grid_2020_1km.gpkg") %>%
    dplyr::select(rcFid_1km, region) %>% st_drop_geometry()
  animalGrid <- animalGrid %>%
    merge(., reg, all = F)
  
  animalGridLong <- animalGrid %>%
    st_drop_geometry() %>%
    pivot_longer(!c(rcFid_1km, region))
  head(animalGridLong)
  # stop("animal grid")
  
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
  ukRegionsPossible <- na.omit(ukRegionsPossible)
  
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
      filter(region == x)
    
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
    if(identical(
      sort(c(names(regionGrid)[3:(ncol(regionGrid))]))
      , sort(c(regionN$animal))
    )){
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
  
  ## check the example and the actual match
  x1 <- sumExcreteN %>%
    filter(rcFid_1km == "487500_215500")
  x2 <- spatialGridN %>%
    filter(rcFid_1km == "487500_215500")
  head(x1)
  head(x2)
  stopifnot(x1$total1km_kgN == x2$totalN)
  
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
  
  ##### 6d - SCENARIOS - calculate the amount of N in manure #####
  if(nManureScens){
    cat(scenarios, sep = "\n")
    # stop("before N scenarios")
    
    # create scenario save path
    scenarioResultsPath <- file.path("scenario", "results")
    
    # load in all the animal grids for the scenarios
    for(i in 1:length(scenarios)){
      tic("scenario run")
      scenName <- gsub(".gpkg", "", basename(scenarios[[i]]))
      cat(paste0("starting N excretion for ... ", scenName, "[", i, "] ..."), "\n")
      # load in full dataset
      fulldf.in <- st_read(scenarios[[i]])
      head(fulldf.in)
      # use original regions (i.e., not county)
      head(reg)
      cat(unique(reg$region), sep = "\n")
      
      ## combine them
      fulldf.in <- fulldf.in %>%
        merge(., reg, by = "rcFid_1km")
      sort(names(fulldf.in))
      
      ### make long
      fulldfLong <- fulldf.in %>%
        # only keep the animals
        dplyr::select(-any_of(c("Name", "lcm", "n"
                                , contains("_ha")))) %>%
        st_drop_geometry() %>%
        pivot_longer(!c(rcFid_1km, region))
      head(fulldfLong)
      
      tic("merging excrete N")
      # merge with the excreted N df, ensuring regions match
      fulldf.ExcreteN <- fulldfLong %>% as.data.frame() %>%
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
      head(fulldf.ExcreteN)
      toc()
      
      table(fulldf.ExcreteN$region)
      table(fulldf.ExcreteN$UK.Region)
      
      # multiply N excreted by value (number of animals in that km2)
      fulldf.ExcreteN2 <- fulldf.ExcreteN %>%
        mutate(totalN = value * kgNExcreteYr_housed)
      head(fulldf.ExcreteN2)
      # sum by pixel
      fulldf.sumExcreteN <- fulldf.ExcreteN2 %>%
        group_by(rcFid_1km) %>%
        summarise(total1km_kgN = sum(totalN, na.rm = T))
      head(fulldf.sumExcreteN)
      
      ###### second matrix multiplication method - SCENARIOS ######
      cat(paste0("starting [exercise 2] for N excretion for ... ", scenName, " ..."), "\n")
      
      # only keep animals
      fulldf.onlyAnimals <- fulldf.in %>%
        # only keep the animals
        dplyr::select(-any_of(c("Name", "lcm", "n", "lcm2007"))) %>%
        dplyr::select(-contains("_ha"))
      
      # add underscores to regions, to make them match output df
      fulldf.onlyAnimals$region <- gsub(" ", "_", fulldf.onlyAnimals$region)
      head(fulldf.onlyAnimals$region)
      # and make London the same as South East
      fulldf.onlyAnimals$region <- ifelse(fulldf.onlyAnimals$region == "London", "South_East"
                                          , fulldf.onlyAnimals$region)
      # convert the humber, to match excrete N df
      fulldf.onlyAnimals$region <- ifelse(fulldf.onlyAnimals$region == "Yorkshire_and_The_Humber"
                                          , "Yorkshire_and_the_Humber"
                                          , fulldf.onlyAnimals$region)
      unique(excreteNcombined$UK.Region)
      unique(fulldf.onlyAnimals$region)
      ## get all different possible regions
      ukRegionsPossible <- unique(fulldf.onlyAnimals$region)
      ukRegionsPossible <- na.omit(ukRegionsPossible)
      
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
        regionGrid <- fulldf.onlyAnimals %>% st_drop_geometry() %>%
          filter(region == x)
        head(regionGrid)
        
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
        if(identical(
          sort(c(names(regionGrid)[3:(ncol(regionGrid))]))
          , sort(c(regionN$animal))
        )){
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
      
      ## check the example and the actual match
      x1 <- fulldf.sumExcreteN %>%
        filter(rcFid_1km == "487500_215500")
      x2 <- spatialGridN %>%
        filter(rcFid_1km == "487500_215500")
      head(x1)
      head(x2)
      stopifnot(x1$total1km_kgN == x2$totalN)
      
      ### save
      st_write(spatialGridN, file.path(scenarioResultsPath
                                       , paste0("spatialN1km_kgN", scenName, ".gpkg"))
               , append = F)
    }
    toc(log = T)
  }
} # end 'part6excretions'

# tidy
rm(DEinside, DEoutside)

#### 7 - part7fertiliserUse ####
if(part7fertiliserUse){
  
  # stop("fertiliser time")
  
  ## ------------ Notes --------------  ##
  ## Fertiliser use emissions are split into lime application and 
  ## other fertiliser application (namely ammonia nitrate and urea)
  
  ## Both of the above factors will be impacted by the type of crop
  ## (and land cover generally) in a pixel.
  
  ## It is assumed that the underlying soils remain the same whether modelling
  ## 2015 data or scenario data. Therefore, 7A deals with 
  ## generally soil characteristics, while later sections assess each data range
  ## in turn (i.e., 2015, and then scenario data)
  ## ------------ ----- --------------  ##
  
  cat("----------- Starting fertiliser use section -----------\n\n")
  
  ##### 7a - Determine fert-relevant soil characteristics #####
  ## ------------ Notes --------------  ##
  ## pH and other soil data come from the P drive, which contains data from 
  ## Cranfield
  ## ------------ ----- --------------  ##
  
  # list the soil data
  # should inlcude: pH, silt, sand, and clay content, and OC.
  soils <- list.files(soilsData
                      , pattern = ".tif$", full.names = T)
  # cat(basename(soils), sep = "\n")
  # combine to a df
  soils <- rast(lapply(soils, rast))
  # print(soils)
  soils <- as.data.frame(soils, xy = T)
  soils <- soils[complete.cases(soils$Soil_TopsoilClayContent_1k), ]
  head(soils)
  # save
  fwrite(soils, "data_in/crops/soil_data.csv", row.names = F)
  
  # read as point data
  soils <- st_as_sf(soils, coords = c("x", "y"), crs = crs(27700))
  st_crs(soils) <- 27700
  
  # calculate possible lime added per hectare based on Table 1.2 in RB209
  # this makes the assumption that farmers are going for optimal yield
  
  # determine what type of soil it is: sandy, silty, clayey
  # base this on the highest proportion in the soil's composition
  soils$soilTyp <- ifelse(soils$Soil_TopsoilSandContent_1k > soils$Soil_TopsoilSiltContent_1k 
                          & soils$Soil_TopsoilSandContent_1k > soils$Soil_TopsoilClayContent_1k
                          , "Sand"
                          , ifelse(soils$Soil_TopsoilSiltContent_1k > soils$Soil_TopsoilSandContent_1k 
                                   & soils$Soil_TopsoilSiltContent_1k > soils$Soil_TopsoilClayContent_1k
                                   , "Silt", "Clay"))
  head(soils)
  
  # save
  cat("*Saving soils here:", file.path("data_in", "crops", "soils.gpkg"), "*\n")
  st_write(soils, file.path("data_in", "crops", "soils.gpkg"), append = F, quiet = T)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path("data_in", "crops", "readme_soils.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-06-05
  Last update:",  format(Sys.Date()), "
  
  Description of 'soils.gpkg':
  Contains soil information at 1 km
  
  Columns of 'soils.gpkg':
    'Soil_Topsoil[type]Content_1k' = percentage of a soil type in the 1 km2
      'type' = either 'sand', 'silt', or 'clay' 
    'Soil_TopsoilKFactor_1k' = soil erodibility factor
    'Soil_TopsoilOrganicCarbonContent_1k' = percentage of carbon content in soil
    'Soil_TopsoilPH_1k' = pH of soil
    'soilTyp' = the type of soil each km is considered to be based on which of the substrates was highest in the clay, silt, or sand columns
    
  Data: soil data
  units: see above
  Spatial resolution: 1 km2
  Spatial extent: UK
  Temporal coverage: 2015
  Native projection: 27700")
  sink(file = NULL)
  
  ## ------------ Notes --------------  ##
  ## From RB209: 'For each field, the amount of lime to apply will depend on the current soil pH,
  ## soil texture, soil organic matter and the target pH, which should be 0.2 pH points above optimum'
  ## Therefore, Table 1.2 in RB209 was used to determine how much lime would be applied
  ## ------------ ----- --------------  ##
  
  # stop("before lime requirements")
  
  # calculate lime requirements, in t/ha, for arable and arable
  ## for sand
  limeReqsSand <- soils %>%
    filter(soilTyp == "Sand") %>%
    # calculate based on starting pH
    mutate(limeReq_tha_arab = if_else(Soil_TopsoilPH_1k <= 5, 10
                                      , if_else(Soil_TopsoilPH_1k <= 5.5, 7
                                                , if_else(Soil_TopsoilPH_1k <= 6.0, 4
                                                          , if_else(Soil_TopsoilPH_1k <= 6.2, 3, 0))))
           , limeReq_tha_gras = if_else(Soil_TopsoilPH_1k <= 5, 5
                                        , if_else(Soil_TopsoilPH_1k <= 5.5, 3
                                                  , if_else(Soil_TopsoilPH_1k <= 6.0, 0
                                                            , if_else(Soil_TopsoilPH_1k <= 6.2, 0, 0)))))
  
  # for silt
  limeReqsSilt <- soils %>%
    filter(soilTyp == "Silt") %>%
    # calculate based on starting pH
    mutate(limeReq_tha_arab = if_else(Soil_TopsoilPH_1k <= 5, 12
                                      , if_else(Soil_TopsoilPH_1k <= 5.5, 8
                                                , if_else(Soil_TopsoilPH_1k <= 6.0, 5
                                                          , if_else(Soil_TopsoilPH_1k <= 6.2, 4, 0))))
           , limeReq_tha_gras = if_else(Soil_TopsoilPH_1k <= 5, 6
                                        , if_else(Soil_TopsoilPH_1k <= 5.5, 4
                                                  , if_else(Soil_TopsoilPH_1k <= 6.0, 0
                                                            , if_else(Soil_TopsoilPH_1k <= 6.2, 0, 0)))))
  # for clay
  limeReqsClay <- soils %>%
    filter(soilTyp == "Clay") %>%
    # calculate based on starting pH
    mutate(limeReq_tha_arab = if_else(Soil_TopsoilPH_1k <= 5, 14
                                      , if_else(Soil_TopsoilPH_1k <= 5.5, 10
                                                , if_else(Soil_TopsoilPH_1k <= 6.0, 6
                                                          , if_else(Soil_TopsoilPH_1k <= 6.2, 4, 0))))
           , limeReq_tha_gras = if_else(Soil_TopsoilPH_1k <= 5, 7
                                        , if_else(Soil_TopsoilPH_1k <= 5.5, 4
                                                  , if_else(Soil_TopsoilPH_1k <= 6.0, 0
                                                            , if_else(Soil_TopsoilPH_1k <= 6.2, 0, 0)))))
  
  # combine - this shows tonnes of lime required per ha based on soil pH
  limeReqs <- rbind(limeReqsSand, limeReqsSilt, limeReqsClay)
  head(limeReqs)
  # tidy
  rm(limeReqsSand, limeReqsSilt, limeReqsClay)
  
  # stop("after lime requirements")
  
  ##### 7b - Calculate factors that will affect fertiliser emissions [rainfall, Soil Nitrogen Supply (SNS)] #####
  ## ------------ Notes --------------  ##
  ## According to RB209, the amount of fertiliser used is dependent on the rainfall amount, 
  ## the previous crop type, and the Soil Nitrogen Supply (SNS)
  ## The steps it recommends are:
  ##  A. Identify the soil category for the field
  ##  B. Identify the previous crop
  ##  C. Select the rainfall range for the field
  ##  D. Identify the provisional SNS Index using the appropriate table
  ## These steps will be followed to create a table to indicate the fertiliser to be used
  ## ------------ ----- --------------  ##
  
  ##### 7b1. Identify the soil category for the field #####
  ## ------------ Notes --------------  ##
  ## The soil type was determined in the 'soils' df above, giving the soil type as 'sand', 'silt', or 'clay'
  ## ------------ ----- --------------  ##
  
  ##### 7b2. Identify the previous crop #####
  ## ------------ Notes --------------  ##
  ## It is assumed that the previous crop will stay the same as the current one
  ## which will depend on whether 2015 or scenario data is used. This is done in
  ## the next section
  ## ------------ ----- --------------  ##
  
  ##### 7b3. Select the rainfall range for the field #####
  ## ------------ Notes --------------  ##
  ## RB209 state "‘low rainfall’ (annual rainfall less than 600 mm ...),
  ## ‘moderate rainfall’ (between 600–700 mm annual rainfall ...), and
  ## ‘high rainfall’ (over 700 mm annual rainfall ...).
  
  # create a rainfall grid from 2012 - 2017 precip data. Use the 'all_mean' data
  rainExtract <- raster("N:/Data/UK/chessmet/precipitation/year_sum_averageikm.tif") %>%
    # convert to df
    as.data.frame(xy=T) %>%
    # using the info on Step 3 from RB209 (pg 9), put rain into three categories
    mutate(rainCat = if_else(layer < 600, "low"
                             , if_else(layer <= 700, "moderate"
                                       , "high")))
  
  
  ##### 7b4. Identify the provisional SNS Index using the appropriate table #####
  ## ------------ Notes --------------  ##
  ## To determine SNS, soil, rainfall, and previous crop need to be combined to use Tables 4.2 - 4.4 in RB209
  ## also need depth - derived using root depth
  ## ------------ ----- --------------  ##
  
  # separate centroid points of the soil data
  separatedSoils <- soils %>%
    mutate(X = as.integer(unlist(map(soils$geometry,1))),
           Y = as.integer(unlist(map(soils$geometry,2)))) %>%
    dplyr::select(X, Y, soilTyp) %>%
    # merge with rain
    merge(., rainExtract %>% 
            mutate(X = as.integer(x), Y = as.integer(y)) %>%
            dplyr::select(X, Y, rainCat))
  
  ## require root depth - data converted from EU scale
  rootDepthUK <- rast(file.path("data_in/crops", "rootdepthuk.tif")) %>%
    # convert to df
    as.data.frame(xy=T) %>%
    # merge with rainfall data
    merge(separatedSoils, .
          , by.x = c("X", "Y")
          , by.y = c("x", "y")
          , all = F) %>%
    # determine RB209 rooting depth category
    mutate(depthRB209 = if_else(soilTyp %in% c("Sand", "Silt") & rootdepthuk > 100, "Deep"
                                , if_else(soilTyp == "Clay" & rootdepthuk > 40, "Deep" 
                                          , if_else(rootdepthuk < 40, "Shallow", "Medium"))))
  head(rootDepthUK)
  
  ## ------------ Notes --------------  ##
  ## The purpose of the next bit of code is to derive the exact values for SNS
  ## from the literature. Two tables will be created 'empty' and will then be
  ## filled in manually and loaded back in
  ## ------------ ----- --------------  ##
  
  ## this code either reads in or creates the initial table ##
  if(file.exists(file.path("./data_in/crops", "snsTableCompleted.csv"))){
    # read in and then replicate all cereals for the cereal types
    snsParameters <- fread(file.path("./data_in/crops", "snsTableCompleted.csv"))
    snsCereals <- snsParameters %>%
      filter(crop_type == "cereal")
    # get nrow
    nrowCereal <- nrow(snsCereals)
    snsCereals <- do.call("rbind", replicate(5, snsCereals, simplify = FALSE))
    snsCereals$crop <- c(rep(c("maize"), nrowCereal)
                         , rep(c("winterwheat"), nrowCereal)
                         , rep(c("winterbarley"), nrowCereal)
                         , rep(c("springwheat"), nrowCereal)
                         , rep(c("springbarley"), nrowCereal)) 
    snsCereals <- snsCereals %>%
      dplyr::select(-"crop_type") %>%
      rename(crop_type = crop) %>%
      relocate(threeParas, crop_type, sns_category)
    # add back to original
    snsParameters <- rbind(snsParameters %>%
                             filter(crop_type != "cereal")
                           , snsCereals)
  } else {
    
    # merge with crop data
    cropSoilRain <- rootDepthUK %>%
      merge(., iscen %>% st_drop_geometry() %>%
              dplyr::select(rcFid_1km, X, Y, scen
                            , ends_with("_ha"))
            , by = c("X", "Y"), all = T) %>%
      filter(!is.na(rcFid_1km)) %>%
      filter(!is.na(soilTyp))
    
    # create a table from the unique crop and categories
    # the three parameters define what SNS category the land is
    ## unique RB209 categories
    uniqueRb209Crops <- c("beans"
                          , "cereal" # cereals are maize, winterwheat, winterbarley, springwheat, springbarley
                          , "forage" # forage crops are oats 
                          , "osr", "potatoes", "sugarbeet"
                          , "uncropped") # uncropped is 'improved grassland' 
    snsParameters <- cbind(threeParas = rep(unique(cropSoilRain$snsCats), length(uniqueRb209Crops))
                           , crop_type = uniqueRb209Crops) %>% as.data.frame()
    fwrite(snsParameters, file.path("data_in/crops", "snsTableEmpty.csv"), row.names = F)
  }
  head(snsParameters)
  
  ## this code either reads in or creates the initial sns crop table ##
  if(file.exists(file.path("data_in/crops", "snsCropTableCompleted.csv"))){
    snsCropParameters <- fread(file.path("data_in/crops", "snsCropTableCompleted.csv"))
  } else {
    
    # merge with crop data
    cropSoilRain <- rootDepthUK %>%
      merge(., iscen %>% st_drop_geometry() %>%
              dplyr::select(rcFid_1km, X, Y, scen
                            , ends_with("_ha"))
            , by = c("X", "Y"), all = T) %>%
      filter(!is.na(rcFid_1km)) %>%
      filter(!is.na(soilTyp))
    # create a table from the unique crop and categories
    ## unique RB209 categories
    uniqueRb209Crops <- c("beans", "maize", "winterwheat", "winterbarley", "springwheat", "springbarley"
                          , "forage" # forage crops are oats 
                          , "osr", "potatoes", "sugarbeet"
                          , "uncropped") # uncropped is lc_401 'improved grassland' 
    snsCropParameters <- cbind(threeParas = rep(unique(cropSoilRain$snsCats), length(uniqueRb209Crops))
                               , crop_type = uniqueRb209Crops
                               , sns0 = 0
                               , sns1 = 1
                               , sns2 = 2
                               , sns3 = 3) %>% as.data.frame()
    fwrite(snsCropParameters, file.path("data_in/crops", "snsCropTableEmpty.csv"), row.names = F)
  }
  head(snsCropParameters)
  
  # make tables wide. These tables contain the three parameters that define what SNS category the land is
  snsParasWide <- snsParameters %>%
    tidyr::pivot_wider(names_from = crop_type, values_from = sns_category)
  colnames(snsParasWide)[2:ncol(snsParasWide)] <- paste0("sns_", colnames(snsParasWide)[2:ncol(snsParasWide)])
  head(snsParasWide)
  # sns parameters
  snsCropParametersWide <- snsCropParameters %>%
    dplyr::select(-ref) %>%
    tidyr::pivot_wider(names_from = crop_type, values_from = c(sns0, sns1, sns2, sns3)) %>%
    rename_at(vars(-(1)), ~ paste0(., '_kgHa'))
  head(snsCropParametersWide)
  names(snsCropParametersWide)
  
  # stop("before 7a - 4575")
  
  #### 7c - determine emissions factors for fertiliser use ####
  ## ------------ Notes --------------  ##
  ## To calculate emissions from fertilisers, several bits of information are required:
  ## soil texture
  ## Drainage
  ## pH
  ## soil  moisture
  ## CEC
  ## SOC
  ## application method
  ## fertiliser type
  ## application rate
  ## crop type
  
  ## these need to be incorporated and derived per cell. With the exception of
  ## crop type, all can be calculated based on soil properties.
  ## Crop type will be added in the scenario modelling section
  ## ------------ ----- --------------  ##
  
  ## add a df here that connects rcFid_1km and X and Y
  ## ------------ Notes --------------  ##
  ## use cropArea which has land cover as part of it
  ## get centroids, then get X and Y
  ## ------------ ----- --------------  ##
  XYrcFid_1km <- ghgGridGB %>% st_centroid() %>%
    mutate(X = as.integer(unlist(map(.$geom, 1))),
           Y = as.integer(unlist(map(.$geom, 2))))
  head(XYrcFid_1km)
  # merge with soils
  XYrcFidSoils <- XYrcFid_1km %>%
    merge(soils %>%
            mutate(X = as.integer(unlist(map(soils$geom, 1))),
                   Y = as.integer(unlist(map(soils$geom, 2)))) %>%
            st_drop_geometry())
  rm(XYrcFid_1km)
  
  ##### 7c1 - soil texture #####
  ## ------------ Notes --------------  ##
  ## there are three categories for soil texture: 'fine', 'medium', or 'course'
  ## this will be based on the soil type that that cell was assigned to
  ## ------------ ----- --------------  ##
  fertiliserFactors <- XYrcFidSoils %>%
    # add soil texture 
    # create list
    mutate(textureCat = if_else(soilTyp == "Sand", "fine"
                                , if_else(soilTyp == "Silt", "medium", "course")))
  # match with the coefficients
  textureCoefs <- c(fine = c(0,0,0)
                    , medium = c(-0.472, 0, 0)
                    , course = c(-0.008, 0, 0))
  
  ##### 7c2 - drainage #####
  ## ------------ Notes --------------  ##
  ## there are two categories, but no data, therefore drainage was presumed to be
  # in between both, which gave values of -0.21, 0.473, 0
  ## ------------ ----- --------------  ##
  drainageCoefs <- c(-0.21, 0.473, 0)
  head(fertiliserFactors)
  
  ##### 7c3 - pH #####
  ## ------------ Notes --------------  ##
  ## pH can be derived from the soil dataframe
  ## ------------ ----- --------------  ##
  fertiliserFactors <- fertiliserFactors %>%
    dplyr::select(rcFid_1km, textureCat, Soil_TopsoilPH_1k
                  , Soil_TopsoilClayContent_1k, Soil_TopsoilOrganicCarbonContent_1k) %>%
    # add soil texture 
    # create list
    mutate(phCat = if_else(Soil_TopsoilPH_1k <= 5.5, "low pH"
                           , if_else(Soil_TopsoilPH_1k <= 7.3, "low-mod pH"
                                     , if_else(Soil_TopsoilPH_1k <= 8.5, "high-mod pH", "high pH"))))
  head(fertiliserFactors)
  
  # get the coefficients from Table 2.4 in CFT
  phCoefs <- c(low_ph = c(0,0,-1.072)
               , low_modph = c(0.109,0,-0.933)
               , high_modph = c(-0.352,0,-0.608)
               , highph = c(-0.352, 0, 0))
  
  ##### 7c4 - soil moisture #####
  ## ------------ Notes --------------  ##
  ## there are data for soil moisture - from Copernicus
  ## ------------ ----- --------------  ##
  # load soil moisture data
  sm <- raster("N:/Data/UK/soil/soil_moisture/yearmean15_20.tif")
  sm <- as.data.frame(sm, xy = T)
  # there are two categories for soil moisture in CFT: moist and dry
  # assuming that half of england fits both these categories, the mean will decide which is which
  moistDryMean <- median(sm$yearmean15_20, na.rm = T)
  moistDryCat <- sm %>% 
    mutate(moistDryCat = if_else(yearmean15_20 <= moistDryMean, "dry", "moist")) %>%
    select(yearmean15_20, moistDryCat, x, y) %>%
    mutate(rcFid_1km = paste(x, y, sep = "_"))
  
  fertiliserFactors <- fertiliserFactors %>%
    merge(., moistDryCat, by = "rcFid_1km"
          , all = F)
  head(fertiliserFactors)
  
  ##### 7c4 - CEC #####
  ## ------------ Notes --------------  ##
  ## Cation-exchange capacity (CEC) is a measure of how many cations can be retained on soil particle surfaces
  ## this cannot be changed. It corresponds to ## the total capacity of a material to hold exchangeable cations
  ## measured in centimoles of charge per kilogram of exchanger (cmol(+)/kg).
  ## according to https://www.soilquality.org.au/factsheets/cation-exchange-capacity,
  ## sand has low cec (~2), while organic matter can range from 250 - 400 
  ## requires: pH, bulk density, SOC, and clay content of soil
  ## ------------ ----- --------------  ##
  
  # bulk density not available, so all soil was given a medium rating of 1.3 g g-3
  fertiliserFactors$BD <- 1.3
  
  # convert soil clay to 1 of 3 categories - Fine: 0.6, Medium: 0.3, Coarse: 0.15
  # divide by 100 to get proportion
  soilClayCEC <- fertiliserFactors %>%
    mutate(soilClayCEC = if_else(Soil_TopsoilClayContent_1k/100 >= 0.6, 0.6
                                 , if_else(Soil_TopsoilClayContent_1k/100 <= 0.15, 0.15, 
                                           0.3))) %>%
    dplyr::select(rcFid_1km, Soil_TopsoilClayContent_1k, soilClayCEC)
  table(as.vector(soilClayCEC$soilClayCEC))
  
  # calculate CEC based on equation 2.3.5 in CFT
  fertiliserFactors$cec <- (51 * fertiliserFactors$Soil_TopsoilPH_1k - 59) * fertiliserFactors$Soil_TopsoilOrganicCarbonContent_1k/3000000/fertiliserFactors$BD +
    (30 + 4.4 * fertiliserFactors$Soil_TopsoilPH_1k) * soilClayCEC$soilClayCEC
  min(fertiliserFactors$cec, na.rm = T)
  max(fertiliserFactors$cec, na.rm = T) # quite a bit lower than maximum of 250 - 400
  head(fertiliserFactors)
  # place resulting data into 1 of 4 categories
  fertiliserFactors <- fertiliserFactors %>% 
    mutate(cecCat = if_else(cec <= 16, "low CEC"
                            , if_else(cec <= 24, "low-mod CEC"
                                      , if_else(cec <= 32, "high-mod CEC", "high CEC"))))
  table(as.vector(fertiliserFactors$cecCat))
  # get the coefficients from Table 2.4 in CFT
  cecCoefs <- c(low_CEC = c(0,0,0.088)
                , low_modCEC = c(0,0,0.012)
                , high_modCEC = c(0,0,0.163)
                , highCEC = c(0, 0, 0))
  
  ##### 7c5 - SOC #####
  ## ------------ Notes --------------  ##
  ## data are available for this
  ## ------------ ----- --------------  ##
  fertiliserFactors <- fertiliserFactors %>% 
    mutate(socCat = if_else(Soil_TopsoilOrganicCarbonContent_1k < 1, "soc_1"
                            , if_else(Soil_TopsoilOrganicCarbonContent_1k <= 3, "soc_1_3"
                                      , if_else(Soil_TopsoilOrganicCarbonContent_1k <= 6, "soc_3_6", "soc_6"))))
  table(as.vector(fertiliserFactors$socCat))
  # get the coefficients from Table 2.4 in CFT
  socCoefs <- c(soc_1 = c(0,0,0)
                , soc_1_3 = c(0.14, 0, 0)
                , soc_3_6 = c(0.58, 2.571, 0)
                , soc_6 = c(1.045, 2.571, 0))
  
  ##### 7c6 - application method #####
  ## ------------ Notes --------------  ##
  ## only affects NH3 emissions
  ## there are 6 options available, but the method is unknown, so an average
  ## value just for NH3 emissions will be used
  ## ------------ ----- --------------  ##
  applMethodNH3 <- mean(c(-1.292, -1.305, -1.844, -2.465, -1.895))
  
  ##### 7c6 - fertiliser type #####
  ## ------------ Notes --------------  ##
  ## the fertiliser type adds different emissions
  ## The type of fertiliser can change, and will be affected by a land 
  ## manager's preference. Here, only three types are used, as they are recommended
  ## in RB209. These are: Ammonium nitrate, urea, and Urea-ammonium nitrate liquid
  
  ## fertiliser type: based on Table 2.4 in CFT, urea ammonia nitrate solution has the 
  ## same coefficients as ammonium nitrate (granulated); therefore, the two types of fertiliser
  ## used for this will be urea ammonia nitrate solution, and urea - 46% N
  ## ------------ ----- --------------  ##
  
  # get the coefficients from Table 2.4 in CFT
  fertiliserCoefs <- c(nitrateSol = c(0.0053, 0.0004, 0)
                       , urea = c(0.0051, 0.0061, 0.666))
  
  ##### 7c7 - application rate #####
  ## ------------ Notes --------------  ##
  ## the application rate is the amount of a fertiliser applied. This is provided
  ## from the SNS tables produced earlier
  ## ------------ ----- --------------  ##
  
  ##### 7c8 - crop type #####
  ## ------------ Notes --------------  ##
  ## using the crop types in Bouwman et al. (2002), the crops present in this study
  ## can be split into three groups, namely grass (uncropped land), legumes (field beans)
  ## and 'other upland' (barley, maize, oats, potato, wheat, sugarbeet, rapeseed)
  ## ------------ ----- --------------  ##
  
  # get the coefficients from Table 1 in Bouwman et al. (2002)
  cropTypeCoefs <- c(grass = c(-1.268, 0, -0.158)
                     , legumes = c(-0.023, 0, -0.045)
                     , other_upland = c(0, 0, -0.04))
  
  # stop("before lime requirements")
  #### 7d - Calculate lime requirements ####
  for(iscen in 1:(length(scenarios) + 1)){
    
    # if 'iscen' is 1, use 2015 data
    if(iscen == 1){
      cat("Starting 2015 data for lime requirements...\n")
      # set input location - for land cover
      currPath <- file.path("./data_in", "land_cover")
      # for crops
      CropPath <- "./data_in/crops"
      # set output location
      savePath <- file.path("./results", "arable")
      
      # import GE 1 km2 crop area for 2015
      cropsMapIn <-st_read(file.path(currPath, "land_cover_table.gpkg"), quiet = T)
      head(cropsMapIn)
      
    } else { # if > 1, use the scenarios
      scenName <- gsub(".gpkg", "", basename(scenarios[iscen-1]))
      cat("Starting", scenName, "data for lime requirements...\n")
      
      # import GE 1 km2 crop area for 2015
      cropsMapIn <-st_read(file.path(scenPath, paste0(scenName, ".gpkg")), quiet = T)
      head(cropsMapIn)
      
    }
    
    # combine with crop area
    cropSoils <- st_join(cropsMapIn %>%
                           dplyr::select(any_of(c("Name", "lcm", "rcFid_1km"))
                                         , matches("_ha$"))
                         , soils) %>% distinct()
    head(cropSoils)
    
    ##### 7d1 - Determine pH, and thus how much lime may be needed #####
    # get names to analyse
    lcCropNames <- c("rcFid_1km", "winterwheat_ha", "winterbarley_ha", "springwheat_ha", "oats_ha", "maize_ha"
                     , "rapeseed_ha", "springbarley_ha", "potato_ha", "fieldbeans_ha", "sugarbeet_ha"
                     , "improved_grass_ha")
    
    # limeRegs are in t/ha, so determine amount for km2 depending on arable and grass areas
    limeReqsAmounts <- cropSoils %>% dplyr::relocate(all_of(lcCropNames)) %>%
      st_drop_geometry() %>%
      # get total arable area
      mutate(total_arab_ha = rowSums(dplyr::select(., c(winterwheat_ha:sugarbeet_ha)))) %>%
      mutate(total_arab_ha = ifelse(is.na(total_arab_ha), 0, total_arab_ha)) %>%
      rename(total_grass_ha = improved_grass_ha) %>%
      dplyr::select(rcFid_1km, total_arab_ha, total_grass_ha) %>%
      # make spatial again
      merge(cropSoils[, "rcFid_1km"], .) %>%
      # join with advised lime amounts
      st_join(., limeReqs %>% dplyr::select(soilTyp:limeReq_tha_gras)) %>%
      # multiply arable cover (in hectares) by amount of lime per ha
      mutate(limeAmount_tonnes_arab = total_arab_ha * limeReq_tha_arab
             # and do the same for grass
             , limeAmount_tonnes_gras = total_grass_ha * limeReq_tha_gras
             # and add them
             , limeAmount_tonnes = limeAmount_tonnes_arab + limeAmount_tonnes_gras) %>%
      filter(!is.na(soilTyp))
    head(limeReqsAmounts)
    
    # calculate the emissions from lime adage
    # https://naei.beis.gov.uk/data/ef-all-results?q=184177 states that the emission factor 
    # for limestone is 0.12 t CO2e t−1
    limeReqsAmounts$limeTco2e <- limeReqsAmounts$limeAmount_tonnes * 0.12
    head(limeReqsAmounts)
    
    # save
    if(iscen == 1){
      cat("   ...saving", file.path(savePath, "liming_eval.gpkg"), "\n")
      st_write(limeReqsAmounts, file.path(savePath, "liming_eval.gpkg")
               , append = F, quiet = T)
      fwrite(st_drop_geometry(limeReqsAmounts), file.path(savePath, "liming_eval.csv")
             , row.names = F)
    } else {
      cat("   ...saving", file.path(scenarioResultsPath
                                    , paste0("liming_eval", scenName, ".gpkg")), "\n")
      st_write(limeReqsAmounts
               , file.path(scenarioResultsPath
                           , paste0("liming_eval", scenName, ".gpkg"))
               , append = F, quiet = T)
      fwrite(st_drop_geometry(limeReqsAmounts)
             , file.path(scenarioResultsPath
                         , paste0("liming_eval", scenName, ".csv")),
             row.names = F)
    }
  }
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(savePath, "readme_liming_eval.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-06-06
  Last update:",  format(Sys.Date()), "
  
  Description of 'liming_eval.csv':
  This file contains information about about quanitity of lime used in agriculture, and its emissions. 

  Columns of 'liming_eval.csv':
    'rcFid_1km' = grid reference of that cell (x and y coordinates of the centroid point separated by a'-')
    'total_[landtype]_ha' = area, in ha, of either arable land 'arab' or grassland 'grass'
    'soilTyp' = the type of soil each km is considered to be based on which of the substrates was highest in the clay, silt, or sand 
    'limeReq_tha_[landtype]' = tonnes per ha of limestone required for landtype to make the soil more alkaline. This is based on RB209 recommendations
    'limeAmount_tonnes_[landtype]' = tonnes per km2 of limestone required for landtype to make the soil more alkaline 
    'limeAmount_tonnes' = tonnes per km2 of limestone required for both landtypes (arable and grassland) to make the soil more alkaline 
    'limeTco2e' = tonnes of CO2e emissions produced by application of limestone
    
  Data: quantity of limestone emissions and the emissions
  units: see above
  Spatial resolution: 1 km2
  Spatial extent: UK
  Temporal coverage: 2015
  Native projection: 27700")
  
  sink(file = NULL)
  
  # stop("after lime scenarios")
  
  ##### 7d2 - determine fertiliser based on SNS #####
  ## ------------ Notes --------------  ##
  ## The below now does the calculations for the original data and the
  ## scenarios  
  ## ------------ ----- --------------  ##
  
  # concatenate the three categories for fertiliser
  rootDepthUK$snsCats <- paste(rootDepthUK$rainCat
                               , rootDepthUK$soilTyp
                               , rootDepthUK$depthRB209
                               , sep = "_") 
  
  iscen.numb <- 0; iscen.list <- list()
  for(iscen in 1:(length(scenarios) + 1)){
    
    tic("   ...one rotation of converting SNS parameters")
    iscen.numb <- iscen.numb + 1
    
    # if 'iscen' is 1, use 2015 data
    if(iscen == 1){
      
      cat("Starting 2015 data for fertiliser parameters...\n")
      # set input location - for land cover
      currPath <- file.path("./data_in", "land_cover")
      # for crops
      CropPath <- "./data_in/crops"
      # set output location
      savePath <- file.path("./results", "arable")
      
      # import GE 1 km2 crop area for 2015
      cropsMapIn <-st_read(file.path(currPath, "land_cover_table.gpkg"), quiet = T) %>%
        st_centroid() %>%
        mutate(X = as.integer(unlist(map(.$geom, 1))),
               Y = as.integer(unlist(map(.$geom, 2)))) %>%
        # add in current scenario
        mutate(scen = "original_data")
      head(cropsMapIn)
      
      # manure for original land cover
      kgNmanure <- st_read(file.path("results", "animals", "spatialN1km_kgN.gpkg")
                           , quiet = T)
      
    } else { # if > 1, use the scenarios
      
      scenName <- gsub(".gpkg", "", basename(scenarios[iscen-1]))
      cat("\nStarting", scenName, "data for fertiliser parameters...\n")
      
      # import GE 1 km2 crop area for 2015
      cropsMapIn <-st_read(file.path(scenPath, paste0(scenName, ".gpkg")), quiet = T) %>%
        st_centroid() %>%
        mutate(X = as.integer(unlist(map(.$geom, 1))),
               Y = as.integer(unlist(map(.$geom, 2)))) %>%
        # add in current scenario
        mutate(scen = scenName)
      head(cropsMapIn)
      
      # manure for scenarios
      ## use the name to get filter total possible fertiliser, and to read the specific manure output
      nmmanure <- file.path("scenario", "results"
                            , paste0("spatialN1km_kgN", scenName, ".gpkg"))
      kgNmanure <- st_read(file.path("results", "animals", "spatialN1km_kgN.gpkg")
                           , quiet = T)
      
    }
    
    # merge with crop data
    cropSoilRain <- rootDepthUK %>%
      merge(., cropsMapIn %>% st_drop_geometry() %>%
              dplyr::select(rcFid_1km, X, Y, scen
                            , ends_with("_ha"))
            , by = c("X", "Y"), all = T) %>%
      filter(!is.na(rcFid_1km)) %>%
      filter(!is.na(soilTyp))
    head(cropSoilRain)
    
    cat("   ...converting SNS parameters...\n")
    # merge to make spatial
    snsParasWide2 <- snsParasWide %>%
      merge(., cropSoilRain
            , by.x = "threeParas"
            , by.y = "snsCats"
            , all = T) %>%
      dplyr::select(-c(soilTyp, rainCat, rootdepthuk, depthRB209)) %>%
      relocate(threeParas, X, Y, rcFid_1km
               # couple others
               , sns_forage, oats_ha
               , sns_sugarbeet, sugarbeet_ha
               
               # cereals
               , sns_maize, maize_ha
               , sns_winterwheat, winterwheat_ha
               , sns_winterbarley, winterbarley_ha
               , sns_springwheat, springwheat_ha
               , sns_springbarley, springbarley_ha
               
               , sns_potatoes, potato_ha
               , sns_beans, fieldbeans_ha
               , sns_osr, rapeseed_ha
               , sns_uncropped, improved_grass_ha) 
    head(snsParasWide2)
    # make NAs 0 
    snsParasWide2[is.na(snsParasWide2)] <- 0
    
    # use matrix calculations - just areas
    cat("   ...matrix calculations...\n")
    snsMat1 <- snsParasWide2 %>%
      st_drop_geometry() %>%
      dplyr::select(oats_ha
                    , sugarbeet_ha
                    , maize_ha, winterwheat_ha, winterbarley_ha, springwheat_ha, springbarley_ha
                    , potato_ha
                    , fieldbeans_ha
                    , rapeseed_ha
                    , improved_grass_ha) %>%
      as.matrix()
    head(snsMat1)
    
    ## just sns categories
    snsMat2 <- snsParasWide2 %>%
      st_drop_geometry() %>%
      dplyr::select(sns_forage
                    , sns_sugarbeet
                    , sns_maize, sns_winterwheat, sns_winterbarley, sns_springwheat, sns_springbarley
                    , sns_potatoes
                    , sns_beans
                    , sns_osr
                    , sns_uncropped) %>%
      as.matrix()
    head(snsMat2)
    
    ## ------------ Notes --------------  ##
    ## the below code determines how much SNS is necessary for each pixel by 
    ## using the crop area and its SNS class
    ## ------------ ----- --------------  ##
    cat("   ...SNS amount calculations...\n")
    # loop through sns categories
    m3 <- list()
    for(i in 0:3){
      if(i == 0){
        cat("      ...Now:", nm <-paste0("sns", i), "...")
        v = 1
      } else {
        v = v+1
        cat(" |", nm <-paste0("sns", i))
      }
      nm <- paste0(nm, "_")
      
      # new instance
      snsNew <- snsMat2
      # print(head(snsNew))
      # make 1 or 0
      M <- ifelse(snsNew==i,1,0)
      # print(head(M))
      
      # multiply by the area
      # print(head(snsMat1))
      m2 <- M * snsMat1
      # print(head(m2))
      
      # sum the rows
      m3[[v]] <- m2 %>%
        as.data.frame() 
      colnames(m3[[v]]) <- paste0(nm, gsub("sns_", "", colnames(m3[[v]])))
      # print(head(m3[[v]]))
      
    }
    
    # bind back together
    snsParasWide3 <- do.call(cbind, list(snsParasWide2, m3)) %>% as.data.frame() %>%
      # merge back to spatial
      merge(cropSoilRain %>% dplyr::select(rcFid_1km, scen), .)
    head(snsParasWide3)
    cat("\n")
    
    # stop("before writing sns_and_crops_scens.gpkg")
    
    # save
    if(iscen == 1){ # (for original data)
      cat("   ...saving", file.path(CropPath, "sns_and_crops.gpkg"), "\n")
      st_write(snsParasWide3 %>% 
                 filter(scen == "original_data") %>%
                 dplyr::select(-scen)
               , file.path(CropPath, "sns_and_crops.gpkg")
               , append = F, row.names = F, quiet = T)
    } else { # save (for scenario data)
      cat("   ...saving", file.path(scenarioResultsPath
                                    , paste0("sns_and_crops_", scenName, ".gpkg"))
          , "\n")
      st_write(snsParasWide3 %>% 
                 filter(scen != "original_data")
               , file.path(scenarioResultsPath
                           , paste0("sns_and_crops_", scenName, ".gpkg"))
               , append = F, row.names = F, quiet = T)
    }
    
    cat("   ...SNS wide...\n")
    # and with the derived SNS values
    snsParas <- merge(snsParasWide3, snsCropParametersWide 
                      , by = "threeParas"
                      , all = T) %>%
      dplyr::select(-c(contains(c("sns_"))))
    head(snsParas)
    unique(snsParas$scen)
    
    # get all unique names
    startCol <- which(colnames(snsParas) == "sns0_forage")
    endCol <- which(colnames(snsParas) == "sns3_uncropped")
    xCurr <- list()
    toc()
    
    tic("   ...xCurr")
    cat("   ...starting xCurr...\n")
    for(i in startCol:endCol){
      
      # cat(i, nm <- colnames(snsParas)[i], "\n")
      nmKg <- paste0(nm <- colnames(snsParas)[i], "_kgN")
      
      # select only columns that match
      # multiple the amount of N suggested to be applied with the amount of hectares of
      # that crop type
      xCurr[[i]] <- snsParas %>%
        st_drop_geometry() %>%
        dplyr::select(rcFid_1km, contains(nm), scen) %>%
        mutate(total = .[, 2] * .[, 3]) %>%
        dplyr::select(rcFid_1km, scen, total) %>%
        rename(!!nmKg := total)
    }
    toc()
    
    cat("   ...starting xCurr2...\n")
    tic("   ...xCurr2")
    xCurr <- Filter(length, xCurr) # only get filled objects
    toc()
    
    # merge all together
    tic("   ...snsParasKgN")
    cat("   ...starting snsParasKgN...\n")
    # Convert data frames in the list to data.table
    xCurr <- lapply(xCurr, as.data.table)
    
    # Use Reduce with data.table::merge
    ## ------------ Notes --------------  ##
    ## The below 
    ## ??
    ## ------------ ----- --------------  ##
    snsParasKgN <- Reduce(function(x, y) merge(x %>% distinct()
                                               , y %>% distinct()
                                               , by = c("rcFid_1km", "scen")
                                               , all = TRUE, allow.cartesian = F), 
                          xCurr)
    head(snsParasKgN)
    toc(log = T)
    
    # ## ------------ Notes --------------  ##
    # ## The next bit is a back-up method in case the 'Reduce' function above does not work efficiently
    # 
    # # Set a chunk size (e.g., 100,000 rows per chunk)
    # chunk_size <- 100000
    # # get the first x amount of rcFids
    # ## Get the length of unique 'rcFid_1km' values
    # uniqueRcFid <- na.omit(unique(snsParas$rcFid_1km))
    # totalLength <- length(uniqueRcFid)
    # ### Create start sequence
    # rcFidSeq <- seq(1, totalLength, by = chunk_size)
    # ### Create end sequence (ensure the last chunk ends at total_length)
    # rcFidSeq2 <- pmin(rcFidSeq + chunk_size - 1, totalLength)
    # 
    # # Loop over rcFidSeq
    # for (i in seq_along(rcFidSeq)) {
    #   
    #   # Precompute the subset of uniqueRcFid for the current chunk
    #   rcFidSubset <- uniqueRcFid[rcFidSeq[i]:rcFidSeq2[i]]
    #   
    #   # Loop over scenarios
    #   for (sc in scensToRun) {
    #     cat("scen:", sc, "|", "i:", i, "\n")
    #     
    #     # Pre-filter xCurr once for the scenario and rcFidSubset
    #     xCurr_filtered <- lapply(xCurr, function(xx) {
    #       xx[scen == sc & rcFid_1km %in% rcFidSubset]
    #     })
    #     head(xCurr_filtered)
    #     
    #     # Merge the filtered results using data.table's merge function
    #     poutyMerged <- Reduce(function(x, y) merge(x, y, by = c("rcFid_1km", "scen")
    #                                                , all = TRUE
    #                                                , allow.cartesian = TRUE)
    #                           , xCurr_filtered)
    #     
    #     # Start with the first element of xCurr_filtered as the initial merged result
    #     poutyMerged <- xCurr_filtered[[1]]
    #     # Loop through the remaining elements in xCurr_filtered (starting from the second element)
    #     for (j in 2:length(xCurr_filtered)) {
    #       cat(j, "|", dim(poutyMerged), "\n")
    #       poutyMerged <- merge(poutyMerged, xCurr_filtered[[j]] %>% distinct(), 
    #                            by = c("rcFid_1km", "scen"), 
    #                            all = T, 
    #                            allow.cartesian = TRUE)
    #     }
    #     
    #     # Save the result as CSV
    #     fwrite(poutyMerged, paste0("scenario/pm", sc, i, ".csv"), row.names = FALSE)
    #   }
    # }
    # ## ------------ end back-up method --------------  ##
    
    tic("   ...snsParasKgN2")
    # merge with crop data
    snsParasKgN <- snsParasKgN %>%
      # get total
      mutate(totalNperKm = rowSums(select(., 3:last_col()), na.rm = T)) %>%
      # merge back to spatial
      merge(ghgGridGB %>% dplyr::select(rcFid_1km), .)
    head(snsParasKgN)
    toc(log = T)
    
    tic("   ...amount N per hectare")
    # determine the amount per ha
    snsParasKgN.aveHa <- snsParasKgN %>%
      mutate(aveNperHa = totalNperKm / 100) %>%
      st_drop_geometry() %>%
      dplyr::select(aveNperHa) %>%
      {quantile(.$aveNperHa, seq(0, 1, 0.1))}
    snsParasKgN.aveHa
    toc(log = T)
    
    ## ------------ Notes --------------  ##
    ## The next bit uses the manure volume calculated from earlier
    ## ------------ ----- --------------  ##
    # stop("before manure into fertiliser")
    tic("   ...manure and N")
    ##### 7d3 - load in manure amounts, and the N you get from it, and reduce from fertiliser cost #####
    snsParasKgNman <- snsParasKgN %>%
      st_drop_geometry() %>%
      merge(., kgNmanure %>% st_drop_geometry()) %>%
      rename(kgN_fromManure = totalN
             , kgNreq_fromFert = totalNperKm) %>%
      distinct()
    head(snsParasKgNman)
    
    ## get total N added by fertiliser (after removing N obtained from manure)
    snsParasKgNman <- snsParasKgNman %>%
      mutate(finalkgNfert = kgNreq_fromFert - kgN_fromManure) %>%
      mutate(finalkgNfert = ifelse(finalkgNfert < 0, 0, finalkgNfert))
    head(snsParasKgNman)
    
    # save
    if(iscen == 1){ # (for original data)
      cat("   ...saving", file.path(CropPath, "KgNapplied_fert.gpkg"), "\n")
      st_write(snsParasKgNman, file.path(CropPath, "KgNapplied_fert.gpkg"), append = F, quiet = T)
      fwrite(snsParasKgNman %>% st_drop_geometry(), file.path(CropPath, "KgNapplied_fert.csv"), row.names = F)
      
      # check how the data fit with national averages
      # Fertiliser usage on farms: Results from the Farm Business Survey,
      # England 2019/20 state that 109 is the average kg N ha-1
      ## remove really low ones first
      length(which(snsParasKgNman$finalkgNfert / 100 < 10))
      finalkgNfertCheck <- snsParasKgNman$finalkgNfert[snsParasKgNman$finalkgNfert > 10]
      hist(finalkgNfertCheck / 100)
      
      range(finalkgNfertCheck / 100, na.rm = T)
      mean(finalkgNfertCheck / 100, na.rm = T)
      median(finalkgNfertCheck / 100, na.rm = T)
      
    } else { # save (for scenario data)
      cat("   ...saving", file.path(scenarioResultsPath
                                    , paste0("KgNapplied_fert_", scenName, ".gpkg"))
          , "\n")
      st_write(snsParasKgNman, file.path(scenarioResultsPath
                                         , paste0("KgNapplied_fert_", scenName, ".gpkg"))
               , append = F, quiet = T)
      fwrite(snsParasKgNman %>% st_drop_geometry(), file.path(scenarioResultsPath
                                                              , paste0("KgNapplied_fert_", scenName, ".csv"))
             , row.names = F)
    }
    toc(log = T)
  }
  
  # tidy
  rm(rootDepthUK, separatedSoils, snsParasWide3.allScens
     , animalGrid, fulldfLong, animalGridLong)
  gc()
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(CropPath, "readme_sns_and_crops.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
   Date created: 2023-05-05
   Last update:",  format(Sys.Date()), "
   
   Description of 'sns_and_crops':
   These files contain the amount of soil nitrogen supply (SNS) required for a pixel, based on RB209 recommendations and different crop and soil types
   
   Columns of 'sns_and_crops':
    'rcFid_1km' = grid reference of that cell (x and y coordinates of the centroid point separated by a'-')
    'threeParas' = type of km relevant for crop growing, separated by '_', in the order of '[rainfall]_[soil type]_[rootdepth]'
        [rainfall] = 'low', 'medium', or 'high', coresponding to annual rainfall of <600 mm, 600 mm - 700 mm, > 700 mm, respectively
        [soil type] = 'clay', 'silt', or 'sand', based on which was most numerous within a pixel
        [rootdepth] = 'deep', 'medium', or 'shallow' based on RB209 definitions, which took soil type into account
    'X' and  'Y' = x and y coordinates for the centroid point, with an EPSG of 27700
    '[crop_type]_ha' = area, in ha, of that crop type within a 1 km2 cell
    'sns[x]_[crop_type]' = the area of each SNS category 'x' required, in ha, for a certain crop type
    'sns[x]_[crop_type]_kgHa' = the quantity, in kg per ha, of N required to meet each  SNS category 'x' requirement, for a certain crop type
    
   Data: soil nitrogen supply (SNS)
   units: hectare, or kg N / ha
   Spatial resolution: 1 km2
   Spatial extent: UK
   Temporal coverage: 2015
   Native projection: 27700")
  sink(file = NULL)
  
  rm(snsParasWide3)
  gc()
  
  ## ------------ Notes --------------  ##
  ## at this stage, crops and their SNSs have been produced for the original
  ## data and the scenario data
  ## ------------ ----- --------------  ##
  # stop("SNS and crops produced")
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(CropPath, "readme_KgNapplied_fert.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
   Date created: 2023-05-05
   Last update:",  format(Sys.Date()), "
   
   Description of 'KgNapplied_fert':
   These files contain the amount of soil nitrogen supply (SNS) required for a pixel for each agricultural land type, based on RB209 recommendations and different crop and soil types. It also includes N produced in that km by manure management. 
   
   Columns of 'KgNapplied_fert':
    'rcFid_1km' = grid reference of that cell (x and y coordinates of the centroid point separated by a'-')
    'X' and  'Y' = x and y coordinates for the centroid point, with an EPSG of 27700
    'sns[x]_[crop_type]_kgN' = the amount of N, in kg, for each SNS category 'x' required for a certain crop type
    'kgNreq_fromFert' = the quantity, in kg, of N required to meet each  SNS category 'x' requirement, for all crop types
    'kgN_fromManure' = the amount of N, in kg, produced from manure management in that km cell
    'finalkgNfert' = the final amount of N, in kg, required from fertilisers, after removing the N received from manure
    
   Data: soil nitrogen supply (SNS)
   units: kg N per km2
   Spatial resolution: 1 km2
   Spatial extent: UK
   Temporal coverage: 2015
   Native projection: 27700")
  sink(file = NULL)
  
  # stop("before 7e")
  # filter, to keep coefs
  coefKeep <- ls()[grep("coef", ls(), ignore.case = T)]
  # rm(list=setdiff(ls(), c("CropPath", "savePath", "snsParasKgNman"
  #                         , "scenarios", "scenarioResultsPath"
  #                         , "fertiliserFactors"
  #                         , "XYrcFid_1km", "gwpN2O","gwpCH4"
  #                         , "cropSoilRain", "soils"
  #                         , objKeep, coefKeep)))
  rm(list = ls(pattern = "^sns|^lime|^XYr")) # Removes all objects whose name starts with strings
  rm(M, m2, m3) # Removes matrices
  rm(cropsMapIn, cropSoils, xCurr) # Removes crops
  gc()
  
  
  #### 7e - calculate emissions from fertiliser use ####
  ## ------------ Notes --------------  ##
  ## All the emissions factors noted in section 7 require to be summed together, following 
  ## Bouwman et al. (2002). The way described above was the original way to
  ## calculate it (i.e. found in the CFT), and is the focus of section '7e'
  ## ------------ ----- --------------  ##
  
  ###### 7e1 - N2O background ######
  # determine L_back_direct from eq. 2.3.3b in CFT
  # these are termed as 'background N2O emissions' in CFT and
  # are based on the emission factors from classes of:
  # crop type, soil texture, soc, cec, pH, drainage, and application method
  # note application method for N2O is 0, so can be ignored
  
  for(iscen in 1:(length(scenarios) + 1)){
    
    # if 'iscen' is 1, use 2015 data
    if(iscen == 1){
      cat("Starting 2015 data for SNS emissions...\n")
      # set input location - for land cover
      currPath <- file.path("./data_in", "land_cover")
      # for crops
      CropPath <- "./data_in/crops"
      # set output location
      savePath <- file.path("./results", "arable")
      
      # import 1 km2 crop area for 2015
      SNSMapIn <-st_read(file.path(CropPath, "sns_and_crops.gpkg"), quiet = T) %>%
        dplyr::select(-c(rcFid_1km, X, Y)) %>%
        st_join(fertiliserFactors)
      head(SNSMapIn)
      
      # load manure
      snsParasKgNman <- fread(file.path(CropPath, "KgNapplied_fert.csv"))
      
    } else { # if > 1, use the scenarios
      scenName <- gsub(".gpkg", "", basename(scenarios[iscen-1]))
      cat("Starting", scenName, "data for SNS emissions...\n")
      
      # import 1 km2 crop area for scens
      SNSMapIn <-st_read(file.path(scenarioResultsPath, paste0("sns_and_crops_", scenName, ".gpkg"))
                         , quiet = T) %>%
        dplyr::select(-c(rcFid_1km, X, Y)) %>%
        st_join(fertiliserFactors)
      head(SNSMapIn)
      
      # load manure
      snsParasKgNman <- fread(file.path(scenarioResultsPath
                                        , paste0("KgNapplied_fert_", scenName, ".csv"))
      )
    }
    
    # group the things that do not change in a 1 km2 cell (i.e. all but crop type)
    fertClasses <- SNSMapIn %>%
      relocate(rcFid_1km) %>%
      st_drop_geometry() %>%
      # soil texture
      mutate(textureValue = if_else(textureCat == "fine", textureCoefs[["fine1"]]
                                    , if_else(textureCat == "medium", textureCoefs[["medium1"]]
                                              , textureCoefs[["course1"]]))) %>%
      # soc
      mutate(socValue = if_else(socCat == "soc_1", socCoefs[["soc_11"]]
                                , if_else(socCat == "soc_1_3", socCoefs[["soc_1_31"]]
                                          , if_else(socCat == "soc_3_6", socCoefs[["soc_3_61"]]
                                                    , socCoefs[["soc_61"]])))) %>%
      # cec
      mutate(cecValue = if_else(cecCat == "low CEC", cecCoefs[["low_CEC1"]]
                                , if_else(cecCat == "low-mod CEC", cecCoefs[["low_modCEC1"]]
                                          , if_else(cecCat == "high-mod CEC", cecCoefs[["high_modCEC1"]]
                                                    , cecCoefs[["highCEC1"]])))) %>%
      # ph
      mutate(phValue = if_else(phCat == "low pH", phCoefs[["low_ph1"]]
                               , if_else(phCat == "low-mod pH", phCoefs[["low_modph1"]]
                                         , if_else(phCat == "high-mod pH", phCoefs[["high_modph1"]]
                                                   , phCoefs[["highph1"]])))) %>%
      # drainage
      mutate(drainValue = drainageCoefs[[1]]) %>%
      # sum all the outcomes
      mutate(sumN2Ofactors = rowSums(dplyr::select(.[,,drop = T], c(textureValue:drainValue)))) %>%
      # add crop class factors
      ## get total crop area of either grass (uncropped), legumes (beans), or other
      mutate(otherUplandHa = rowSums(dplyr::select(.[,,drop = T], c(oats_ha, sugarbeet_ha
                                                                    , maize_ha, winterwheat_ha, winterbarley_ha
                                                                    , springwheat_ha, springbarley_ha
                                                                    , potato_ha, rapeseed_ha))))
    # Replace NA values with 0 only in numeric columns
    fertClasses <- fertClasses %>%
      mutate(across(where(is.numeric), ~ replace_na(., 0)))
    head(fertClasses)
    
    cat("   ...starting fertClassesSplit...\n")
    # deal with the different crops in the same cell separately
    fertClassesSplit <- fertClasses %>%
      # where otherUplandHa is present, calculate summed factors
      mutate(summedOU = if_else(otherUplandHa > 0, cropTypeCoefs[["other_upland1"]] + sumN2Ofactors, 0)
             # where beans are present, calculate summed factors
             , summedLegumes = if_else(fieldbeans_ha > 0, cropTypeCoefs[["legumes1"]] + sumN2Ofactors, 0)
             # where grass is present, calculate summed factors
             , summedGrass = if_else(improved_grass_ha > 0, cropTypeCoefs[["grass1"]] + sumN2Ofactors, 0)
      ) %>%
      # extract just important columns (i.e. ha of crops, and summed factors)
      dplyr::select(rcFid_1km, otherUplandHa, summedOU, fieldbeans_ha, summedLegumes
                    , improved_grass_ha, summedGrass) %>%
      # add constant, k, and then exponential
      mutate(summedOUtotal = exp(summedOU + -0.41)
             , summedLegumestotal = exp(summedLegumes + -0.41)
             , summedGrasstotal = exp(summedGrass + -0.41)) %>%
      ## ------------ Notes --------------  ##
      ## The above gives background emission values in kg N2O-N, per ha
      ## Each of the crop areas are in hectares, so multiply to get value for 
      ## each km cell
      ## ------------ ----- --------------  ##
      mutate(backgN20_kgkm2 = ((summedOUtotal *  otherUplandHa) +
                                 (summedLegumestotal * fieldbeans_ha) +
                                 (summedGrasstotal * improved_grass_ha))) %>%
      dplyr::select(rcFid_1km, otherUplandHa, fieldbeans_ha
                    , improved_grass_ha, backgN20_kgkm2) 
    head(fertClassesSplit)
    
    ###### 7e2 - N2O direct ######
    cat("   ...starting N2O direct [2.3.3a]...\n")
    # determine L_direct from eq. 2.3.3a in CFT
    # these are termed as 'direct N2O emissions' in CFT and
    # are based on the emission factors from classes of (in addition to the ones for
    # background emissions): application rate, and fertiliser type
    # Bouwman et al. (2002) states that the fertiliser type and the application rate must be used.
    # for simplicity, the ammonium nitrate urea solution will be the only one used in this study
    # due to generally lower emissions
    
    # deal with the different crops and amount of N in the same cell separately
    # this is an interaction term
    fertClassesDirect <- fertClasses %>%
      # fertiliser type - nitrate solution
      mutate(fertValueSolution = fertiliserCoefs[["nitrateSol1"]]) %>%
      # application rate
      merge(., snsParasKgNman %>% st_drop_geometry()
            , by = c("rcFid_1km")) %>%
      # get combined N for 'otherUplandHa' crops - for a Ha
      mutate(NforOtherUp = rowSums(dplyr::select(.[,,drop = T]
                                                 , c(contains(c("maize_kgN", "winterwheat_kgN", "springwheat_kgN"
                                                                , "winterbarley_kgN", "springbarley_kgN"
                                                                , "potatoes_kgN", "sugarbeet_kgN", "osr_kgN"
                                                                , "forage_kgN" # oats
                                                 ))))) / 100
             # get combined N for beans - for a Ha
             , NforBeans = rowSums(dplyr::select(.[,,drop = T]
                                                 , c(contains(c("beans_kgN"))))) / 100
             # get combined N for grass - for a Ha
             , NforGrass = rowSums(dplyr::select(.[,,drop = T]
                                                 , c(contains(c("uncropped_kgN"))))) / 100
             , ) %>%
      # calculate the interaction term
      mutate(NrateFertInterOU = fertValueSolution * NforOtherUp
             , NrateFertInterBean = fertValueSolution * NforBeans
             , NrateFertInterGras = fertValueSolution * NforGrass) %>%
      ## where otherUplandHa is present, calculate summed factors
      mutate(summedOUSol = if_else(otherUplandHa > 0, cropTypeCoefs[["other_upland1"]] + 
                                     sumN2Ofactors + NrateFertInterOU, 0)
             ## where beans are present, calculate summed factors
             , summedLegumesSol = if_else(fieldbeans_ha > 0, cropTypeCoefs[["legumes1"]] + 
                                            sumN2Ofactors + NrateFertInterBean, 0)
             ## where grass is present, calculate summed factors
             , summedGrassSol = if_else(improved_grass_ha > 0, cropTypeCoefs[["grass1"]] + 
                                          sumN2Ofactors + NrateFertInterGras, 0)
      ) %>%
      # extract just important columns (i.e. ha of crops, and summed factors)
      dplyr::select(rcFid_1km, otherUplandHa, summedOUSol, NforOtherUp
                    , fieldbeans_ha, summedLegumesSol, NforBeans
                    , improved_grass_ha, summedGrassSol, NforGrass) %>%
      # add constant, k, and then exponential
      mutate(summedOUtotal = exp(summedOUSol + -0.41)
             , summedLegumestotal = exp(summedLegumesSol + -0.41)
             , summedGrasstotal = exp(summedGrassSol + -0.41)) %>%
      
      ## ------------ Notes --------------  ##
      ## The above gives direct emission values in kg N2O-N, per ha
      ## Each of the crop areas are in hectares, so multiply to get value for 
      ## each km cell
      ## ------------ ----- --------------  ##
      
      mutate(directN20_kgkm2 = ((summedOUtotal *  otherUplandHa) +
                                  (summedLegumestotal * fieldbeans_ha) +
                                  (summedGrasstotal *  improved_grass_ha))) %>%
      dplyr::select(rcFid_1km, otherUplandHa, summedOUtotal
                    , fieldbeans_ha, summedLegumestotal
                    , improved_grass_ha, summedGrasstotal
                    , directN20_kgkm2) 
    head(fertClassesDirect)
    
    ###### 7e3 - N2O direct and background ######
    cat("   ...starting N2O direct and background [2.3.3b]...\n")
    fertN2ODirBackg <- fertClassesDirect %>%
      dplyr::select(rcFid_1km, directN20_kgkm2) %>%
      merge(., fertClassesSplit %>% dplyr::select(rcFid_1km, backgN20_kgkm2) %>%
              st_drop_geometry())
    head(fertN2ODirBackg)
    
    ###### 7e4 - NO direct ######
    cat("   ...starting NO direct [2.3.3a]...\n")
    # NO is calculated in the same way N2O direct, but requires fewer factors: only
    # application rate, fertiliser type, soc, and drainage
    # group the things that do not change in a 1 km2 cell
    
    # From CFT: 'Emissions of NO are also modelled using the same functional form as in
    # equation 2.3.3b and a conversion factor of 0.01 is used to calculate the
    # resulting N2O.'
    
    # separate SNSMapIn into SNs and other factors required for NO calculation
    ## ensure correct order
    SNSMapIn <- SNSMapIn %>%
      relocate(rcFid_1km, x, y) %>%
      # remove sns
      dplyr::select(-c(matches("^sns")))
    
    fertClasses <- SNSMapIn %>% distinct() %>%
      # Replace NA values with 0 only in numeric columns
      mutate(across(where(is.numeric), ~ replace_na(., 0))) %>%
      st_drop_geometry() %>%
      
      # Calculate socValue with case_when for readability and efficiency
      mutate(
        socValue = case_when(
          socCat == "soc_1" ~ socCoefs[["soc_12"]],
          socCat == "soc_1_3" ~ socCoefs[["soc_1_32"]],
          socCat == "soc_3_6" ~ socCoefs[["soc_3_62"]],
          TRUE ~ socCoefs[["soc_62"]]
        )
        , drainValue = drainageCoefs[[2]]
        , fertValueSolution = fertiliserCoefs[["nitrateSol2"]]
        , .keep = "unused") %>%
      
      # Calculate sumN2Ofactors using a consolidated select statement
      mutate(sumN2Ofactors = rowSums(dplyr::select(., socValue, drainValue)))
    head(fertClasses)
    
    fc2 <- fertClasses %>%
      # Join snsParasKgNman and calculate combined N values
      left_join(snsParasKgNman %>% st_drop_geometry() %>%
                  # remove sns
                  dplyr::select(-c(matches("^sns")))
                , by = c("rcFid_1km")) %>%
      # sum all N required for certain categories of plant (crop)
      ## and then remove the columns that were included in the calculations
      mutate(
        NforOtherUp = rowSums(dplyr::select(., contains(c("maize_kgN", "winterwheat_kgN"
                                                          , "springwheat_kgN", "winterbarley_kgN"
                                                          , "springbarley_kgN"
                                                          , "potatoes_kgN", "sugarbeet_kgN"
                                                          , "osr_kgN", "forage_kgN")))) / 100,
        NforBeans = rowSums(dplyr::select(., contains("beans_kgN"))) / 100,
        NforGrass = rowSums(dplyr::select(., contains("uncropped_kgN"))) / 100
        , .keep = "unused") %>%
      # sum total area (ha) for 'other upland' crops
      ## and then remove the columns that were included in the calculations
      mutate(
        otherUplandHa = rowSums(dplyr::select(., oats_ha, sugarbeet_ha, maize_ha
                                              , winterwheat_ha, winterbarley_ha, springwheat_ha
                                              , springbarley_ha, potato_ha, rapeseed_ha))
        , .keep = "unused") %>%
      # Calculate interaction terms
      mutate(
        NrateFertInterOU = fertValueSolution * NforOtherUp,
        NrateFertInterBean = fertValueSolution * NforBeans,
        NrateFertInterGras = fertValueSolution * NforGrass
      ) %>%
      
      # Calculate summed factors based on conditions
      mutate(
        summedOUSol = if_else(otherUplandHa > 0, cropTypeCoefs[["other_upland2"]] + 
                                sumN2Ofactors + NrateFertInterOU, 0),
        summedLegumesSol = if_else(fieldbeans_ha > 0, cropTypeCoefs[["legumes2"]] + 
                                     sumN2Ofactors + NrateFertInterBean, 0),
        summedGrassSol = if_else(improved_grass_ha > 0, cropTypeCoefs[["grass2"]] + 
                                   sumN2Ofactors + NrateFertInterGras, 0)
      ) %>%
      
      # Final selection and exponential transformations
      dplyr::select(rcFid_1km, otherUplandHa, summedOUSol
                    , NforOtherUp, fieldbeans_ha, summedLegumesSol
                    , NforBeans, improved_grass_ha, summedGrassSol, NforGrass) %>%
      
      mutate(
        summedOUtotal = exp(summedOUSol - 1.527),
        summedLegumestotal = exp(summedLegumesSol - 1.527),
        summedGrasstotal = exp(summedGrassSol - 1.527),
        directN0_kgkm2 = (summedOUtotal * otherUplandHa) + 
          (summedLegumestotal * fieldbeans_ha) + 
          (summedGrasstotal * improved_grass_ha),
        N2O_kgkm_fromNO = directN0_kgkm2 * 0.01
      ) %>% 
      dplyr::select(rcFid_1km, otherUplandHa, summedOUtotal
                    , fieldbeans_ha, summedLegumestotal
                    , improved_grass_ha, summedGrasstotal
                    , directN0_kgkm2, N2O_kgkm_fromNO)
    
    # View the result
    head(fc2)
    
    ###### 7e5 - N2O direct and background and from NO ######
    cat("   ...starting NO direct and background [2.3.3b]...\n")
    fertN2Os <- fertN2ODirBackg %>%
      merge(., fc2 %>% dplyr::select(rcFid_1km, N2O_kgkm_fromNO) %>%
              st_drop_geometry())
    head(fertN2Os)
    
    ###### 7e6 - NH3 from volatisation ######
    cat("   ...starting NH3 from volatisation [2.3.4]...\n")
    # NH3 is calculated in a similar way to N2O direct, but requires fewer factors: 
    # crop type, soil texture, cec, pH, drainage, fertiliser type, and application method
    
    # group the things that do not change in a 1 km2 cell (i.e. all but crop type)
    fertClassesNH3 <- SNSMapIn %>% st_drop_geometry() %>% relocate(rcFid_1km) %>%
      # Replace NA values with 0 only in numeric columns
      mutate(across(where(is.numeric), ~ replace_na(., 0))) %>%
      # soil texture
      mutate(textureValue = if_else(textureCat == "fine", textureCoefs[["fine3"]]
                                    , if_else(textureCat == "medium", textureCoefs[["medium3"]]
                                              , textureCoefs[["course3"]]))
             , .keep = "unused") %>%
      # cec
      mutate(cecValue = if_else(cecCat == "low CEC", cecCoefs[["low_CEC3"]]
                                , if_else(cecCat == "low-mod CEC", cecCoefs[["low_modCEC3"]]
                                          , if_else(cecCat == "high-mod CEC", cecCoefs[["high_modCEC3"]]
                                                    , cecCoefs[["highCEC3"]])))
             , .keep = "unused") %>%
      # ph
      mutate(phValue = if_else(phCat == "low pH", phCoefs[["low_ph3"]]
                               , if_else(phCat == "low-mod pH", phCoefs[["low_modph3"]]
                                         , if_else(phCat == "high-mod pH", phCoefs[["high_modph3"]]
                                                   , phCoefs[["highph3"]])))
             , .keep = "unused") %>%
      # drainage
      mutate(drainValue = drainageCoefs[[3]]) %>%
      # application method
      mutate(applValue = applMethodNH3) %>%
      # fertiliser type
      mutate(fertValue = fertiliserCoefs[["nitrateSol3"]]) %>%
      # sum all the outcomes
      mutate(sumNH3factors = rowSums(dplyr::select(., c(textureValue:fertValue))))
    head(fertClassesNH3)
    
    # deal with the different crops in the same cell separately
    fertClassesNH3crop <- fertClassesNH3 %>%
      ## get total crop area of either grass (uncropped), legumes (beans), or other
      mutate(otherUplandHa = rowSums(dplyr::select(., c(oats_ha, sugarbeet_ha
                                                        , maize_ha, winterwheat_ha
                                                        , winterbarley_ha, springwheat_ha
                                                        , springbarley_ha
                                                        , potato_ha, rapeseed_ha)))) %>%
      # where otherUplandHa is present, calculate summed factors
      mutate(summedOU = if_else(otherUplandHa > 0, cropTypeCoefs[["other_upland3"]] + sumNH3factors, 0)
             # where beans are present, calculate summed factors
             , summedLegumes = if_else(fieldbeans_ha > 0, cropTypeCoefs[["legumes3"]] + sumNH3factors, 0)
             # where grass is present, calculate summed factors
             , summedGrass = if_else(improved_grass_ha > 0, cropTypeCoefs[["grass3"]] + sumNH3factors, 0)) %>%
      # add constant, k, and then exponential
      mutate(summedOUtotal = exp(summedOU)
             , summedLegumestotal = exp(summedLegumes)
             , summedGrasstotal = exp(summedGrass)) %>%
      ## ------------ Notes --------------  ##
      ## unlike N2O and NO, NH3 needs to be multiplied by application rate *after* the
      ## sum of the factors has taken place
      ## Therefore, the exact amount of N per ha needs to be determined. This can be done
      ## by back-calculating the SNS classes
      ## ------------ ----- --------------  ##
      merge(., snsParasKgNman) %>%
      dplyr::select(c(rcFid_1km, threeParas
                      , oats_ha, maize_ha, winterwheat_ha, winterbarley_ha
                      , springwheat_ha, springbarley_ha, otherUplandHa
                      , potato_ha, rapeseed_ha, fieldbeans_ha, improved_grass_ha, sugarbeet_ha
                      , contains(c("sns0_", "sns1_", "sns2_", "sns4_")))
                    , summedOUtotal, summedLegumestotal, summedGrasstotal) %>%
      # get rowsums of different crops in terms of N application
      mutate(sugarbeet_kgN = rowSums(dplyr::select(.[,,drop = T], c(contains("sugarbeet_kgN"))))
             , winterbarley_kgN = rowSums(dplyr::select(.[,,drop = T], c(contains("winterbarley_kgN"))))
             , springbarley_kgN = rowSums(dplyr::select(.[,,drop = T], c(contains("springbarley_kgN"))))
             , winterwheat_kgN = rowSums(dplyr::select(.[,,drop = T], c(contains("winterwheat_kgN"))))
             , springwheat_kgN = rowSums(dplyr::select(.[,,drop = T], c(contains("springwheat_kgN"))))
             , maize_kgN = rowSums(dplyr::select(.[,,drop = T], c(contains("maize_kgN"))))
             , potatoes_kgN = rowSums(dplyr::select(.[,,drop = T], c(contains("potatoes_kgN"))))
             , beans_kgN = rowSums(dplyr::select(.[,,drop = T], c(contains("beans_kgN"))))
             , osr_kgN = rowSums(dplyr::select(.[,,drop = T], c(contains("osr_kgN"))))
             , oats_kgN = rowSums(dplyr::select(.[,,drop = T], c(contains("forage_kgN"))))
             , uncropped_kgN = rowSums(dplyr::select(.[,,drop = T], c(contains("uncropped_kgN"))))
      ) %>%
      # divide amount of N by hectare of each crop
      mutate(sugarbeet_kgN_ha = sugarbeet_kgN / sugarbeet_ha
             , winterbarley_kgN_ha = winterbarley_kgN / winterbarley_ha
             , springbarley_kgN_ha = springbarley_kgN / springbarley_ha
             , winterwheat_kgN_ha = winterwheat_kgN / winterwheat_ha
             , springwheat_kgN_ha = springwheat_kgN / springwheat_ha
             , maize_kgN_ha = maize_kgN / maize_ha
             , potatoes_kgN_ha = potatoes_kgN / potato_ha
             , beans_kgN_ha = beans_kgN / fieldbeans_ha
             , osr_kgN_ha = osr_kgN / rapeseed_ha
             , oats_kgN_ha = oats_kgN / oats_ha
             , uncropped_kgN_ha = uncropped_kgN / improved_grass_ha
      ) %>% 
      # Replace NA values with 0 only in numeric columns
      mutate(across(where(is.numeric), ~ replace_na(., 0))) %>%
      dplyr::select(c(rcFid_1km, threeParas
                      , oats_ha, maize_ha, winterwheat_ha, winterbarley_ha, springwheat_ha, springbarley_ha
                      , potato_ha, rapeseed_ha, fieldbeans_ha, improved_grass_ha, sugarbeet_ha
                      , contains(c("kgN_ha")))
                    , summedOUtotal, summedLegumestotal, summedGrasstotal) %>%
      # calculate NH3 emissions for all different crops based on the N input
      # area (in ha) * emissions from factors * N amount per ha
      mutate(sugarbeet_NH3_km2 = sugarbeet_ha * summedOUtotal * sugarbeet_kgN_ha
             , winterbarley_NH3_km2 = winterbarley_ha * summedOUtotal * winterbarley_kgN_ha
             , springbarley_NH3_km2 = springbarley_ha * summedOUtotal * springbarley_kgN_ha
             , winterwheat_NH3_km2 = winterwheat_ha * summedOUtotal * winterwheat_kgN_ha
             , springwheat_NH3_km2 = springwheat_ha * summedOUtotal * springwheat_kgN_ha
             , maize_NH3_km2 = maize_ha * summedOUtotal * maize_kgN_ha
             , potatoes_NH3_km2 = potato_ha * summedOUtotal * potatoes_kgN_ha
             , beans_NH3_km2 = fieldbeans_ha * summedLegumestotal * beans_kgN_ha
             , osr_NH3_km2 = rapeseed_ha * summedOUtotal * osr_kgN_ha
             , oats_NH3_km2 = oats_ha * summedOUtotal * oats_kgN_ha
             , uncropped_NH3_km2 = improved_grass_ha * summedGrasstotal * uncropped_kgN_ha
      ) %>%
      # sum together, to get overall value for km2
      mutate(NH3_km2 = rowSums(dplyr::select(.[,,drop = T], c(contains("NH3_km2"))))) %>%
      # convert to N2O using a 0.01 conversion
      mutate(N2O_kgkm2_fromNH3 = NH3_km2 * 0.01) %>%
      dplyr::select(c(rcFid_1km
                      , oats_ha, maize_ha, winterwheat_ha, winterbarley_ha, springwheat_ha, springbarley_ha
                      , potato_ha, rapeseed_ha, fieldbeans_ha, improved_grass_ha, sugarbeet_ha
                      , NH3_km2, N2O_kgkm2_fromNH3))
    head(fertClassesNH3crop)
    
    ###### 7d7 - N2O direct and background and from NO and from NH3 (volatisation) ######
    fertN2Os <- fertN2Os %>% distinct() %>%
      merge(., fertClassesNH3crop %>% dplyr::select(rcFid_1km, N2O_kgkm2_fromNH3) %>%
              st_drop_geometry())
    head(fertN2Os)
    
    ###### 7d8 - NH3 from leaching ######
    cat("   ...starting NH3 from leaching [2.3.6]...\n")
    # EQUATION 11.11 in IPCC (2019)
    # sum(A (fertiliser rate [kg N]) * FracGASF (leaching fraction [Volatilisation from synthetic fertiliser]))
    # FracGASF is 0.08, according to Table 11.3 in IPCC (2019)
    # multiplied by emission factor for N2O emissions from atmospheric deposition of N on soils and water, which
    # is 0.014 in IPCC (2019)
    
    leachingNH3 <- snsParasKgNman %>%
      dplyr::select(rcFid_1km, finalkgNfert) %>%
      # calculate leaching
      mutate(kgN2O_N_leach = (finalkgNfert * 0.08) * 0.014)
    head(leachingNH3)
    
    ###### 7d9 - N2O from NO, NH3 (volatisation), leaching ######
    cat("   ...combining fertiliser emissions...\n")
    fertN2Os <- fertN2Os %>%
      merge(., leachingNH3 %>% dplyr::select(rcFid_1km, kgN2O_N_leach) %>%
              st_drop_geometry())
    # all of the above are in the N2O-N form. To convert these, the equation is:
    # N2O-N * 44/28 = N2O
    fertN2Os_N2O <- fertN2Os %>%
      mutate(across(directN20_kgkm2:kgN2O_N_leach, ~ . * 44/28, .names = "finN2O_{.col}")) %>%
      # get final value for combination of background and direct emissions, from volatisation, and from leaching
      mutate(totalN20_kgkm2 = finN2O_directN20_kgkm2 + finN2O_N2O_kgkm_fromNO + 
               finN2O_N2O_kgkm2_fromNH3 + finN2O_kgN2O_N_leach) %>%
      # to get the amount in CO2e, times by 265
      mutate(totalCO2e_kgkm2 = totalN20_kgkm2 * gwpN2O)
    
    # add liming CO2e values in
    # if 'iscen' is 1, use 2015 data
    if(iscen == 1){
      liming <- st_read(file.path(savePath, "liming_eval.gpkg"), quiet = T)
    } else {
      liming <- st_read(file.path(scenarioResultsPath
                                  , paste0("liming_eval", scenName, ".gpkg")), quiet = T)
    }
    fertCO2e <- fertN2Os_N2O %>%
      merge(., liming %>% st_drop_geometry()) %>%
      # add lime to total
      mutate(fertliser_CO2e_kgkm2 = totalCO2e_kgkm2 + (limeTco2e * 1000) # to make kg
      ) %>%
      dplyr::select(c(rcFid_1km: kgN2O_N_leach
                      , limeAmount_tonnes, limeTco2e
                      , fertliser_CO2e_kgkm2))
    head(fertCO2e)
    
    # write
    # if 'iscen' is 1, use 2015 data
    if(iscen == 1){
      st_write(fertCO2e, xSave <- file.path(savePath, "fertiliser_emissions.gpkg"), append = F, quiet = T)
      fwrite(fertCO2e %>% st_drop_geometry(), file.path(savePath, "fertiliser_emissions.csv"), row.names = F)
    } else { # scenario data
      st_write(fertCO2e, xSave <- file.path(scenarioResultsPath
                                            , paste0("fertiliser_emissions", scenName, ".gpkg"))
               , quiet = T
               , append = F)
      fwrite(fertCO2e %>% st_drop_geometry()
             , file.path(scenarioResultsPath
                         , paste0("fertiliser_emissions", scenName, ".csv"))
             , row.names = F)
    }
    cat("      ...saved the final fertiliser CO2e product:"
        , xSave, "\n\n")
  } # end of scenarios
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(savePath, "readme_fertiliser_emissions.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2023-05-03
Last update:",  format(Sys.Date()), "

The files containing 'fertiliser_emissions' contain information the emissions, in the form of CO2 equivalents, produced by the use of potential fertilisers per km2.
These are potential emissions, and were based on the guidance in RB209, which stated what one would need to do to get the highest potential yield. However, it should be noted that this does not take in financial considerations
The final values in the fertiliser_emissions files are the combination of direct (i.e. fertiliser-induced) and background (soil) emissions, emissions from volatisation, leaching, and addition of limestone
The values where derived in the form of N2O equivalents, but have been transformed to CO2e.

The columns:
      'rcFid_1km' = 1 km2 pixel identifier, x and y centroid coordinates of 1 km pixel (EPSG: 27700)
      'directN20_kgkm2' = emissions, in kg N2O-N, from direct (i.e. fertiliser-induced) and soil sources
      'backgN20_kgkm2' = emissions, in kg N2O-N, from soil sources
      'directN0_N2O_kgkm2' = emissions, in kg N2O-N, from direct (i.e. fertiliser-induced) sources. This is the NO portion, which has then been converted to the equivalent N2O emissions (using a 0.01 conversion)
      'N2O_kgkm2_fromNH3' = emissions, in kg N2O-N, from volatilisation of ammonia. These were originally calculated in terms of NH3, which was then converted to the equivalent N2O emissions (using a 0.01 conversion)
      'kgN2O_N_leach' = emissions, in kg N2O-N, from leaching of N
      'limeAmount_tonnes' = amount, in metric tonnes, of lime applied to a km2 to get highest potential yield, as calculated based on RB209 advice
      'limeTco2e' = emissions, in t co2e, from the application of limestone
      'fertliser_CO2e_kgkm2' = total co2e emissions, in kg co2e, from all of the aspects of fertilisers that create emissions
          
Data:
    units: kg CO2e emissions per km2 (for 'fertliser_CO2e_kgkm2') 
    Spatial resolution: 1 km2
    Spatial extent: GB
    Temporal coverage: 2015/2016
    Native projection: 27700")
  
  sink(file = NULL)
  
  # tidy
  rm(part7fertiliserUse)
  
} # end of fertiliser (part 7)

#### 8 - part8landuse ####
if(part8landuse){
  # stop("start of part 8")
  # set input location - for land cover
  currPath <- file.path("./data_in", "land_cover")
  # and for crops
  CropPath <- "./data_in/crops"
  # results
  savePath <- "./results/land_use"
  
  # import 1 km2 crop area for 2015
  landArea.2015 <- st_read(file.path(currPath, "land_cover_table.gpkg")) %>%
    relocate(improved_grass_ha, .before=winterwheat_ha) 
  head(landArea.2015)
  
  # import 1 km2 crop area for scenarios
  landArea.names <- list.files(file.path("scenario", "scen_maps", "finer_detail")
                               , pattern = "_land.gpkg"
                               , full.names = T)
  landArea.scen <- pblapply(landArea.names, function(x) {
    st_read(x) %>%
      relocate(improved_grass_ha, .before=winterwheat_ha) 
  })
  head(landArea.scen[[1]])
  
  ##### 8a - land cover for 2015 #####
  # make all crops into 'arable'
  # Identify the columns ending with '_ha'
  hectareCols1 <- which(colnames(landArea.2015) == "winterwheat_ha")
  hectareCols2 <- which(colnames(landArea.2015) == "sugarbeet_ha")
  
  cropTableHa.2015 <- landArea.2015 %>%
    replace(is.na(.), 0) %>%
    mutate(Arable_ha = rowSums(.[, hectareCols1:hectareCols2, drop=TRUE], na.rm = TRUE)) %>%
    dplyr::select(-c(winterwheat_ha:sugarbeet_ha)) %>%
    # assign correct order
    relocate(broadleaf_ha, conifer_ha, Arable_ha, improved_grass_ha
             , neutral_grass_ha, calc_grass_ha, acid_grass_ha, fen_marsh_ha, heather_ha               
             , heather_grass_ha, bog_ha, inland_rock_ha, saltwater_ha, freshwater_ha, sup_lit_rock_ha
             , sup_lit_sed_ha, lit_rock_ha, lit_sed_ha, saltmarsh_ha, urban_ha, suburban_ha)
  head(cropTableHa.2015)
  
  # save
  st_write(cropTableHa.2015, file.path(currPath, "land_cover_area.gpkg"), append = F)
  fwrite(st_drop_geometry(cropTableHa.2015)
         , file.path(currPath, "land_cover_area.csv")
         , row.names = F)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(currPath, "readme_land_cover_area.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2023-05-03
Last update:",  format(Sys.Date()), "

The files containing 'land_cover_area' contain the area, in hectares, of different land covers, derived from CEH's LCM2015 and Land Cover® plus: Crops 2016 (Rowland et al., 2017)
The columns:
      'rcFid' = 1 km reference id
      '[land cover type]' = area of the land covered, in ha, by the title land cover type in a 1 km2 pixel
          land covers:
           arable = winter wheat + winter barley + spring wheat + oats + maize + rapeseed + spring barley + potato + field beans + sugat beet
           improved_grass = improved grassland
           broadleaf = broadleaved woodland
           conifer = coniferous woodland
           neutral_grass = neutral grassland
           calc_grass = calcareous grassland
           acid_grass = acid grassland
           fen_marsh = fen, marsh, and swamp
           heather = heather
           heather_grass = heather grassland
           bog = bog
           inland_rock = inland rock
           saltwater = saltwater
           freshwater = freshwater
           sup_lit_rock = supra-littoral rock
           sup_lit_sed = supra-littoral sediment
           lit_rock = littoral rock
           lit_sed = littoral sediment
           saltmarsh = saltmarsh
           urban = urban
           suburban = suburban
units: ha/km2 
Spatial resolution: 1 km2
Spatial extent: England
Temporal coverage: 2016
Native projection: 27700
      
citations: Rowland, C.S.; Morton, R.D.; Carrasco, L.; McShane, G.; O'Neil, A.W.; Wood, C.M. (2017). Land Cover Map 2015 (1km dominant aggregate class, GB). NERC Environmental Information Data Centre. https://doi.org/10.5285/711c8dc1-0f4e-42ad-a703-8b5d19c92247 ")
  
  sink(file = NULL)
  
  ##### 8a2 - land cover for scenarios #####
  # create new land cover directory
  dir.create(file.path("scenario", "scen_maps", "land_cover")
             , showWarnings = F)
  
  # make all crops into 'arable'
  # Identify the columns ending with '_ha'
  hectareCols1 <- which(colnames(landArea.scen[[1]]) == "winterwheat_ha")
  hectareCols2 <- which(colnames(landArea.scen[[1]]) == "sugarbeet_ha")
  
  for(x in 1:length(landArea.scen)){
    # name
    nm <- gsub("_land.gpkg", "", basename(landArea.names[[x]]))
    x2 <- landArea.scen[[x]] %>%
      replace(is.na(.), 0) %>%
      mutate(Arable_ha = rowSums(.[, hectareCols1:hectareCols2, drop=TRUE], na.rm = TRUE)) %>%
      dplyr::select(-c(winterwheat_ha:sugarbeet_ha)) %>%
      # assign correct order
      relocate(broadleaf_ha, conifer_ha, Arable_ha, improved_grass_ha
               , neutral_grass_ha, calc_grass_ha, acid_grass_ha, fen_marsh_ha, heather_ha               
               , heather_grass_ha, bog_ha, inland_rock_ha, saltwater_ha, freshwater_ha, sup_lit_rock_ha
               , sup_lit_sed_ha, lit_rock_ha, lit_sed_ha, saltmarsh_ha, urban_ha, suburban_ha)
    # save
    st_write(x2, file.path("scenario", "scen_maps", "land_cover"
                           , paste0("land_cover_area", nm, ".gpkg"))
             , append = F)
  }
  
  # write readme
  sink(file = NULL)
  sink(file = file.path("scenario", "scen_maps", "land_cover", "readme_land_cover_area.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2023-05-03
Last update:",  format(Sys.Date()), "

The files containing 'land_cover_area' contain the area, in hectares, of different land covers, derived from Redhead et al. (2020).
The prefixes of the file names relate to specific scenarios taken from Redhead et al. (2020):
  'Ag_15' and 'Ag_30' = 15 and 30% agricultural expansion from the 2007 baseline
  'Gr_15' and 'Gr_30' = 15 and 30% grassland restoration from the 2007 baseline
  'baseline2007' = baseline land cover from 2007
The maps were specifically used from the baseline, 15% and 30% increase in grassland, and 15% increase in agriculutral expansion.
The columns:
      'rcFid' = 1 km reference id
      '[land cover type]' = area of the land covered, in ha, by the title land cover type in a 1 km2 pixel
          land covers:
           arable = winter wheat + winter barley + spring wheat + oats + maize + rapeseed + spring barley + potato + field beans + sugat beet
           improved_grass = improved grassland
           broadleaf = broadleaved woodland
           conifer = coniferous woodland
           neutral_grass = neutral grassland
           calc_grass = calcareous grassland
           acid_grass = acid grassland
           fen_marsh = fen, marsh, and swamp
           heather = heather
           heather_grass = heather grassland
           bog = bog
           inland_rock = inland rock
           saltwater = saltwater
           freshwater = freshwater
           sup_lit_rock = supra-littoral rock
           sup_lit_sed = supra-littoral sediment
           lit_rock = littoral rock
           lit_sed = littoral sediment
           saltmarsh = saltmarsh
           urban = urban
           suburban = suburban
units: ha/km2 
Spatial resolution: 1 km2
Spatial extent: England
Temporal coverage: 2016
Native projection: 27700
      
citations: Redhead, J. W., Powney, G. D., Woodcock, B. A., & Pywell, R. F. (2020). Effects of future agricultural change scenarios on beneficial insects. Journal of Environmental Management, 265, 110550. https://doi.org/10.1016/j.jenvman.2020.110550")
  
  sink(file = NULL)
  
  #### 8b - add in margins ####
  ## ------------ Notes --------------  ##
  ## each km2 has different amounts of field margins.
  ## for ease of interpretation, these are will be subtracted from 'arable' land
  ## areas, if possible. If not, they will come off improved grassland.
  ## ------------ ----- --------------  ##
  
  # if ccs and ees data has not already been created, import
  if (!file.exists(file.path(currPath, "margin_area.gpkg"))){
    
    # the countryside survey has field margin data, so import it
    ccs_aes <- read.csv("N:\\Data\\UK\\countryside_survey/Countryside_Stewardship_Scheme_2016_Management_Options_(England).csv")
    # and environmental stewardship schemes AES 
    ees_aes <- read.csv("N:\\Data\\UK\\AES/Environmental_Stewardship_Scheme_Options_(England).csv")
    
    # check
    str(ccs_aes)
    # opt_desc column important
    str(ees_aes)
    # 'opttitle' column important
    
    # reduce to just margin data
    ccs_margin <- ccs_aes[grepl("margin|strip", ccs_aes$opt_desc), ]
    ees_margin <- ees_aes[grepl("margin|strip", ees_aes$opttitle), ]
    # tidy
    rm(ccs_aes, ees_aes)
    
    # determine unique margins
    # and remove any that are not relevant
    print(list(unique(ccs_margin$opt_desc)))
    print(list(unique(ees_margin$opttitle)))
    ees_margin <- ees_margin[!ees_margin$opttitle %in% c("Cultivated fallow plots or margins for arable plants"
                                                         , "Uncropped, cultivated margins for rare plants on arable land"
                                                         , "Uncropped, cultivated margins for rare plants"), ]
    
    # make point data
    ees_mar_point <- st_as_sf(ees_margin, coords = c("X", "Y"), crs = 27700)
    ccs_mar_point <- st_as_sf(ccs_margin, coords = c("X", "Y"), crs = 27700)
    # tidy
    # rm(ees_margin, ccs_margin)
    
    # only keep necessary cols
    ees_mar_point <- ees_mar_point[, c("OBJECTID", "optcode", "opttitle", "areaha")]
    ccs_mar_point <- ccs_mar_point[, c("OBJECTID", "opt_code", "opt_desc", "quantity")]
    # make both a common format
    colnames(ees_mar_point) <- colnames(ccs_mar_point)
    # convert
    ccs_mar_point$quantity <- as.numeric(ccs_mar_point$quantity)
    
    # spatial join, removing geometry first
    marginPoints <- rbind(ees_mar_point %>% as.data.frame(), ccs_mar_point %>% as.data.frame())
    # make spatial again
    marginPoints <- st_as_sf(marginPoints) %>%
      rename(area_ha = quantity)
    
    # type of margin
    marginPoints$type <- ifelse(grepl("tree", marginPoints$opt_desc), "tree", "flower")
    head(marginPoints)
    st_crs(marginPoints)
    
    # using cropTable, add the points to a polygon, keeping the information from both
    margPolyInter <- st_intersection(marginPoints, cropTableHa) %>%
      st_as_sf() %>%
      # drop spatial
      st_drop_geometry() %>%
      # keep only necessary
      dplyr::select(rcFid_1km, area_ha, type)
    head(margPolyInter)
    unique(margPolyInter$type)
    
    # sum for each rcFid based on area and type
    marginSum <- margPolyInter %>%
      group_by(rcFid_1km, type) %>%
      # make tonnes from kg
      summarise(total_area = sum(area_ha, na.rm = T)) %>% as.data.frame() %>%
      # pivot wider, to separate flower and tree
      pivot_wider(names_from = type, values_from = total_area) %>%
      replace(is.na(.), 0)
    
    # merge back to spatial
    marginSum <- marginSum %>%
      rename(flowerMargin_ha = flower
             , treeMargin_ha = tree) %>%
      merge(cropTableHa %>% dplyr::select(rcFid_1km), .)
    str(marginSum)
    head(marginSum)
    
    # save
    st_write(marginSum, file.path(currPath, "margin_area.gpkg"), append = F)
    fwrite(st_drop_geometry(marginSum)
           , file.path(currPath, "margin_area.csv")
           , row.names = F)
    
    # write readme
    sink(file = NULL)
    sink(file = file.path(currPath, "readme_margin_area.md"))
    cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2023-05-03
Last update:",  format(Sys.Date()), "

The files containing 'margin_area' contain the area, in hectares, of tree and flower field margins. This information was derived from countryside survey data, and environmental stewardship schemes AES data. Each margin was classified as either 'tree' or 'flower', indicating a field margin composed of trees and flowers, respectively.
The columns:
      'rcFid_1km' = 1 km reference id
      'flowerMargin_ha' = area of the land covered by flower field margins in a 1 km2 pixel (ha)
      'treeMargin_ha' = area of the land covered by tree field margins in a 1 km2 pixel (ha)
units: ha/km2 
Spatial resolution: 1 km2
Spatial extent: England
Temporal coverage: 2016
Native projection: 27700
      
the 'ccs_ees_margins.gpkg' file includes data that were obtained from Countryside Stewardship Scheme 2016 Management Options (England), which were provided by Natural England [url: https://naturalengland-defra.opendata.arcgis.com/datasets/countryside-stewardship-scheme-2016-management-options-england/explore?location=52.688216%2C-2.434032%2C6.68], under an Open Government Licence - https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3.
  It contains information about the position and type of field margins in England.")
    
    sink(file = NULL)
  } else {
    # read in, if it already exists
    marginSum <- st_read(file.path(currPath, "margin_area.gpkg"))
  }
  
  #### 8c - add in hedgerows ####
  ## ------------ Notes --------------  ##
  ## each km2 has different amounts of hedgerows.
  ## for ease of interpretation, these are will be subtracted from 'improved grassland' 
  ## areas, if possible. If not, they will come off 'arable' lands.
  ## ------------ ----- --------------  ##
  
  # first, the hedge area needs to be determined
  hedge1km <- raster("N:/Data/UK/hedge/HedgeLength.grd")
  crs(hedge1km) <- CRS('+init=EPSG:27700')
  # convert to points
  hedgePoints <- rasterToPoints(hedge1km, spatial = T) %>%
    # convert to sf
    st_as_sf(crs = 27700) %>%
    # give column name
    rename(hedgeLengthMperKm2 = 1)
  
  # To get ground area of hedgerow, multiply width by total length in a km
  hedgePoints$hedgeAreaMperKm2 <- hedgePoints$hedgeLengthMperKm2 * 3.4
  
  # to get volume
  hedgePoints$hedgeVolumeMperKm2 <- hedgePoints$hedgeLengthMperKm2 * 3.4 * 3.5
  
  str(hedgePoints)
  head(hedgePoints)
  
  # save as vector
  st_write(hedgePoints, file.path(currPath, "hedge_dims.gpkg"), append = F)
  
  # using cropTable, add the points to a polygon, keeping the information from both
  hedgePolyInter <- st_intersection(hedgePoints, cropTableHa.2015) %>%
    st_as_sf() %>%
    # remove points
    st_drop_geometry() %>%
    # keep only necessary
    dplyr::select(rcFid_1km, hedgeAreaMperKm2, hedgeVolumeMperKm2) %>%
    # merge into polygons
    merge(cropTableHa.2015 %>% dplyr::select(rcFid_1km), .)
  head(hedgePolyInter)
  
  # convert m2 to ha
  hedgePolyInter <- hedgePolyInter %>%
    mutate(hedgeAreaHaperKm2 = hedgeAreaMperKm2 / 10000)
  
  str(hedgePolyInter)
  head(hedgePolyInter)
  
  # save
  st_write(hedgePolyInter, file.path(currPath, "hedgerow_area.gpkg"), append = F)
  fwrite(st_drop_geometry(hedgePolyInter)
         , file.path(currPath, "hedgerow_area.csv")
         , row.names = F)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(currPath, "readme_hedgerow_area.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2023-05-03
Last update:",  format(Sys.Date()), "

The files containing 'hedgerow_area' contain the area and derived volume, in hectares and m3 respectively, of hedgerow area and volume in GB. This information was derived from CEH's Woody Linear Features dataset [url: https://catalogue.ceh.ac.uk/documents/d7da6cb9-104b-4dbc-b709-c1f7ba94fb16] (Scholefield et al. (2016))
The dataset provides the length of woody features. The area and volume were derived from average hedgerow measurements in England.

The columns:
      'rcFid_1km' = 1 km reference id
      'hedgeAreaMperKm2' = area of the land covered by hedgerows in a 1 km2 pixel (m2)
      'hedgeAreaHaperKm2' = area of the land covered by hedgerows in a 1 km2 pixel (ha)
      'hedgeVolumeMperKm2' = volume of hedgerows in a 1 km2 pixel (m3)
units: ha/km2 (for 'hedgeAreaHaperKm2')
Spatial resolution: 1 km2
Spatial extent: GB
Temporal coverage: 2016
Native projection: 27700
      
the 'HedgeLength.grd' file includes data that were obtained from Scholefield et al. (2016)
      
citation: Scholefield, P.A.; Morton, R.D.; Rowland, C.S.; Henrys, P.A.; Howard, D.C.; Norton, L.R. (2016). Woody linear features framework, Great Britain v.1.0. NERC Environmental Information Data Centre. https://doi.org/10.5285/d7da6cb9-104b-4dbc-b709-c1f7ba94fb16")
  sink(file = NULL)
  
  ##### 8d - combine land cover data into one dataframe #####
  landCoverHa <- merge(merge(cropTableHa.2015
                             , marginSum %>% st_drop_geometry()
                             , by = "rcFid_1km", all = T)
                       , hedgePolyInter %>% dplyr::select(rcFid_1km, hedgeAreaHaperKm2) %>% st_drop_geometry()
                       , by = "rcFid_1km", all = T) %>%
    replace(is.na(.), 0)
  head(landCoverHa)
  
  # save
  st_write(landCoverHa, file.path(currPath, "landCover_margins.gpkg"), append = F)
  fwrite(st_drop_geometry(landCoverHa)
         , file.path(currPath, "landCover_margins.csv")
         , row.names = F)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(currPath, "readme_landCover_margins.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2023-05-03
Last update:",  format(Sys.Date()), "

The files containing 'landCoverHa' contain the area, in hectares, of different land covers, derived from CEH's LCM2015 and Land Cover® plus: Crops 2016 (Rowland et al., 2017), Natural England's Countryside Survey data, and Scholefield et al. (2016)
The columns:
      'rcFid' = 1 km reference id
      '[land cover type]' = area of the land covered the title land cover type in a 1 km2 pixel (ha)
          land covers:
           arable = winter wheat + winter barley + spring wheat + oats + maize + rapeseed + spring barley + potato + field beans + sugat beet
           improved_grass = improved grassland
           broadleaf = broadleaved woodland
           conifer = coniferous woodland
           neutral_grass = neutral grassland
           calc_grass = calcareous grassland
           acid_grass = acid grassland
           fen_marsh = fen, marsh, and swamp
           heather = heather
           heather_grass = heather grassland
           bog = bog
           inland_rock = inland rock
           saltwater = saltwater
           freshwater = freshwater
           supra_lit_rock = supra-littoral rock
           supra_lit_sed = supra-littoral sediment
           lit_rock = littoral rock
           lit_sed = littoral sediment
           saltmarsh = saltmarsh
           urban = urban
           suburban = suburban
           flowerMargin_ha = field margins (composed of flowers) 
           treeMargin_ha = field margins (composed of trees) 
           hedgeAreaHaperKm2 = hedgerow
           
units: ha/km2 
Spatial resolution: 1 km2
Spatial extent: England
Temporal coverage: 2016
Native projection: 27700
      
citations:  Rowland, C.S.; Morton, R.D.; Carrasco, L.; McShane, G.; O'Neil, A.W.; Wood, C.M. (2017). Land Cover Map 2015 (1km dominant aggregate class, GB). NERC Environmental Information Data Centre. https://doi.org/10.5285/711c8dc1-0f4e-42ad-a703-8b5d19c92247)
            Scholefield, P.A.; Morton, R.D.; Rowland, C.S.; Henrys, P.A.; Howard, D.C.; Norton, L.R. (2016). Woody linear features framework, Great Britain v.1.0. NERC Environmental Information Data Centre. https://doi.org/10.5285/d7da6cb9-104b-4dbc-b709-c1f7ba94fb16")
  sink(file = NULL)
  
  ##### 8e - adjust area based on margins #####
  # reduce area of arable and improved grassland for field margins and hedgerows respectively
  landCoverHaAdapt <- landCoverHa %>%
    mutate(arableNew_ha = Arable_ha - (flowerMargin_ha + treeMargin_ha)
           # if the above lead to negative arable, take off IG instead
           , IGnew_ha = ifelse(arableNew_ha < 0, improved_grass_ha - abs(arableNew_ha)
                               , improved_grass_ha)
           , arableNew_ha = ifelse(arableNew_ha < 0, 0, arableNew_ha)
           , IGnew_ha = ifelse(IGnew_ha < 0, 0, IGnew_ha)) %>%
    # reorder
    relocate(rcFid_1km, broadleaf_ha, conifer_ha, arableNew_ha, IGnew_ha) %>%
    # remove unneeded
    dplyr::select(-c(Arable_ha, improved_grass_ha))
  head(landCoverHaAdapt)
  
  ##### 8f - use 2015 data to get 2007 margins etc #####
  # get the average margins and hedges per pixel, based on land cover types
  landCoverHa.analysis <- landCoverHa %>% st_drop_geometry() %>%
    # calculate all grasses
    mutate(all_grass_ha = rowSums(select(., neutral_grass_ha:heather_grass_ha))) %>%
    # keep only grasses, arable, and improved grassland
    dplyr::select(rcFid_1km, Arable_ha, improved_grass_ha, all_grass_ha
                  , flowerMargin_ha: hedgeAreaHaperKm2)
  head(landCoverHa.analysis)
  
  ## model flower margin
  fmMod <- lm(flowerMargin_ha ~ Arable_ha + improved_grass_ha + all_grass_ha
              , data = landCoverHa.analysis)
  ## model tree margin
  tmMod <- lm(treeMargin_ha ~ Arable_ha + improved_grass_ha + all_grass_ha
              , data = landCoverHa.analysis)
  ## model hedges
  hedgeMod <- lm(hedgeAreaHaperKm2 ~ Arable_ha + improved_grass_ha + all_grass_ha
                 , data = landCoverHa.analysis)
  
  ### use the above models to determine the area of margins and hedges for the scenarios
  scen2007 <- list.files(file.path("scenario", "scen_maps", "land_cover")
                         , pattern = ".gpkg$", full.names = T)
  landCoverHa.scens <- pblapply(scen2007, function(i) {
    # read in
    xInSpat <- st_read(i)
    xIn <- xInSpat %>% st_drop_geometry() %>%
      # calculate all grasses
      mutate(all_grass_ha = rowSums(select(., neutral_grass_ha:heather_grass_ha))) %>%
      # keep only grasses, arable, and improved grassland
      dplyr::select(rcFid_1km, Arable_ha, improved_grass_ha, all_grass_ha) %>%
      ## predict the margins and hedges
      mutate(hedgeAreaHaperKm2 = predict(hedgeMod, .)
             , flowerMargin_ha = predict(fmMod, .)
             , treeMargin_ha = predict(tmMod, .)) %>%
      ### if any below 0, make 0
      mutate(hedgeAreaHaperKm2 = ifelse(hedgeAreaHaperKm2<0, 0, hedgeAreaHaperKm2)
             , flowerMargin_ha = ifelse(flowerMargin_ha<0, 0, flowerMargin_ha)
             , treeMargin_ha = ifelse(treeMargin_ha<0, 0, treeMargin_ha))
    head(xIn)
    # get nm
    nm <- gsub(".*area(.+).gpkg", "\\1", i)
    print(nm)
    # range
    print(range(xIn$hedgeAreaHaperKm2))
    
    ### adjust area based on margins and hedges
    # reduce area of arable and improved grassland for field margins and hedgerows respectively
    landCoverHaAdapt <- xIn %>%
      dplyr::select(rcFid_1km, hedgeAreaHaperKm2 : treeMargin_ha) %>%
      # merge with the original spatial df
      merge(xInSpat, .) %>%
      mutate(arableNew_ha = Arable_ha - (flowerMargin_ha + treeMargin_ha)
             # if the above lead to negative arable, take off IG instead
             , IGnew_ha = ifelse(arableNew_ha < 0, improved_grass_ha - abs(arableNew_ha)
                                 , improved_grass_ha)
             , arableNew_ha = ifelse(arableNew_ha < 0, 0, arableNew_ha)
             , IGnew_ha = ifelse(IGnew_ha < 0, 0, IGnew_ha)) %>%
      # reorder
      relocate(rcFid_1km, broadleaf_ha, conifer_ha, arableNew_ha, IGnew_ha) %>%
      # remove unneeded
      dplyr::select(-c(Arable_ha, improved_grass_ha))
    head(landCoverHaAdapt)
    return(landCoverHaAdapt)
  })
  head(landCoverHa.scens[[1]])
  
  ##### 8g - get emissions (for original data and scenarios) #####
  # list all data
  ## original = 1 + scens
  landCoverHa.length <- 1 + length(landCoverHa.scens)
  
  # read in table with CO2e coefficients
  lcCoefficients <- fread("data_in/land_cover/lcm_land_cover_co2.csv") %>% as.data.frame() %>%
    rename(tCO2_ha = 3) %>%
    # get value per 25 m2
    mutate(tCO2_25m2 = tCO2_ha / 10000)
  head(lcCoefficients)
  
  for(i in 1:landCoverHa.length){
    
    if(i == 1){
      lcIn <- landCoverHa
      nm <- "original_data"
      saveName <- file.path(savePath, "land_use_emissions.gpkg")
    } else {
      lcIn <- landCoverHa.scens[[i-1]]
      nm <- gsub("_land.gpkg", "", basename(landArea.names[[i-1]]))
      saveName <- file.path("scenario", "results"
                            , paste0("land_use_emissions", nm, ".gpkg"))
    }
    cat("... Calculating emissions from land use for", nm, "... \n")
    # print(head(lcIn))
    
    # ensure the land use headers, and coefficients are in the same order
    stopifnot(grepl("rable", lcCoefficients$`LCM2015 target class`[[3]]) & grepl("rable", names(lcIn)[4]))
    stopifnot(grepl("uburban", lcCoefficients$`LCM2015 target class`[[21]]) & grepl("uburban", names(lcIn)[22]))
    
    ##### non-margin land use #####
    # multiply the area (in ha) by the per ha values
    cropTimesCoef <- lcIn %>% st_drop_geometry() %>%
      dplyr::select(-c(rcFid_1km, flowerMargin_ha, treeMargin_ha, hedgeAreaHaperKm2)) %>%
      as.matrix %*% diag(lcCoefficients$tCO2_ha) %>% as.data.frame()
    # rename
    names(cropTimesCoef) <- paste0(names(lcIn)[2:22], "_tCO2")
    names(cropTimesCoef) <- sub("_ha", "", (names(cropTimesCoef)))
    head(cropTimesCoef)
    # get totals
    cropTimesCoef <- cropTimesCoef %>%
      mutate(total_landuse_tco2 = rowSums(select(., c(broadleaf_tCO2:suburban_tCO2))))
    
    # add ref id back in
    cropTimesCoef <- bind_cols(rcFid_1km = lcIn$rcFid_1km
                               , cropTimesCoef) 
    stopifnot(names(cropTimesCoef)[1] == "rcFid_1km")
    head(cropTimesCoef)
    
    ##### field margins #####    
    # get co2e emissions/ uptake from size of field margins - grass strips
    # using the middle value from Yang et al. (2019),
    #  = 1.834 Mg CO2 ha-1 yr-1
    # note: lots of different estimates exist - from uptake to emissions
    # multiply by area
    
    # get co2e emissions/ uptake from size of field margins - tree strips
    # use the value from Falloon et al. (2004) for natural BL regen:
    # 10276 kg co2 (see Google Doc)
    
    landCoverHaMargin <- lcIn %>%
      mutate(fm_kgCo2Ha = (
        # converts Mg to kg, and uptake
        (lcIn$flowerMargin_ha * 1.834) / 1000 * -1) +  # if flower
          (lcIn$treeMargin_ha * 10276) * -1) # if tree
    head(landCoverHaMargin)
    
    ##### hedgerow #####   
    # using the values from Blair (2021), convert ground area into emissions
    landCoverHaHedge <- lcIn %>%
      mutate(emisHedge_tco2e = hedgeAreaHaperKm2 * -9.305)
    
    ##### combine ##### 
    combineEmis <- Reduce(function(x, y) {
      merged_data <- merge(
        x, y
        , all = TRUE
      )
      return(merged_data)
    }, list(cropTimesCoef, landCoverHaMargin, landCoverHaHedge)) %>%
      dplyr::select(-c(broadleaf_ha: suburban_ha, hedgeAreaHaperKm2)) %>%
      # create t for kg for fm_kgCo2Ha
      mutate(fm_tco2 = fm_kgCo2Ha/1000) %>%
      st_drop_geometry() %>%
      dplyr::select(-c(fm_kgCo2Ha, flowerMargin_ha, treeMargin_ha)) %>%
      # final sum
      mutate('landCoverTotalEmis_tco2' = rowSums(select(., c(total_landuse_tco2, emisHedge_tco2e, fm_tco2)), na.rm = T)) %>%
      st_as_sf()
    # head(combineEmis)
    # str(combineEmis)
    attributes(combineEmis$total_landuse_tco2) <- NULL
    
    cat("           Saving here:", saveName, "\n")
    st_write(combineEmis, saveName, append = F)
    fwrite(st_drop_geometry(combineEmis)
           , gsub(".gpkg", ".csv", saveName)
           , row.names = F)
  }
  
  # write readme (for original data)
  sink(file = NULL)
  sink(file = file.path(savePath, "readme_land_use_emissions.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2023-05-03
Last update:",  format(Sys.Date()), "

The files containing 'land_use_emissions' contain information the annual amount of emissions or uptake CO2 per km2, based on underlying land covers. The land cover areas were derived from the 2015 land cover map, countryside survey (field margins), and woody linear features

The columns:
      'rcFid_1km' = 1 km reference id
      '[land cover type]_tCO2' = emissions or sequestration (represented by a '-' value) of CO2 that each land cover contributes across a 1 km2 (units: tonnes CO2/km2/yr)
      'total_landuse_tco2' = total emissions or sequestration (represented by a '-' value) of CO2 for a 1 km2, based on the composition of land covers at 25 m2 (units: tonnes CO2/km2/yr)
      'emisHedge_tco2e' = total emissions or sequestration (represented by a '-' value) of CO2 for a 1 km2, based on area of hedges (units: tonnes CO2/km2/yr)
      'fm_tco2' = total emissions or sequestration (represented by a '-' value) of CO2 for a 1 km2, based on area of field margins (units: tonnes CO2/km2/yr)
      'landCoverTotalEmis_tco2' = total emissions or sequestration (represented by a '-' value) of CO2 for a 1 km2, based on broad land cover categories, field margins, and hedges (units: tonnes CO2/km2/yr)
Data:
    units: tonnes CO2/km2/yr 
    Spatial resolution: 1 km2
    Spatial extent: GB
    Temporal coverage: 2015
    Native projection: 27700")
  
  sink(file = NULL)
  
  # write readme (for scenarios)
  sink(file = NULL)
  sink(file = file.path("scenario", "results", "readme_land_use_emissions.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2023-05-03
Last update:",  format(Sys.Date()), "

The files containing 'land_use_emissions' contain information the annual amount of emissions or uptake CO2 per km2, based on underlying land covers. The land cover areas were derived from the 2007 land cover map, and the scenarios found in Redhead et al. (2020).
The suffixes of the file names relate to specific scenarios taken from Redhead et al. (2020):
  'Ag_15' and 'Ag_30' = 15 and 30% agricultural expansion from the 2007 baseline
  'Gr_15' and 'Gr_30' = 15 and 30% grassland restoration from the 2007 baseline
  'baseline2007' = baseline land cover from 2007
  
The columns:
      'rcFid_1km' = 1 km reference id
      '[land cover type]_tCO2' = emissions or sequestration (represented by a '-' value) of CO2 that each land cover contributes across a 1 km2 (units: tonnes CO2/km2/yr)
      'total_landuse_tco2' = total emissions or sequestration (represented by a '-' value) of CO2 for a 1 km2, based on the composition of land covers at 25 m2 (units: tonnes CO2/km2/yr)
      'emisHedge_tco2e' = total emissions or sequestration (represented by a '-' value) of CO2 for a 1 km2, based on area of hedges (units: tonnes CO2/km2/yr)
      'fm_tco2' = total emissions or sequestration (represented by a '-' value) of CO2 for a 1 km2, based on area of field margins (units: tonnes CO2/km2/yr)
      'landCoverTotalEmis_tco2' = total emissions or sequestration (represented by a '-' value) of CO2 for a 1 km2, based on broad land cover categories, field margins, and hedges (units: tonnes CO2/km2/yr)
Data:
    units: tonnes CO2/km2/yr 
    Spatial resolution: 1 km2
    Spatial extent: GB
    Temporal coverage: 2015
    Native projection: 27700")
  
  sink(file = NULL)
}

#### 9 - part9residue ####
if(part9residue){
  
  # stop("residue time - part 9")
  
  ## ------------ Notes --------------  ##
  ## The residue amounts come from WOFOST model output for the most part
  ## (with the exception of oats and field beans)
  
  ## There are six crop residue management options:
  ## Removed; left untreated in heaps or pits
  ## Removed; non-forced-aeration compost
  ## Removed; forced-aeration compost
  ## Left on field; incorporated or mulch
  ## Burned (no longer done in the UK)
  ## Exported off farm
  ## ------------ ----- --------------  ##
  
  # set input location - for land cover
  currPath <- file.path("./data_in", "land_cover")
  # results
  savePath <- "./results/arable"
  
  residuePath <- "N:/Projects/AgLand/models/crop_production/potential_production/combined_outputs"
  # create path to save the means - yearly
  dir.create(meanResPath <- file.path(residuePath, "means_crop")
             , showWarnings = F)
  # create path to save the means - variety
  dir.create(resPathBroadCrop <- file.path(residuePath, "means_crop", "broad_crop")
             , showWarnings = F)
  CropPath <- "./data_in/crops"
  
  # import 1 km2 crop area for 2015
  landArea <- st_read(file.path(currPath, "land_cover_table.gpkg"))
  head(landArea)
  
  ##### 9a - create residue functions #####  
  # create function to produce yearly means from input residue data
  # from each unique crop variety, get an average value for the five-year span the
  # data came from (2012 - 2016)
  aveCropYearFunc <- function(x){
    
    if(file.exists(file.path(meanResPath, paste0(x, ".csv")))){
      cat("already exists:", file.path(meanResPath, paste0(x, ".csv")), "\n")
    } else {
      
      # select all files with that crop variety identifier
      xResList <- resiList[grepl(x, resiList)]
      # check it is five
      stopifnot(length(xResList) == 5)
      
      # read in
      varList <- pblapply(xResList, function (y) { fread(y) %>%
          as.data.frame() %>%
          dplyr::select(c(pixsn, X, Y
                          , TAGP,  TWSO, TWLV, TWST, TWRT))})
      
      # get the different columns
      # TAGP
      TAGPmean <- sapply(varList, function(x) x$TAGP) %>%
        as.data.frame() %>%
        mutate(mean = rowMeans(.[]))
      # TWSO
      TWSOmean <- sapply(varList, function(x) x$TWSO) %>%
        as.data.frame() %>%
        mutate(mean = rowMeans(.[]))
      # TWLV
      TWLVmean <- sapply(varList, function(x) x$TWLV) %>%
        as.data.frame() %>%
        mutate(mean = rowMeans(.[]))
      # TWST
      TWSTmean <- sapply(varList, function(x) x$TWST) %>%
        as.data.frame() %>%
        mutate(mean = rowMeans(.[]))
      # TWRT
      TWRTmean <- sapply(varList, function(x) x$TWRT) %>%
        as.data.frame() %>%
        mutate(mean = rowMeans(.[]))
      
      # combine back into one df
      meanYearDf <- varList[[1]] %>%
        dplyr::select(pixsn, X, Y) %>%
        # bind to obtain yearly means
        bind_cols(., list(TAGPmean$mean
                          , TWSOmean$mean
                          , TWLVmean$mean
                          , TWSTmean$mean
                          , TWRTmean$mean)) %>%
        rename(TAGP = 4
               , TWSO = 5
               , TWLV = 6
               , TWST = 7
               , TWRT = 8)
      
      # save to file
      fwrite(meanYearDf
             , file.path(meanResPath, paste0(x, ".csv"))
             , row.names = F)
      cat("saved at:", file.path(meanResPath, paste0(x, ".csv")), "\n")
    }
  }
  
  # create function to produce all-variety means from input average annual residue data
  # from each unique broad crop type, get the average. This will encompass different varieties
  # for each broad crop type
  aveCropVarFunc <- function(x){
    
    if(file.exists(file.path(resPathBroadCrop, paste0(x, "_residues.csv")))){
      cat("already exists:", file.path(resPathBroadCrop, paste0(x, "_residues.csv")), "\n")
    } else {
      
      # select all files with that crop variety identifier
      xResList <- resiList[grepl(x, resiList)]
      # check names
      cat(xResList, sep = "\n")
      
      # read in
      varList <- pblapply(xResList, function (y) { fread(y) %>%
          as.data.frame() %>%
          dplyr::select(c(pixsn, X, Y
                          , TAGP,  TWSO, TWLV, TWST, TWRT))})
      
      # get the different columns
      # TAGP
      TAGPmean <- sapply(varList, function(x) x$TAGP) %>%
        as.data.frame() %>%
        mutate(mean = rowMeans(.[]))
      # TWSO
      TWSOmean <- sapply(varList, function(x) x$TWSO) %>%
        as.data.frame() %>%
        mutate(mean = rowMeans(.[]))
      # TWLV
      TWLVmean <- sapply(varList, function(x) x$TWLV) %>%
        as.data.frame() %>%
        mutate(mean = rowMeans(.[]))
      # TWST
      TWSTmean <- sapply(varList, function(x) x$TWST) %>%
        as.data.frame() %>%
        mutate(mean = rowMeans(.[]))
      # TWRT
      TWRTmean <- sapply(varList, function(x) x$TWRT) %>%
        as.data.frame() %>%
        mutate(mean = rowMeans(.[]))
      
      # combine back into one df
      meanVarietyDf <- varList[[1]] %>%
        dplyr::select(pixsn, X, Y) %>%
        # bind to obtain yearly means
        bind_cols(., list(TAGPmean$mean
                          , TWSOmean$mean
                          , TWLVmean$mean
                          , TWSTmean$mean
                          , TWRTmean$mean)) %>%
        rename(TAGP = 4
               , TWSO = 5
               , TWLV = 6
               , TWST = 7
               , TWRT = 8) %>%
        # calculate total weight of residues for a pixel
        mutate(total_res_kgha = TWLV + TWST + TWRT)
      
      # save to file
      fwrite(meanVarietyDf
             , file.path(resPathBroadCrop, paste0(x, "_residues.csv"))
             , row.names = F)
      cat("saved at:", file.path(resPathBroadCrop, paste0(x, "_residues.csv")), "\n")
    }
  }
  
  #### 9b - load in crop values ####
  # list all crop tables available
  resiList <- list.files(residuePath
                         , pattern = ".csv$"
                         , full.names = T)
  cat(resiList, sep = "\n")
  
  # determine all unique crop variety types
  cropVarUnique <- unique(gsub("^(.+)_2.*", "\\1", basename(resiList)))
  cat(cropVarUnique, sep = "\n")
  
  # create yearly averages using function
  pblapply(cropVarUnique, aveCropYearFunc)
  
  # write readme
  sink(file = file.path(meanResPath, "readme.txt"))
  cat("Description of csv files
  -----------------------------------------------------------
 
  Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-05-03
  Last update:",  format(Sys.Date()), "

  Files in this directory contain yearly average properties of a certain crop variety type for a certain pixel. 
  Five years were used to get the yearly average: 2012 - 2016, with each individual year initially providing a single value for each property  
  The names of each file indicates the broad crop type and its variety number, in the form '[cropType]_[VarietyType].csv' 
      where:
        - [cropType] indicates the broad crop type
        - [VarietyType] indicates a crop's variety number

  The columns:
        'pixsn' = pixel identifier (running number)
        'X' = X-coordinate of 1 km pixel centre (EPSG: 27700)
        'Y' = Y-coordinate of 1 km pixel centre (EPSG: 27700)
        'TAGP' = Total above ground production (dry weight (Kg ha-1))
        'TWSO' = Total storage organ weight (dry weight (Kg ha-1))
        'TWLV' = Total weight of leaves (dry weight (Kg ha-1))
        'TWST' = Total weight of stems (dry weight (Kg ha-1))
        'TWRT' = Total weight of roots (dry weight (Kg ha-1))
        
  Data:
      Units: dry weight (Kg ha-1) (for the non-coordinate columns)
      Spatial resolution: 1 km2
      Spatial extent: GB
      Temporal coverage: 2012 - 2016
      Native projection: 27700")
  
  sink(file = NULL)
  
  #### 9c - get average for a particular crop type per pixel, per ha ####
  
  # list all crop tables available - these have already been averaged across five years
  resiList <- list.files(meanResPath
                         , pattern = ".csv$"
                         , full.names = T)
  
  # determine all unique broad crop types
  cropVarUnique <- unique(gsub("^(.+)_.*", "\\1", basename(resiList)))
  cat(cropVarUnique, sep = "\n")
  
  # create yearly averages using function
  pblapply(cropVarUnique, aveCropVarFunc)
  
  # write readme
  sink(file = file.path(resPathBroadCrop, "readme.txt"))
  cat("Description of csv files
  -----------------------------------------------------------
 
  Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
  Date created: 2023-05-03
  Last update:",  format(Sys.Date()), "

  Files in this directory contain average broad crop properties for a certain pixel. 
  Different number of crop varieties were used to get the average for each property  
  The names of each file indicates the broad crop type and its amount of material/residue, in the form '[cropType]_residues.csv' 
      where:
        - [cropType] indicates the broad crop type

  The columns:
        'pixsn' = pixel identifier (running number)
        'X' = X-coordinate of 1 km pixel centre (EPSG: 27700)
        'Y' = Y-coordinate of 1 km pixel centre (EPSG: 27700)
        'TAGP' = Total above ground production (dry weight (Kg ha-1))
        'TWSO' = Total storage organ weight (dry weight (Kg ha-1))
        'TWLV' = Total weight of leaves (dry weight (Kg ha-1))
        'TWST' = Total weight of stems (dry weight (Kg ha-1))
        'TWRT' = Total weight of roots (dry weight (Kg ha-1))
        'total_res_kgha' = total weight of residues (i.e. all material minus the storage organ weight) (dry weight (Kg ha-1))
        
  Data:
      Units: dry weight (Kg ha-1) (for the non-coordinate columns)
      Spatial resolution: 1 km2
      Spatial extent: GB
      Temporal coverage: 2012 - 2016
      Native projection: 27700")
  
  sink(file = NULL)
  
  #### load residue/crop parameters, and residue amounts ####
  ## ------------ Notes --------------  ##
  ## the residue parameters come from the CFT values currently
  ## ------------ ----- --------------  ##
  
  residParam <- read_excel(file.path(CropPath, "residue_parameters.xlsx"), sheet = "resid_data")
  head(residParam)
  cropParam <- read_excel(file.path(CropPath, "residue_parameters.xlsx"), sheet = "crop_data")
  head(cropParam)
  
  #### 9d - use area to determine total residue in a 1 km2 ####
  ## ------------ Notes --------------  ##
  ## do this for both original data and the scenarios
  ## ------------ ----- --------------  ##
  # list all data
  ## original = 1 + scens
  landArea.length <- 1 + length(scenarios)
  ## get scen names
  sNames <- gsub(".gpkg", "", basename(scenarios))
  ### and use them to load appropriate land cover
  landCoverHa.scens <- list.files(file.path("scenario", "scen_maps", "finer_detail")
                                  , pattern = ".*_land.gpkg"
                                  , full.names = T)
  landCoverHa.scens
  
  for(iscen in 1:landArea.length){
    
    if(iscen == 1){
      lcIn <- landArea
      nm <- "original_data"
      saveName <- file.path(savePath, "residEmisComb_kgCO2e.gpkg")
    } else {
      lcIn <- st_read(landCoverHa.scens[[iscen-1]])
      nm <- gsub("_land.gpkg", "", basename(landCoverHa.scens[[iscen-1]]))
      saveName <- file.path("scenario", "results"
                            , paste0("residEmisComb_kgCO2e", nm, ".gpkg"))
    }
    cat("... Calculating emissions from residues for", nm, "... \n")
    
    # check row sums
    cropTableSum <- lcIn %>%
      mutate(crop_area = rowSums(dplyr::select(.[,,drop = T], c(winterwheat_ha:sugarbeet_ha)), na.rm = T))
    head(cropTableSum)
    
    # select just crops
    arableSelection <- lcIn  %>%
      st_drop_geometry() %>%
      dplyr::select(c(winterwheat_ha:sugarbeet_ha))
    
    # list all residue tables available - these have already been averaged 
    resiList <- list.files(resPathBroadCrop
                           , pattern = ".csv$"
                           , full.names = T)
    # cat(resiList, sep = "\n")
    
    # loop through the files, extracting only the total residue amount
    # resave, in a new list
    cat("   ... Extracting residue amount", "\n")
    newResiList <- list()
    for(i in 1:length(resiList)){
      
      # get crop name
      cropNm <- gsub("^(.+)_r.*", "\\1", basename(resiList[[i]]))
      
      newResiList[[i]] <- fread(resiList[[i]]) %>% as.data.frame() %>%
        dplyr::select(c(pixsn, X, Y
                        , total_res_kgha)) %>%
        rename(!!cropNm := total_res_kgha)
    }
    
    cat("   ... Merging dataframes", "\n")
    # merge all the final residue dataframes together
    allResidues <- purrr::reduce(newResiList, merge)
    # rename, to indicate residues
    names(allResidues)[4:ncol(allResidues)] <- paste0(names(allResidues)[4:ncol(allResidues)], "_resKgHa")
    allResidues <- allResidues %>%
      mutate(X = as.numeric(X),
             Y = as.numeric(Y)) %>%
      mutate(rcFid_1km = paste(round(X, 0), round(Y, 0), sep = "_")) %>%
      dplyr::select(-c(X, Y, pixsn))
    head(allResidues)
    
    # merge residues with crop composition 
    ## get centroid of crop table
    cropCentroid <- lcIn %>%
      st_centroid() 
    head(cropCentroid)
    
    ## merge the centroid point crop df with residues df
    cropAndResidues <- merge(cropCentroid
                             , allResidues
                             , by = "rcFid_1km"
                             , all = T) %>%
      st_drop_geometry()
    head(cropAndResidues)
    class(cropAndResidues)
    
    # get the total amount of residues per pixel
    # based on area of crop, and amount of average residues from the crop model
    names(cropAndResidues)
    
    ## get the positions of the columns
    startCol <- which(names(cropAndResidues) == "winterwheat_ha")
    endCol <- which(names(cropAndResidues) == "sugarbeet_ha")
    
    ## find improved grassland, so it is not included
    igCol <- which(names(cropAndResidues) == "improved_grass_ha")
    
    ### create a list, ignoring IG
    colList <- seq(startCol, endCol, 1)
    colList <- colList[colList != igCol]
    
    cat("   ... Listing crops and residues", "\n")
    # create empty list
    resDfList <- list()
    
    i.count <- 0
    for(i in colList){
      i.count <- i.count + 1
      
      # remove "_ha", to get name
      readyName <- sub("_ha", "", names(cropAndResidues)[[i]])
      rdName <- paste0(readyName, "_totalResid")
      
      cat(paste0("iscen = ", iscen)
          , paste0("i.count = ", i.count)
          , paste0("crop column = ", i)
          , paste0("crop = ", readyName, "\n")
          , sep = " | ")
      
      # use the column name to get residue and area
      resDf <- cropAndResidues[, grepl(readyName, names(cropAndResidues))] %>%
        as.data.frame() %>% st_drop_geometry() %>%
        ## make NAs 0
        replace(is.na(.), 0)
      head(resDf)
      if(ncol(resDf)>1){
        cat(names(resDf), sep = " | ")
        cat("\n")
      }
      
      # only continue if there are two cols: crop area and residue
      if(ncol(resDf) == 2){
        
        cat(readyName, "| extracted; OK\n")
        
        # multiply the two columns
        resDfList[[i.count]] <- apply(resDf, 1, prod, na.rm=TRUE) %>%
          as.data.frame() %>% rename(!!rdName := 1)
        head(resDf)
        head(resDfList[[i.count]])
        xMax <- which(resDfList[[i.count]][, 1] == max(resDfList[[i.count]][, 1]))
        
        stopifnot(resDfList[[i.count]][10, 1] == resDf[10, 1] * resDf[10, 2])
        stopifnot(resDfList[[i.count]][100, 1] == resDf[100, 1] * resDf[100, 2])
        stopifnot(resDfList[[i.count]][10000, 1] == resDf[10000, 1] * resDf[10000, 2])
        stopifnot(resDfList[[i.count]][xMax, 1] == resDf[xMax, 1] * resDf[xMax, 2])
        
      } else {
        
        # if not from the crop model, equation from Table 2.1 (Eq. 2.2.2) in CFT is required
        # yield for the equation comes from https://ahdb.org.uk/news/top-of-the-crops - the top 25%
        
        cat(readyName, "| needed another value; OK\n")
        
        if(readyName == "oats"){
          
          # rResiduesKgHa Total amount of residues [kg]
          # rBelow Amount of below-ground residues [kg]
          # rAbove Amount of above-ground residues [kg]
          
          rAbove <- 7.5 * 0.89 * 0.91 + 0.89
          rBelow <- rAbove * 0.25
          rResiduesKgHa <- rAbove + rBelow
          
          # combine area and average residue amount
          resDf <- bind_cols(resDf, rResiduesKgHa)
          # multiply the two columns
          resDfList[[i.count]] <- apply(resDf, 1, prod, na.rm=TRUE) %>%
            as.data.frame() %>% rename(!!rdName := 1)
          
        }
        
        if(readyName == "fieldbeans"){
          
          # rResiduesKgHa Total amount of residues [kg]
          # rBelow Amount of below-ground residues [kg]
          # rAbove Amount of above-ground residues [kg]
          
          rAbove <- 4.4 * 0.9 * 0.36 + 0.68
          rBelow <- rAbove * 0.19
          rResiduesKgHa <- rAbove + rBelow
          
          # combine area and average residue amount
          resDf <- bind_cols(resDf, rResiduesKgHa)
          # multiply the two columns
          resDfList[[i.count]] <- apply(resDf, 1, prod, na.rm=TRUE) %>%
            as.data.frame() %>% rename(!!rdName := 1)
        }
      }
    }
    
    # if(iscen == 2){
    #   stop("iscensore")
    # }
    str(resDfList)
    
    cat("   ... Combining the result", "\n")
    # combine the outcome - after testing to ensure that crop and outcome match
    testy <- bind_cols(resDfList) %>% as.data.frame()
    # for(tt in 1:ncol(testy)){
    #   # ensure the highest places match
    #   tt3 <- which(testy[, tt] == max(testy[, tt]))
    #   tt4 <- which(cropAndResidues[, colList[[tt]]] == max(cropAndResidues[, colList[[tt]]], na.rm = T))
    #   stopifnot(identical(tt3, tt4))
    #   }
    stopifnot(identical(which(cropAndResidues[10000, colList] > 0)
                        , which(testy[10000, ] > 0)))
    stopifnot(identical(which(cropAndResidues[100000, colList] > 0)
                        , which(testy[100000, ] > 0)))
    stopifnot(identical(which(cropAndResidues[500000, colList] > 0)
                        , which(testy[500000, ] > 0)))
    
    ## merge all the final residue dataframes together
    allResidueskgKm2 <- bind_cols(resDfList) %>% as.data.frame() %>%
      bind_cols(cropAndResidues %>% dplyr::select(c(rcFid_1km)), .)
    names(allResidueskgKm2) <- sub("_ha", "", names(allResidueskgKm2))
    head(allResidueskgKm2)
    head(cropAndResidues[, c(1, colList)])
    
    ## ------------ Notes --------------  ##
    ## the residue amounts were calculated
    ## ------------ ----- --------------  ##
    
    # amounts
    residAmounts <- allResidueskgKm2 
    
    ## get aboveground residues by using the ratio of below and above ground (i.e.'ratio_blw:abg')
    ## this is due to only aboveground residues 
    ## loop through crops individually
    rAbove <- list()
    rBelow <- list()
    
    nContentinResid <- list()
    
    for(i in 2:ncol(residAmounts)){
      
      # extract column and name
      newCol <- residAmounts[, i] %>% as.data.frame()
      newNm <- gsub("_totalResid", "", names(residAmounts)[i])
      # use name to get parameters
      newParas <- cropParam %>%
        filter(grepl(newNm, Crop)) %>%
        dplyr::select("ratio_blw:abg")
      
      # multiply column by amount
      rAbove[[i]] <- newCol * as.numeric(1 - newParas)
      rBelow[[i]] <- newCol * as.numeric(newParas)
      
      # get N content of all residues (required if mulching)
      # use name for N parameters
      # see Eq 2.2.6 in CFT
      nContentParas <- cropParam %>%
        filter(grepl(newNm, Crop)) %>%
        dplyr::select(c(Abg_N_frac, Blw_N_frac))
      nContentAbove <- rAbove[[i]] * as.numeric(nContentParas[1])
      nContentBelow <- rBelow[[i]] * as.numeric(nContentParas[2])
      nContentinResid[[i]] <- nContentAbove + nContentBelow
      print(dim(nContentinResid[[i]]))
    }
    
    # merge all the final residue dataframes together
    ## above
    rAbovekgKm2 <- do.call(cbind, rAbove[2:11]) %>% as.data.frame() %>%
      bind_cols(residAmounts %>% dplyr::select(c(rcFid_1km)), .)
    ## below
    rBelowkgKm2 <- do.call(cbind, rBelow[2:11]) %>% as.data.frame() %>%
      bind_cols(residAmounts %>% dplyr::select(c(rcFid_1km)), .)
    # change names
    names(rAbovekgKm2) <- names(rBelowkgKm2) <- names(residAmounts)
    head(rAbovekgKm2)
    
    # merge N content residues together
    nContentinResid <- do.call(cbind, nContentinResid[2:11]) %>% as.data.frame() %>%
      bind_cols(residAmounts %>% dplyr::select(c(rcFid_1km)), .) %>%
      replace(is.na(.), 0)
    names(nContentinResid) <- gsub("_totalResid", "_residN_kg", names(residAmounts))
    head(nContentinResid)
    
    #### 9f - calculate emissions based on different management ####
    # use the rAbove values and multiply them by management type coefficient
    ## get total km2 residue amounts
    rAbovekgKm2Totals <- rAbovekgKm2 %>%
      mutate(totalResidkgKm = rowSums(.[, 2:11], na.rm = T))
    head(rAbovekgKm2Totals)
    
    # create empty list to store different emissions
    residEmisCH4 <- list()
    residEmisN2O <- list()
    residEmisCombCO2e <- list()
    
    # loop through the possible treatments, ignoring burning (which is no longer done in the UK)
    for(i in 1:nrow(residParam)){
      # get name
      nm <- residParam$treatment[[i]]
      print(nm)
      
      # get CH4 and N2O coefs
      ch4 <- as.numeric(residParam[i, "ch4"])
      n2o <- as.numeric(residParam[i, "n2o"])
      
      residEmisCH4[[i]] <- rAbovekgKm2Totals$totalResidkgKm * ch4
      residEmisN2O[[i]] <- rAbovekgKm2Totals$totalResidkgKm * n2o
      residEmisCombCO2e[[i]] <- residEmisCH4[[i]] + residEmisN2O[[i]]
    }
    
    # merge all the final residue dataframes together
    ## ch4 (already converted to CO2e)
    residEmisCH4 <- do.call(cbind, residEmisCH4) %>% as.data.frame() %>%
      bind_cols(rAbovekgKm2Totals %>% dplyr::select(c(rcFid_1km, totalResidkgKm)), .)
    ## n2o (already converted to CO2e)
    residEmisN2O <- do.call(cbind, residEmisN2O) %>% as.data.frame() %>%
      bind_cols(rAbovekgKm2Totals %>% dplyr::select(c(rcFid_1km, totalResidkgKm)), .)
    ## combined n2o and ch4 (already converted to CO2e)
    residEmisCombCO2e <- do.call(cbind, residEmisCombCO2e) %>% as.data.frame() %>%
      bind_cols(rAbovekgKm2Totals %>% dplyr::select(c(rcFid_1km, totalResidkgKm)), .)
    # change names
    names(residEmisCH4)[3:8] <- names(residEmisN2O)[3:8] <- names(residEmisCombCO2e)[3:8] <- paste0(residParam$treatment, " (100%)")
    head(residEmisCH4)
    
    # save
    tic("written to file")
    cat("   ... writing to file ...\n")
    fwrite(residEmisCombCO2e
           , gsub(".gpkg", ".csv", saveName)
           , row.names = F)
    cat("   ...... here:", saveName, "\n")
    # combine back to spatial
    residEmisCombCO2eSpat <- landArea %>% as.data.frame() %>%
      dplyr::select(rcFid_1km) %>%
      merge(., residEmisCombCO2e, by = "rcFid_1km", all = T)
    st_write(residEmisCombCO2eSpat
             , saveName
             , append = F)
    toc()
  }
  
  # write readme
  sink(file = file.path(savePath, "readme_residEmisComb_kgCO2e.txt"))
  cat("Description of 'residEmisCombCO2e' files
Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2023-05-03
Last update:",  format(Sys.Date()), "

Files in this directory containing 'residEmisComb_kgCO2e' show emissions (in kg CO2e) produced from different residue managements.

The columns:
      'rcFid_1km' = 1 km2 pixel identifier. Coordinate of 1 km pixel centre (EPSG: 27700)
      '[managementType] (100%)' indicates the emissions produced, in kg CO2e, when all residues within a km2 are all treated the same way, with the management type indicated in the column header
        
Data:
    Units: kg CO2e (produced by residue management)
    Spatial resolution: 1 km2
    Spatial extent: GB
    Temporal coverage: 2015/2016
    Native projection: 27700")
  sink(file = NULL)
  
  # write readme - for scenarios
  sink(file = file.path("scenario", "results", "readme_residEmisComb_kgCO2e.txt"))
  cat("Description of 'residEmisCombCO2e' files
Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2023-05-03
Last update:",  format(Sys.Date()), "

Files in this directory containing 'residEmisComb_kgCO2e' show emissions (in kg CO2e) produced from different residue managements.

The columns:
      'rcFid_1km' = 1 km2 pixel identifier. Coordinate of 1 km pixel centre (EPSG: 27700)
      '[managementType] (100%)' indicates the emissions produced, in kg CO2e, when all residues within a km2 are all treated the same way, with the management type indicated in the column header
        
Data:
    Units: kg CO2e (produced by residue management)
    Spatial resolution: 1 km2
    Spatial extent: GB
    Temporal coverage: 2015/2016
    Native projection: 27700")
  sink(file = NULL)
}

#### 10 - part10onfarmenergy ####
if(part10onfarmenergy){
  
  # stop("start of energy - part 10")
  
  currPath <- file.path("data_in", "land_cover")
  savePath <- file.path("results", "land_use")
  
  ##### 10a - on-farm fuel use #####
  
  # load in fuel use
  fuelUse <- fread("data_in/energy/fuel_use.csv") %>%
    # remove first row
    slice(-1) %>%
    # remove columns that being with 'v' - these are the confidence intervals
    dplyr::select(-c(starts_with('v'), Observations)) %>%
    # rename for clarity
    rename(road_fuel_l = "Road fuel (l)*"
           , red_diesel_l = "Red diesel (l)"
           , red_diesel_con_l = "Red diesel used by contractors (l)"
           , lpg_kg = "LPG (kg)"
           , kerosene_l = "Kerosene (l)"
           , elec_units = "Electricity (units)"
           , heating_oil_l = "Heating oil (l)") %>%
    # convert amounts to numeric
    mutate(across(2:8, ~as.numeric(.))) %>%
    replace(is.na(.), 0) %>%
    # sum red diesel
    mutate(red_dies_l = red_diesel_l + red_diesel_con_l) %>%
    # add heating oil to kerosene
    mutate(kerosene_l = kerosene_l + heating_oil_l) %>% 
    dplyr::select(-c(red_diesel_l, red_diesel_con_l, heating_oil_l)) %>% as.data.frame()
  head(fuelUse)
  
  # get a list of all the relevant fuel types
  names(fuelUse)[2:ncol(fuelUse)]
  fuels <- c("Diesel"
             , "Petrol"
             , "LPG"
             , "Burning Oil (kerosene)"
             , "electricity"
             , "Gas Oil (red diesel)")
  
  # make long
  fuelUse <- fuelUse %>%
    pivot_longer(cols = !Farm_type)
  head(fuelUse)
  
  ## ------------ Notes --------------  ##
  ## 'Road fuel' consists of derv (diesel oil for road vehicles) and petrol
  ## ------------ ----- --------------  ##
  
  #### 10b - import fuel emissions table ####
  
  ## ------------ Notes --------------  ##
  ## the amount of CO2e emitted when using different fuels comes from 
  ## Greenhouse gas reporting: Conversion Factors 2019
  ## (https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/806025/Conversion-Factors-2019-Condensed-set-for-most-users.xls)
  ## the imported table list the different fuels, and kg CO2e emissions based on the use
  ## of one unit. For example, as seen in 'head(energyUseEFs)', one tonne of CNG being used
  ## would emit 2542.04 kg CO2e
  ## electricity use was derived from https://webarchive.nationalarchives.gov.uk/ukgwa/20130403061708mp_/http://archive.defra.gov.uk/environment/business/reporting/pdf/101006-guidelines-ghg-conversion-factors.xls (Hill et al., 2011)
  ## ------------ ----- --------------  ##
  
  energyUseEFs <- fread("data_in/energy/emis_factors.csv") %>% as.data.frame()
  head(energyUseEFs)
  energyUseEFs <- energyUseEFs %>%
    # interested in kg co2e
    dplyr::select(Fuel, Unit, "kg CO2e")
  head(energyUseEFs)
  
  # refine to only fuels listed in fuelUse.csv
  ## check which are match
  names(fuelUse)
  unique(energyUseEFs$Fuel)
  ## get fuel factors (non-elec)
  energyEFsShort <- energyUseEFs %>%
    filter(Fuel %in% c("LPG", "Petrol (average biofuel blend)", "Diesel (average biofuel blend)"
                       , "Gas oil", "Burning oil"))
  ### for diesel, petrol, burning oil (kerosene), heating oil, and gas oil (red diesel), get litres
  energyEFs <- bind_rows(energyEFsShort %>%
                           filter(Fuel %in% c("Petrol (average biofuel blend)", "Diesel (average biofuel blend)"
                                              , "Gas oil", "Burning oil"), 
                                  Unit == "litres")
                         ### and kg (convert from tonnes) for LPG
                         , energyEFsShort %>%
                           filter(Fuel == "LPG" 
                                  , Unit == "tonnes")) %>%
    # convert amounts to numeric
    mutate(across(3, ~as.numeric(.)))
  energyEFs$kgCO2e <- ifelse(energyEFs$Fuel == "LPG", energyEFs$`kg CO2e` / 1000, energyEFs$`kg CO2e`)
  # change the unit from tonnes to kg
  energyEFs$Unit[energyEFs$Unit == "tonnes"] <- "kg"
  energyEFs <- energyEFs %>% dplyr::select(-"kg CO2e")
  head(energyEFs)
  ## for road fuel, take an average of diesel and petrol
  energyEFs <- bind_rows(energyEFs, data.frame(Fuel = "road_fuel_l"
                                               , Unit = "litres"
                                               , kgCO2e = mean(c(2.5941, 2.2090))))
  
  ## get fuel factors (electricity) - 2008 value in Hill et al., 2011
  energyEFs <- bind_rows(energyEFs, data.frame(Fuel = "Electiricty"
                                               , Unit = "kWh"
                                               , kgCO2e = 0.54522))
  head(energyEFs)
  str(energyEFs)
  # tidy
  rm(energyEFsShort, energyUseEFs)
  
  ## rename to match fuel use
  energyEFs$fuelMatch <- c("kerosene_l", "diesel", "red_dies_l", "petrol", "lpg_kg", "road_fuel_l", "elec_units")
  head(energyEFs, 7)
  
  #### 10c - multiply volume used per hectare by emissions factors ####
  # merge both tables
  fuelMerge <- merge(energyEFs, fuelUse
                     , by.y = "name"
                     , by.x = "fuelMatch"
                     , all = F) %>%
    # multiply volume and factor
    mutate(kgCO2eHa = kgCO2e * value) %>%
    # sum each farm type to get a fuel use total CO2e per farm type
    group_by(Farm_type) %>%
    summarise(totalkgCO2eHa = sum(kgCO2eHa)) %>%
    as.data.frame()
  head(fuelMerge)
  
  #### 10d - determine hectarage of all animals ####
  
  ## ------------ Notes --------------  ##
  ## 'cereals' will account for all arable hectares
  ## livestock farming values will be divided based on the number of animals
  ## in a km2, and the recommended stocking rate from Nix (2021)
  ## ------------ ----- --------------  ##
  
  ##### original data and scenarios #####
  ## ------------ Notes --------------  ##
  ## do this for both original data and the scenarios
  ## ------------ ----- --------------  ##
  # list all data
  ## original = 1 + scens
  landArea.length <- 1 + length(scenarios)
  ## get scen names
  sNames <- gsub(".gpkg", "", basename(scenarios))
  ### and use them to load appropriate land cover
  landCoverHa.scens <- list.files(file.path("scenario", "scen_maps", "finer_detail")
                                  , pattern = ".*_land.gpkg"
                                  , full.names = T)
  landCoverHa.scens
  
  for(iscen in 1:landArea.length){
    
    if(iscen == 1){
      lcm <- st_read("data_in/land_cover/land_cover_table.gpkg")
      nm <- "original_data"
      saveName <- file.path(savePath, "fuel_emissions.gpkg")  
      animalLand <- st_read("data_in/animals/all_animals_1km.gpkg")
    } else {
      lcm <- st_read(landCoverHa.scens[[iscen-1]])
      nm <- gsub("_land.gpkg", "", basename(landCoverHa.scens[[iscen-1]]))
      saveName <- file.path("scenario", "results"
                            , paste0("fuel_emissions", nm, ".gpkg"))
      # get animal data from scenario
      animalLand <- list.files(file.path("scenario", "scen_maps", "finer_detail")
                               , pattern = paste0(nm, "_animal.gpkg")
                               , full.names = T)
      stopifnot(length(animalLand) == 1)
      animalLand <- st_read(animalLand)
    }
    cat("... Calculating emissions from on-farm energy use for", nm, "... \n")
    
    # read in animal numbers
    animals <- animalLand
    head(animals)
    names(animals)
    ## sum them, according to categories in Nix
    animals2 <- animals %>% st_drop_geometry() %>%
      ### dairy
      mutate(diary = rowSums(select(., contains("airy")))) %>%
      dplyr::select(-contains("airy")) %>% # drop summed
      ### In—calf heifers
      mutate(In_calf_heifs = rowSums(select(., contains("eifers")))) %>%
      dplyr::select(-contains("eifers")) %>%# drop summed
      ### In—calf heifers
      mutate(bull = rowSums(select(., contains("ulls")))) %>%
      dplyr::select(-contains("ulls")) %>% # drop summed
      ### beef
      mutate(beuf = rowSums(select(., contains("eef")))) %>%
      dplyr::select(-contains("eef")) %>% # drop summed
      ### Other cattle 0–1 years old
      mutate(cow0_1years = rowSums(select(., c(contains(c("calves"
                                                          , "6...9.months", "0...3.months", "9...12.months", "3...6.months")))))) %>%
      dplyr::select(-c(contains(c("calves"
                                  , "6...9.months", "0...3.months", "9...12.months", "3...6.months")))) %>% # drop summed
      ### Other cattle 1–2 years old
      mutate(cow1_2years = rowSums(select(., c(contains(c("12...15.months", "15...18.months", "18...21.months", "21...24.months")))))) %>%
      dplyr::select(-c(contains(c("12...15.months", "15...18.months", "18...21.months", "21...24.months")))) %>% # drop summed
      ### Other cattle >2 years old
      mutate(cowabove2years = rowSums(select(., c(contains("..."))))) %>%
      dplyr::select(-c(contains("..."))) %>%
      ### poultry
      mutate(poultry = rowSums(select(., c(contains("poultry"))))) %>%
      ### pigs
      mutate(all_pigs = rowSums(select(., c(contains("pigs"))))) %>%
      ### sheep
      mutate(sheep = rowSums(select(., contains(c("ewes", "rams", "lambs"))))) %>%
      dplyr::select(-contains(c("ewes", "rams", "lambs"))) %>% # drop summed
      ### remove non-totals for sheep, and remove horses
      dplyr::select(-c(contains(c("lambs", "rams", "ewe", "total.poultry..c48.", ".pigs.", "horse", "region")))) %>%
      rename(dairy = diary)
    head(animals2)
    names(animals2)
    
    # create Grazing Livestock Unit (GLU) table per animal group
    GLUanimals <- c("In_calf_heifs", "bull", "dairy", "beuf", "cow0_1years"
                    , "cow1_2years", "cowabove2years", "sheep", "poultry", "all_pigs")
    GLUvalue <- c(In_calf_heifs = 0.8
                  , bull = 0.65
                  , dairy = 1
                  , beuf = 0.75
                  , cow0_1years = 0.34
                  , cow1_2years = 0.65
                  , cowabove2years = 0.8
                  , sheep = 0.08
                  , poultry = 1/2500 # gov.uk welfare poultry recommendations https://www.gov.uk/government/publications/poultry-on-farm-welfare/poultry-welfare-recommendations
                  , all_pigs = 1/30 # according to https://redtractorassurance.org.uk/standards/outdoor-pigs/
    )
    animalGLU <- bind_cols(animal_type = GLUanimals, GLU = GLUvalue)
    animalGLU
    
    # determine hectares of each animal group
    ## rearrange columns to match animal list categories
    animals3 <- setcolorder(animals2, as.character(animalGLU$animal_type))
    ## check column headers and row names match
    stopifnot(names(animals3)[1:(ncol(animals3)-1)] == animalGLU$animal_type)
    animalGLU
    
    ## use GLUs to multiply - transpose before and after to match vector
    gluRecommHect <- t(t(animals3[, c(1:(ncol(animals3)-1))]) %>%
                         replace(is.na(.), 0) * animalGLU$GLU) %>%
      as.data.frame() %>%
      bind_cols(rcFid_1km = animals3$rcFid_1km, .)
    head(gluRecommHect)    
    head(animalGLU)
    head(animals3)
    
    ##### 10e - determine hectares of different farm types #####
    ## ------------ Notes --------------  ##
    ## the hectare area required for each animals category and their quantity
    ## needs to be separated into Robust Farm Types
    ## ------------ ----- --------------  ##
    
    gluRecommHect2 <- gluRecommHect %>%
      # all non-dairy cattle will be classed as 'Lowland Grazing Livestock' 
      mutate(lowlandGrazing = rowSums(select(., c(bull, beuf, cow0_1years, cow1_2years, cowabove2years)))
             , dairy = rowSums(select(., c(dairy, In_calf_heifs)))) %>%
      dplyr::select(-c(bull, beuf, cow0_1years, cow1_2years, cowabove2years, In_calf_heifs)) %>%
      # rename the rest
      rename(
        ## sheep become 'LFA Grazing Livestock'
        LFAGrazing = sheep
        , pigs = all_pigs
        , poultry = poultry)
    head(gluRecommHect2)
    
    ## adjust total hectares by comparison with lcm grassland area 
    gluRecommHect3 <- gluRecommHect2 %>%
      # all non-dairy cattle will be classed as 'Lowland Grazing Livestock' 
      mutate(totalHectares = rowSums(select(., c(LFAGrazing:lowlandGrazing))))
    head(gluRecommHect3)
    
    # compare against actual grassland from Land cover data
    names(lcm)
    lcmGrass <- lcm %>%
      mutate(totalGrassHa = rowSums(select(.[,,drop = T]
                                           , c(neutral_grass_ha, calc_grass_ha, acid_grass_ha
                                               , heather_grass_ha, improved_grass_ha)))) 
    
    # make all arable into 'cereals' - except potatoes
    lcmGrass <- lcmGrass %>%
      replace(is.na(.), 0) %>%
      dplyr::select(-potato_ha) %>%
      mutate(cereal_ha = rowSums(.[, which(colnames(lcmGrass) == "winterwheat_ha"):which(colnames(lcmGrass) == "sugarbeet_ha")
                                   , drop=TRUE], na.rm = TRUE))
    
    # merge animals and grass
    lcmGrassAnim <- merge(gluRecommHect3, lcmGrass %>% st_drop_geometry(), by = "rcFid_1km")
    head(lcmGrassAnim)
    
    # if the total area required by all the animals is higher that the grass land area identified in the lcm, 
    # adjust all hectarages by the multiplied difference
    hectTooHigh <- lcmGrassAnim %>%
      filter(totalHectares > totalGrassHa) %>%
      # get magnitude difference
      mutate(magDiff = totalGrassHa/totalHectares) %>%
      # just each animal hectares by that proportion
      mutate(across(LFAGrazing:lowlandGrazing, ~ . * magDiff)) %>%
      # recalculate total hectares - it should match grass hectares
      mutate(totalHectares = rowSums(select(., c(LFAGrazing:lowlandGrazing))))
    head(hectTooHigh)
    
    # combine back with rest
    finalHect <- bind_rows(hectTooHigh
                           , lcmGrassAnim %>%
                             filter(totalHectares <= totalGrassHa) %>% 
                             dplyr::select(c(rcFid_1km:totalHectares, cereal_ha, totalGrassHa))) 
    head(finalHect)
    
    ##### 10f - get emissions #####
    ## ------------ Notes --------------  ##
    ## the adjusted hectarage of each animal group can now be multipled by the emission
    ## factors that that type of farm typical uses (from 'fuelMerge')
    ## ------------ ----- --------------  ##
    
    # change fuelMerge names
    fuelMerge$farm <- ifelse(grepl("LFA", fuelMerge$Farm_type), "LFAGrazing"
                             , ifelse(grepl("Lowland", fuelMerge$Farm_type), "lowlandGrazing"
                                      , ifelse(grepl("ereal", fuelMerge$Farm_type), "cereal_ha"
                                               , tolower(fuelMerge$Farm_type))))
    # match column and row names
    finalHect2 <- finalHect %>% 
      relocate(fuelMerge$farm)
    head(finalHect2)
    # ensure fuelmerge and final hectarage in same order
    stopifnot(names(finalHect2)[1:nrow(fuelMerge)] == fuelMerge$farm)
    
    ## use fuelMerge to multiply - transpose before and after to match vector
    emissionsFuel <- t(t(finalHect2[, c(1:nrow(fuelMerge))]) %>%
                         replace(is.na(.), 0) * fuelMerge$totalkgCO2eHa) %>%
      bind_cols(., rcFid_1km = finalHect2$rcFid_1km) %>%
      as.data.frame() %>%
      # rename, making units known
      rename_with(~paste0(., "_kgCO2e"), 1:6) %>%
      # get total in kg CO2e
      mutate(totalkgCO2e = rowSums(select(., c(1:6)))) %>%
      # convert to tCO2e
      mutate(Fuel_tCO2e = totalkgCO2e/1000)
    head(emissionsFuel)
    
    # make spatial again
    emissionsFuelSpat <- merge(lcmGrass %>% dplyr::select(rcFid_1km)
                               , emissionsFuel
                               , by = "rcFid_1km")
    str(emissionsFuelSpat)
    
    ##### 10g - energy use for processing potatoes #####
    # Farmers Weekly suggests (based on AHDB data) that average crop yields
    # for potatoes are 41.7 t/ha, with a maximum yield of 80 t/ha. 
    # According to AHDB, C footprints are: 
    # 45.4, 39.2, 30.77 depending on the type of potato storage (per t)
    # The mean average value equals: 38.46 kgCO2e/tonne
    # assuming each ha produces 41.7 t, each ha growing potatoes = 1603.782 kg CO2e
    
    # get potato column
    crop_pot <- lcm %>%
      dplyr::select(rcFid_1km, potato_ha)
    # times by above emissions amount
    crop_pot <- crop_pot %>%
      mutate(potatoEmis_kgco2e = potato_ha * 1603.782)
    
    # add to final df
    emissionsFuelSpat <- emissionsFuelSpat %>%
      merge(., crop_pot %>% st_drop_geometry() %>% dplyr::select(-potato_ha)) %>%
      # get total in kg CO2e
      mutate(totalkgCO2e = totalkgCO2e + potatoEmis_kgco2e) %>%
      # convert to tCO2e
      mutate(Fuel_tCO2e = totalkgCO2e/1000)
    head(emissionsFuelSpat)
    
    # save totals file
    st_write(emissionsFuelSpat
             , saveName
             , append = F)
    fwrite(st_drop_geometry(emissionsFuelSpat)
           , gsub(".gpkg", ".csv", saveName)
           , row.names = F)
  }
  
  # write readme
  sink(file = NULL)
  sink(file = file.path(savePath, "readme_fuel_emissions.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2023-05-03
Last update:",  format(Sys.Date()), "

The files containing 'fuel_emissions' contain information the annual amount of emissions of CO2 per km2 due to on-farm fuel use, based on underlying land covers, and animal numbers and their recommended hectarage. The land cover areas were derived from the 2015 land cover map.
The amount required to process potatoes is also included.

The columns:
      'rcFid_1km' = 1 km reference id
      '[land cover type]_kgCO2e' = emissions of CO2 that each land cover contributes across a 1 km2 based on its fuel use (units: kg CO2/km2/yr)
      'totalkgCO2e' = total emissions of CO2 for a 1 km2, based on on-farm fuel use (units: kg CO2/km2/yr)
      'Fuel_tCO2e' = total emissions of CO2 for a 1 km2, based on on-farm fuel use (units: t CO2/km2/yr)
      
Data:
    units:  tonnes CO2/km2/yr ('Fuel_tCO2e')
            kg CO2/km2/yr ('[land cover type]_kgCO2e')
    Spatial resolution: 1 km2
    Spatial extent: GB
    Temporal coverage: 2015
    Native projection: 27700")
  sink(file = NULL)
  
  # write readme - for scenarios
  sink(file = NULL)
  sink(file = file.path("scenario", "results", "readme_fuel_emissions.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2023-05-03
Last update:",  format(Sys.Date()), "

The files containing 'fuel_emissions' contain information the annual amount of emissions of CO2 per km2 due to on-farm fuel use, based on underlying land covers, and animal numbers and their recommended hectarage. The land cover areas were derived from the 2007 scenarios land cover maps.
The amount required to process potatoes is also included.

The columns:
      'rcFid_1km' = 1 km reference id
      '[land cover type]_kgCO2e' = emissions of CO2 that each land cover contributes across a 1 km2 based on its fuel use (units: kg CO2/km2/yr)
      'totalkgCO2e' = total emissions of CO2 for a 1 km2, based on on-farm fuel use (units: kg CO2/km2/yr)
      'Fuel_tCO2e' = total emissions of CO2 for a 1 km2, based on on-farm fuel use (units: t CO2/km2/yr)
      
Data:
    units:  tonnes CO2/km2/yr ('Fuel_tCO2e')
            kg CO2/km2/yr ('[land cover type]_kgCO2e')
    Spatial resolution: 1 km2
    Spatial extent: GB
    Temporal coverage: 2007
    Native projection: 27700")
  sink(file = NULL)
}

#### 11 - part11combineemissions ####]
## ------------ Notes --------------  ##
## In this section, all emissions are combined.
## Where the reference input file is a '.csv', the calculations are carried 
## out below based on regional animal values
## Where the suffix is '.gpkg' spatial values were calculated above
## ------------ ----- --------------  ##

if(part11combineemissions){
  
  # stop("part11")
  
  # griddf = gridAnimals
  # rSel = entFerm
  # finalCol = "kgCO2e_entferm"
  # colofInterest = "kgCO2.50pcScen"
  
  # create a tranposed multiply function
  ## run it through a function, extracting one region at a time
  regionsTranspose <- function(griddf = "x"
                               , rSel = "x"
                               , finalCol = "y"
                               , colofInterest = "z"){
    
    ## get all different possible regions
    ukRegionsPossible <- unlist(unique(griddf %>% st_drop_geometry() %>%
                                         dplyr::select(region)) %>% as.vector()) %>% as.character()
    ### ignore NAs
    ukRegionsPossible <- ukRegionsPossible[!is.na(ukRegionsPossible)]
    
    # print(ukRegionsPossible)
    
    x <- ukRegionsPossible[[1]]
    # print(x)
    cat("      ...ensuring names and columns match\n")
    
    # split by regions current region
    regSplit <- pblapply(ukRegionsPossible, function(x){
      cat(x, "| ")
      currRegion <- x
      # use region to extract from over all grid
      regionGrid <- griddf %>% st_drop_geometry() %>%
        filter(region == currRegion)
      head(regionGrid)
      
      # use the same information for the other grid, but add 'All' as well
      regionSelect <- rSel %>% 
        filter(UK.Region %in% c(currRegion, "All"))
      head(regionSelect)
      
      # ensure lengths match
      stopifnot((ncol(regionGrid)-2) == nrow(regionSelect))
      
      # ensure colnames of the spatial grid and row names of the amounts match
      ## ensure that across and down are in the same positions,
      ## otherwise change their positions
      if(identical(
        c(names(regionGrid)[3:(ncol(regionGrid))])
        , c(regionSelect$animal)
      )
      ){
        regionalMix <- regionSelect
        cat("all names match\n")
        
      } else {
        
        # see if any names are not present
        wx <- colnames(regionGrid)[!which(colnames(regionGrid) %in% regionSelect$animal)]
        stopifnot(length(wx) == 0)
        
        regionalMix <- regionSelect[match(names(regionGrid)[3:(ncol(regionGrid))]
                                          , regionSelect$animal), ]   
        
        # recheck - identical?
        stopifnot(identical(
          c(names(regionGrid)[3:(ncol(regionGrid))])
          , c(regionalMix$animal)
        ))
        cat("names match after reshuffle\n")
      }
      head(regionalMix)
      
      # multiply the animals in that region by their emissions
      ## shorten regional animals emissions
      regionalMix0 <- regionalMix %>% as.data.frame()
      regionalMix0[is.na(regionalMix0)] <- 0
      head(regionalMix0)
      
      # shorten the animal numbers
      anr <- regionGrid %>%
        dplyr::select(-c(rcFid_1km, region)) %>%
        st_drop_geometry()
      names(anr)
      
      # remove horses / goats
      anr <- anr %>%
        dplyr::select(-c(contains(c("horses", "goats"))))
      names(anr)
      regionalMix0 <- regionalMix0 %>%
        filter(!grepl(paste(c("horses", "goats"), collapse = "|")
                      , animal))
      
      stopifnot(length(regionalMix0$animal[which(!regionalMix0$animal %in% names(anr))]) == 0)
      
      # get index of column of interest
      coltoAnalyse <- which(names(regionalMix0) == colofInterest)
      # multiply for result
      trya <- as.matrix(anr)
      dim(trya)
      dim(regionalMix0 %>%
            dplyr::select(all_of(coltoAnalyse)) %>%
            as.vector() %>%
            unlist() %>%
            as.data.frame())
      anrTimes <- as.matrix(anr) %*% diag(regionalMix0 %>%
                                            dplyr::select(all_of(coltoAnalyse)) %>%
                                            as.vector() %>%
                                            unlist()) %>%
        as.data.frame() %>%
        mutate(!!finalCol := rowSums(.[], na.rm = T)) 
      anrTimes[22005:22010, ]
      anrTimes[91:95, ]
      
      length(c(names(regionGrid)[3:(ncol(regionGrid))]))
      length(names(anrTimes)[1:(ncol(anrTimes)-1)])
      names(anrTimes)[1:(ncol(anrTimes)-1)] <- c(names(regionGrid)[3:(ncol(regionGrid))])
      head(anrTimes)
      
      # run some checks
      stopifnot(anr[20, 20] * regionalMix0[20, colofInterest] == anrTimes[20, 20])
      stopifnot(anr[100, 100] * regionalMix0[100, colofInterest] == anrTimes[100, 100])
      stopifnot(anr[200, 200] * regionalMix0[200, colofInterest] == anrTimes[200, 200])
      
      anrTimes <- anrTimes %>%
        # include id col back in
        bind_cols(regionGrid %>% dplyr::select(rcFid_1km), .)
      head(anrTimes)
      return(anrTimes)
      
    })
    
    # combine back together
    anrTimesbind <- bind_rows(regSplit)
    
    return(anrTimesbind)
  } # end of multiplying transposing function
  
  
  ##### 11.0 - select original data or scenario #####
  for(iscen in 1:(length(scenarios) + 1)){
    # for(iscen in 1:2){
    
    # stop("part 11 - combining emissions")
    
    # get the regions by themselves
    reg <- st_read(file.path("data_in", "animals", "nonBovine_grid_2020_1km.gpkg"), quiet = T) %>% 
      dplyr::select(rcFid_1km, region) %>% st_drop_geometry()
    ## convert non-england regions to 'other' - to match original cow regions
    reg$region <- ifelse(reg$region %in% 
                           c("Highlands and Islands", "Mid and West Wales"
                             , "West Scotland", "South Scotland", "North Wales"
                             , "Mid Scotland and Fife", "South Wales West"
                             , "Glasgow", "Central Scotland", "South Wales Central", "Lothian"
                             , "North East Scotland", "South Wales East")
                         , "other", reg$region)
    ## add underscores to regions, to make them match output df
    reg$region <- gsub(" ", "_", reg$region)
    ## and make London the same as South East
    reg$region <- ifelse(reg$region == "London", "South_East"
                         , reg$region)
    ## convert the humber, to match excrete N df
    reg$region <- ifelse(reg$region == "Yorkshire_and_The_Humber"
                         , "Yorkshire_and_the_Humber"
                         , reg$region)
    unique(reg$region)
    
    # if 'iscen' is 1, use 2015 data
    if(iscen == 1){
      cat("Starting 2015 combining emissions...\n")
      # set input location - for land cover
      currPath <- file.path("data_in", "animals")
      # set output location
      resultsPath <- file.path("results")
      # load grids - animals and land cover - with regions 
      gridAnimals.original <- st_read(file.path(currPath, "all_animals_1km.gpkg"), quiet = T) %>%
        st_drop_geometry()
      gridLand.original <- st_read(file.path("data_in", "land_cover", "land_cover_table.gpkg"), quiet = T) %>%
        st_drop_geometry() %>%
        merge(ghgGridGB, .)
      head(gridLand.original)
      
      gridAnimals <- gridAnimals.original
      gridLand <- gridLand.original
      
      rm(gridAnimals.original, gridLand.original)
      
    } else { # if > 1, use the scenarios
      scenName <- gsub(".gpkg", "", basename(scenarios[iscen-1]))
      cat("Starting", scenName,  "combining emissions...\n")
      # create saving directory
      dir.create(savePath <- file.path("scenario", "final_results"), showWarnings = F)
      resultsPath <- file.path("scenario", "results")
      # load grids - animals and land cover - with regions 
      landArea.scens <- st_read(file.path("scenario", "scen_maps", "finer_detail"
                                          , paste0(scenName, "_land.gpkg"))
                                , quiet = T) %>%
        st_drop_geometry()
      animalArea.scens <- st_read(file.path("scenario", "scen_maps", "finer_detail"
                                            , paste0(scenName, "_animal.gpkg"))
                                  , quiet = T) %>%
        st_drop_geometry()
      
      gridAnimals <- animalArea.scens
      gridLand <- landArea.scens
      
      rm(animalArea.scens, landArea.scens)
    }
    
    # add the appropriate regions
    gridAnimals <- gridAnimals %>%
      merge(., reg) %>% # merge, to get region
      relocate(rcFid_1km, region)
    # print(head(gridAnimals))
    
    gridLand <- gridLand %>%
      merge(., reg) %>% # merge, to get region
      relocate(rcFid_1km, region)
    # print(head(gridLand))
    
    #### 11b - enteric fermentation ####
    cat("   ...starting enteric fermentation...\n")
    # load in enteric fermentation - a per-animal emission value
    entFerm <- fread(file.path("results", "animals", "entericFermCH4CO2e.csv"))
    head(entFerm)
    
    # run through the transposing multiplication function
    totalEntFerm <- regionsTranspose(griddf = gridAnimals
                                     , rSel = entFerm
                                     , finalCol = "kgCO2e_entferm"
                                     , colofInterest = "kgCO2.50pcScen")
    head(totalEntFerm)
    class(totalEntFerm)
    
    ###### add to final dataframe ######
    cat("   ...adding enteric fermentation...\n")
    ghgGridResults <- totalEntFerm %>% dplyr::select(rcFid_1km, kgCO2e_entferm) %>%
      merge(ghgGridGB, .
            , by = "rcFid_1km")
    head(ghgGridResults)
    plot(ghgGridResults[1])
    
    #### 11c - manure management ####
    cat("   ...starting manure management...\n")
    # load in manure management - a per-animal emission value - it is a value for 6 months (indoors)
    manMgmt <- fread(file.path("results", "animals", "MMCO2e.csv"))
    head(manMgmt)
    
    # run through the transposing multiplication function
    totalmanMgmt <- regionsTranspose(griddf = gridAnimals
                                     , rSel = manMgmt
                                     , finalCol = "kgCO2e_manMgmt"
                                     , colofInterest = "CO2e_kg_yr")
    head(totalmanMgmt)
    
    ###### add to final dataframe ######
    cat("   ...adding manure management...\n")
    ghgGridResults <- totalmanMgmt %>% dplyr::select(rcFid_1km, kgCO2e_manMgmt) %>%
      merge(ghgGridResults, .
            , by = "rcFid_1km")
    head(ghgGridResults)
    
    #### 11d - animal excretions ####
    cat("   ...starting animal excretions...\n")
    # load in animal excretions - a per-animal emission value
    excrete <- fread(file.path("results", "animals", "annual_grazing_emissions.csv"))
    head(excrete)
    
    # run through the transposing multiplication function
    totalexcrete <- regionsTranspose(griddf = gridAnimals
                                     , rSel = excrete
                                     , finalCol = "kgCO2e_excrete"
                                     , colofInterest = "grazEmisYr_kgCO2e")
    head(totalexcrete)
    
    ###### add to final dataframe ######
    cat("   ...adding animal excretions...\n")
    ghgGridResults <- totalexcrete %>% dplyr::select(rcFid_1km, kgCO2e_excrete) %>%
      merge(ghgGridResults, .
            , by = "rcFid_1km")
    head(ghgGridResults)
    
    #### 11e - fertiliser use ####
    cat("   ...starting fertiliser use...\n")
    # load in fertiliser - an amount per km2
    # if 'iscen' is 1, use 2015 data
    if(iscen == 1){
      fertUse <- st_read(file.path(resultsPath, "arable", "fertiliser_emissions.gpkg"), quiet = T)
    } else {
      fertUse <- st_read(file.path(resultsPath
                                   , paste0("fertiliser_emissions", scenName, ".gpkg")), quiet = T)
    }
    head(fertUse)
    
    ###### add to final dataframe ######
    cat("   ...adding fertiliser use...\n")
    ghgGridResults <- fertUse %>% dplyr::select(fertliser_CO2e_kgkm2, rcFid_1km) %>%
      st_drop_geometry() %>%
      merge(ghgGridResults, .
            , by = "rcFid_1km"
            , all = T) %>%
      replace(is.na(.), 0)
    head(ghgGridResults)
    
    #### 11f - land use ####
    # stop("before land use")
    cat("   ...starting land use...\n")
    # load in land use emissions - an amount per km2 for non-ag areas
    if(iscen == 1){
      landUse <- st_read(file.path(resultsPath, "land_use", "land_use_emissions.gpkg")
                         , quiet = T) %>%
        # get a kg co2e value
        mutate(landCoverTotalEmis_kgco2 = landCoverTotalEmis_tco2 * 1000) %>%
        st_drop_geometry() %>% distinct()
    } else {
      landUse <- st_read(file.path(resultsPath
                                   , paste0("land_use_emissions", scenName, ".gpkg"))
                         , quiet = T) %>%
        # get a kg co2e value
        mutate(landCoverTotalEmis_kgco2 = landCoverTotalEmis_tco2 * 1000) %>%
        st_drop_geometry() %>% distinct()
    }
    head(landUse)
    
    ###### add to final dataframe ######
    cat("   ...adding land use...\n")
    ghgGridResults <- landUse %>% dplyr::select(landCoverTotalEmis_kgco2, rcFid_1km) %>%
      st_drop_geometry() %>%
      merge(ghgGridResults, .
            , by = "rcFid_1km"
            , all = T) %>%
      replace(is.na(.), 0)
    head(ghgGridResults)
    
    #### 11g - residue management ####
    cat("   ...starting residue management...\n")
    # load in residue management emissions - an amount per km2 
    if(iscen == 1){
      residues <- st_read(file.path(resultsPath, "arable", "residEmisComb_kgCO2e.gpkg")
                          , quiet = T) %>%
        # assume all burned
        mutate(totalResid_kgKm2 = Burned..100..)
    } else {
      residues <- st_read(file.path(resultsPath
                                    , paste0("residEmisComb_kgCO2e", scenName, ".gpkg"))
                          , quiet = T) %>%
        # assume all burned
        mutate(totalResid_kgKm2 = Burned..100..)
    }
    
    head(residues)
    
    ###### add to final dataframe ######
    cat("   ...adding residue management...\n")
    ghgGridResults <- residues %>% dplyr::select(totalResid_kgKm2, rcFid_1km) %>%
      st_drop_geometry() %>%
      merge(ghgGridResults, .
            , by = "rcFid_1km"
            , all = T) %>%
      replace(is.na(.), 0)
    head(ghgGridResults)
    
    #### 11h - on-farm energy ####
    cat("   ...starting on-farm energy...\n")
    # load in on-farm energy emissions - an amount per km2 
    if(iscen == 1){
      onFarmEnergy <- st_read(file.path(resultsPath, "land_use", "fuel_emissions.gpkg")
                              , quiet = T) %>%
        rename(fuel_kgco2e = totalkgCO2e)
    } else {
      onFarmEnergy <- st_read(file.path(resultsPath
                                        , paste0("fuel_emissions", scenName, ".gpkg"))
                              , quiet = T) %>%
        rename(fuel_kgco2e = totalkgCO2e)
    }
    
    head(onFarmEnergy)
    
    ###### add to final dataframe ######
    cat("   ...adding on-farm energy...\n")
    ghgGridResults <- onFarmEnergy %>% dplyr::select(fuel_kgco2e, rcFid_1km) %>%
      st_drop_geometry() %>%
      merge(ghgGridResults, .
            , by = "rcFid_1km"
            , all = T) %>%
      replace(is.na(.), 0)
    head(ghgGridResults)
    # str(ghgGridResults)
    
    #### 12i - sum all together ####
    cat("   ...summing all emissions together...\n")
    ghgGridResultsSum <- ghgGridResults %>% st_drop_geometry() %>%
      mutate(total_kgco2e = rowSums(dplyr::select(.[,,], c(kgCO2e_entferm:fuel_kgco2e)))) %>%
      # convert to tonnes
      mutate(total_Tco2e = total_kgco2e / 1000) %>%
      distinct()
    head(ghgGridResultsSum)
    # str(ghgGridResultsSum)
    max(ghgGridResultsSum$total_Tco2e, na.rm = T)
    
    attr(ghgGridResultsSum$total_kgco2e, "ATT") <- NULL
    attr(ghgGridResultsSum$total_Tco2e, "ATT") <- NULL
    # str(ghgGridResultsSum)
    
    #### 12j - inverse ####
    cat("   ...creating inverse emissions...\n")
    ghgGridResultsSum <- ghgGridResultsSum %>%
      mutate(total_Tco2e_inv = total_Tco2e * -1)
    head(ghgGridResultsSum)
    # str(ghgGridResultsSum)
    
    # make spatial... again
    ghgGridResultsSpat <- ghgGridResultsSum %>%
      distinct() %>%
      merge(ghgGridGB, .) 
    head(ghgGridResultsSpat)
    # str(ghgGridResultsSpat)
    
    #### 12k - save final output ####
    cat("   ...saving final emissions...\n")
    which(ghgGridResultsSpat$total_Tco2e == max(ghgGridResultsSpat$total_Tco2e, na.rm = T))
    ghgGridResultsSpat[197730, ]
    
    # make spatial (raster data)
    vpRast <- st_rasterize(ghgGridResultsSpat %>% dplyr::select(total_Tco2e_inv, geometry))
    plot(vpRast)
    empty_raster <- raster(extent(ghgGridResultsSpat), res = c(1000, 1000))
    vpRast <- fasterize::fasterize(ghgGridResultsSpat, empty_raster, field = "total_Tco2e_inv")
    crs(vpRast) <- "epsg:27700"
    
    rasterDf <- st_rasterize(ghgGridResultsSpat %>% dplyr::select(total_Tco2e, geometry)
                             , dx = 1000, dy = 1000
                             , crs = 27700)
    
    plot(rasterDf)
    
    # save outputs
    ## ------------ Notes --------------  ##
    ## Order of saves are:
    ## Geopackage of all layers of emissions
    ## all layers of emissions in csv form
    ## raster of final emissions in its inverse form (i.e., making emissions negative)
    ## raster of final emissions (i.e., making emissions positive)
    ## ------------ ----- --------------  ##
    
    if(iscen == 1){ # original data
      st_write(ghgGridResultsSpat, file.path(resultsPath, "ghg_final_emissions_CO2e.gpkg")
               , append = F, quiet = T)
      fwrite(st_drop_geometry(ghgGridResultsSpat)
             , file.path(resultsPath, "ghg_final_emissions_CO2e.csv")
             , row.names = F)
      cat("      ...saved", file.path(resultsPath, "ghg_final_emissions_CO2e.gpkg"), "\n")
      writeRaster(vpRast, file.path(resultsPath, "ghgfin_emissions_CO2e_inv.tif"), overwrite = T)
      cat("      ...saved", file.path(resultsPath, "ghgfin_emissions_CO2e_inv.tif"), "\n")   
      write_stars(rasterDf,  file.path(resultsPath, "ghg_final_emissions_CO2e.tif"))
      
    } else { # scenario data
      
      # if(scenName == "Ag_15"){stop("Ag_15")}
      
      st_write(ghgGridResultsSpat, file.path(savePath
                                             , paste0("ghg_final_emissions_CO2e", scenName, ".gpkg"))
               , append = F, quiet = T)
      fwrite(st_drop_geometry(ghgGridResultsSpat)
             , file.path(savePath
                         , paste0("ghg_final_emissions_CO2e", scenName, ".csv"))
             , row.names = F)
      cat("      ...saved", file.path(savePath
                                      , paste0("ghg_final_emissions_CO2e", scenName, ".gpkg")), "\n")
      writeRaster(vpRast, file.path(savePath
                                    , paste0("ghgfin_emissions_CO2e_inv", scenName, ".tif"))
                  , overwrite = T)
      cat("      ...saved", file.path(savePath
                                      , paste0("ghgfin_emissions_CO2e_inv", scenName, ".tif")), "\n")   
      
      write_stars(rasterDf,  file.path(savePath
                                       , paste0("ghg_final_emissions_CO2e", scenName, ".tif"))
      )
      
    }
    
    #### 12l - plot final output ####
    # set colour ramp
    rasterDf2 <- rasterDf %>%
      as.data.frame(xy = T)
    head(rasterDf2)
    
    breaks <- c(-1000, 0, 1000, 2000)  # Define your custom breaks
    colors <- c("#40B0A6", "white", "orange1", "#E66100")
    
    gg <- ggplot() +
      geom_stars(data = rasterDf, aes(x = x, y = y, fill = total_Tco2e)) +
      theme_void() +
      coord_equal() +
      scale_fill_gradientn(breaks = breaks, colours = colors
                           , na.value="white") +
      labs(fill = expression(paste("GHG emissions (t ", CO[2]*"e)"))) +
      theme(legend.text = element_text(size = 8))
    
    cat("   ...saving final emissions [as images]...\n")
    
    # stop("before images")
    
    if(iscen == 1){ # original data
      png("images/ghg_balance_wider.png"
          , res = 200
          , width = 1000, height = 1200)
      print(gg)
      dev.off()
      
      png("images/ghg_balance_narrower.png"
          , res = 200
          , width = 900, height = 1200)
      print(gg)
      dev.off()
    } else { # scenario data
      png(paste0("images/ghg_balance_wider", scenName, ".png")
          , res = 200
          , width = 1000, height = 1200)
      print(gg)
      dev.off()
      
      png(paste0("images/ghg_balance_narrower", scenName, ".png")
          , res = 200
          , width = 900, height = 1200)
      print(gg)
      dev.off()
    }
    
    #### 12m - replace entFerm emissions with UK-specific equations
    cat("   ...creating UK emissions...\n")
    # load in enteric fermentation - a per-animal emission value (uk-specific equation)
    entFermUK <- fread(file.path("results", "animals", "entericFerm_uk_CH4CO2e.csv"))
    head(entFermUK)
    
    # run through the transposing multiplication function
    totalEntFermUK <- regionsTranspose(griddf = gridAnimals
                                       , rSel = entFermUK
                                       , finalCol = "kgCO2e_uk_entferm"
                                       , colofInterest = "co2e_kg_year")
    head(totalEntFermUK)
    class(totalEntFermUK)
    
    head(ghgGridResultsSum)
    ghgGridResultsSum.uk <- ghgGridResultsSum %>%
      dplyr::select(-c(kgCO2e_entferm, total_kgco2e: total_Tco2e_inv)) %>%
      # add uk Enteric fermentation in
      merge(totalEntFermUK %>%
              dplyr::select(rcFid_1km, kgCO2e_uk_entferm))
    head(ghgGridResultsSum.uk)
    
    ## total calculation
    ghgGridResultsSum.uk <- ghgGridResultsSum.uk %>%
      mutate(total_kgco2e = rowSums(dplyr::select(.[,,], c(kgCO2e_manMgmt:kgCO2e_uk_entferm)))) %>%
      # convert to tonnes
      mutate(total_Tco2e = total_kgco2e / 1000) 
    head(ghgGridResultsSum.uk)
    # str(ghgGridResultsSum.uk)
    max(ghgGridResultsSum.uk$total_Tco2e, na.rm = T)
    ghgGridResultsSum.uk <- ghgGridResultsSum.uk %>%
      merge(ghgGridGB %>% dplyr::select(rcFid_1km), .)
    head(ghgGridResultsSum.uk)
    # str(ghgGridResultsSum.uk)
    
    ### create the inverse
    ghgGridResultsSum.uk <- ghgGridResultsSum.uk %>%
      mutate(total_Tco2e_inv = total_Tco2e * -1) %>%
      relocate(rcFid_1km:total_Tco2e, total_Tco2e_inv)
    head(ghgGridResultsSum.uk)
    # str(ghgGridResultsSum.uk)
    
    ### write
    if(iscen == 1){
      st_write(ghgGridResultsSum.uk, file.path("results", "ghg_emissions_CO2e_uk.gpkg"), append = F)
      fwrite(st_drop_geometry(ghgGridResultsSum.uk)
             , file.path("results", "ghg_emissions_CO2e_uk.csv")
             , row.names = F)
    } else { # scenario data
      
      st_write(ghgGridResultsSum.uk
               , file.path(savePath, paste0("ghg_emissions_CO2e_uk", scenName, ".gpkg"))
               , append = F)
      fwrite(st_drop_geometry(ghgGridResultsSum.uk)
             , file.path(savePath, paste0("ghg_emissions_CO2e_uk", scenName, ".csv"))
             , row.names = F)
    }
  }
  # stop("elephant")
}

if(runAllScenarios){
  maps <- list.files(file.path("scenario", "scen_maps")
                     , pattern = ".gpkg"
                     , full.names = T)
  for(i in 1:5){
    
    qq <- st_read(maps[[i]]) %>% 
      st_drop_geometry()
    names(qq)
    qq2 <- colSums(qq[4:ncol(qq)], na.rm = T)
    qname <- gsub(".*ps/(.+).gpkg", "\\1", maps[[i]])
    
    if(i == 1){
      tq <- rbind(scen = qname
                  , qq2 %>% as.data.frame())
    } else {
      tq <- bind_cols(tq, rbind(scen = qname
                                , qq2 %>% as.data.frame()))
    }
    
  }
  
  tqt <- t(tq) %>% as.data.frame() %>% dplyr::select(scen
                                                     , contains(c("Lowland_60...240.months_Beef.cows"
                                                                  , "Medium.dairy_Dairy.calves"
                                                                  , "fattening"
                                                                  , "winterwheat_ha"
                                                                  , "improved_grass_ha")))
  fwrite(tqt, "tqt.csv", row.names = F)
  
  #### 13 - part13combineemissions for scenarios ####
  # create saving directory
  dir.create(file.path("scenario", "final_results"), showWarnings = F)
  
  # stop("part 13 - scenario end")
  
  # create a tranposed multiply function
  ## run it through a function, extracting one region at a time
  regionsTranspose <- function(griddf = "x"
                               , rSel = "x"
                               , finalCol = "y"
                               , colofInterest = "z"){
    
    ## get all different possible regions
    ukRegionsPossible <- unlist(unique(griddf %>% st_drop_geometry() %>%
                                         dplyr::select(region)) %>% as.vector()) %>% as.character()
    ### ignore NAs
    ukRegionsPossible <- ukRegionsPossible[!is.na(ukRegionsPossible)]
    
    # print(ukRegionsPossible)
    
    x <- ukRegionsPossible[[1]]
    # print(x)
    
    # split by regions current region
    # regSplit <- pblapply(ukRegionsPossible, function(x){
    regSplit <- lapply(ukRegionsPossible, function(x){
      # cat(x, "| ")
      currRegion <- x
      # use region to extract from over all grid
      regionGrid <- griddf %>% st_drop_geometry() %>%
        filter(region == currRegion)
      head(regionGrid)
      
      # use the same information for the other grid, but add 'All' as well
      regionSelect <- rSel %>% 
        filter(UK.Region %in% c(currRegion, "All"))
      head(regionSelect)
      
      # ensure lengths match
      stopifnot((ncol(regionGrid)-2) == nrow(regionSelect))
      
      # ensure colnames of the spatial grid and row names of the amounts match
      ## ensure that across and down are in the same positions,
      ## otherwise change their positions
      if(identical(
        c(names(regionGrid)[3:(ncol(regionGrid))])
        , c(regionSelect$animal)
      )
      ){
        regionalMix <- regionSelect
        # cat("all names match\n")
        
      } else {
        
        # see if any names are not present
        wx <- colnames(regionGrid)[!which(colnames(regionGrid) %in% regionSelect$animal)]
        stopifnot(length(wx) == 0)
        
        regionalMix <- regionSelect[match(names(regionGrid)[3:(ncol(regionGrid))]
                                          , regionSelect$animal), ]   
        
        # recheck - identical?
        stopifnot(identical(
          c(names(regionGrid)[3:(ncol(regionGrid))])
          , c(regionalMix$animal)
        ))
        # cat("names match after reshuffle\n")
      }
      head(regionalMix)
      
      # multiply the animals in that region by their emissions
      ## shorten regional animals emissions
      regionalMix0 <- regionalMix %>% as.data.frame()
      regionalMix0[is.na(regionalMix0)] <- 0
      head(regionalMix0)
      
      # shorten the animal numbers
      anr <- regionGrid %>%
        dplyr::select(-c(rcFid_1km, region)) %>%
        st_drop_geometry()
      names(anr)
      
      # remove horses / goats
      anr <- anr %>%
        dplyr::select(-c(contains(c("horses", "goats"))))
      names(anr)
      regionalMix0 <- regionalMix0 %>%
        filter(!grepl(paste(c("horses", "goats"), collapse = "|")
                      , animal))
      
      stopifnot(length(regionalMix0$animal[which(!regionalMix0$animal %in% names(anr))]) == 0)
      
      # get index of column of interest
      coltoAnalyse <- which(names(regionalMix0) == colofInterest)
      # multiply for result
      trya <- as.matrix(anr)
      dim(trya)
      dim(regionalMix0 %>%
            dplyr::select(all_of(coltoAnalyse)) %>%
            as.vector() %>%
            unlist() %>%
            as.data.frame())
      anrTimes <- as.matrix(anr) %*% diag(regionalMix0 %>%
                                            dplyr::select(all_of(coltoAnalyse)) %>%
                                            as.vector() %>%
                                            unlist()) %>%
        as.data.frame() %>%
        mutate(!!finalCol := rowSums(.[], na.rm = T)) 
      anrTimes[22005:22010, ]
      anrTimes[91:95, ]
      
      length(c(names(regionGrid)[3:(ncol(regionGrid))]))
      length(names(anrTimes)[1:(ncol(anrTimes)-1)])
      names(anrTimes)[1:(ncol(anrTimes)-1)] <- c(names(regionGrid)[3:(ncol(regionGrid))])
      head(anrTimes)
      
      # run some checks
      stopifnot(anr[20, 20] * regionalMix0[20, colofInterest] == anrTimes[20, 20])
      stopifnot(anr[100, 100] * regionalMix0[100, colofInterest] == anrTimes[100, 100])
      stopifnot(anr[200, 200] * regionalMix0[200, colofInterest] == anrTimes[200, 200])
      
      anrTimes <- anrTimes %>%
        # include id col back in
        bind_cols(regionGrid %>% dplyr::select(rcFid_1km), .)
      head(anrTimes)
      return(anrTimes)
      
    })
    
    # combine back together
    anrTimesbind <- bind_rows(regSplit)
    
    return(anrTimesbind)
  } # end of multiplying transposing function
  
  ## ------------ Notes --------------  ##
  ## run final calculations on all the final scenarios
  ## ------------ ----- --------------  ##
  scenNames <- gsub(".gpkg", "", basename(scenarios))
  
  #### 13aS - load grids - animals and land cover - with regions ####
  # list the animal and land cover grids for each scenario
  landArea.scens <- list.files(file.path("scenario", "scen_maps", "finer_detail")
                               , pattern = "_land.gpkg"
                               , full.names = T)
  animalArea.scens <- list.files(file.path("scenario", "scen_maps", "finer_detail")
                                 , pattern = "_animal.gpkg"
                                 , full.names = T)
  
  ## ------------ Notes --------------  ##
  ## go through each scenario, one at a time
  ## ------------ ----- --------------  ##
  
  resultsPath <- file.path("scenario", "results")
  currPath <- file.path("data_in", "animals")
  
  # get general region
  reg <- st_read(file.path(currPath, "nonBovine_grid_2020_1km.gpkg")) %>%
    dplyr::select(rcFid_1km, region) %>% st_drop_geometry()
  head(reg)
  
  for(sn in scenNames){
    
    cat("\n\n##############################################"
        , "\nstarting", sn, "\n")
    
    # reduce grids to specific scenario
    landArea.curr <- landArea.scens[grepl(sn, landArea.scens)]
    animalArea.curr <- animalArea.scens[grepl(sn, animalArea.scens)]
    stopifnot(length(landArea.curr) == 1)
    stopifnot(length(animalArea.curr) == 1)
    
    cat("   animal grid ...\n")
    ## animal grid
    gridAnimals.scenario <- st_read(animalArea.curr, quiet = T) %>%
      st_drop_geometry() 
    gridAnimals.scenario <- gridAnimals.scenario %>%
      merge(., reg) %>% # merge, to get region
      relocate(rcFid_1km, region) %>%
      # remove NAs from regions
      filter(!is.na(region))
    head(gridAnimals.scenario)
    head(reg)
    
    # convert non-england regions to 'other' - to match original cow regions
    gridAnimals.scenario$region <- ifelse(gridAnimals.scenario$region %in% 
                                            c("Highlands and Islands", "Mid and West Wales"
                                              , "West Scotland", "South Scotland", "North Wales"
                                              , "Mid Scotland and Fife", "South Wales West"
                                              , "Glasgow", "Central Scotland", "South Wales Central", "Lothian"
                                              , "North East Scotland", "South Wales East")
                                          , "other", gridAnimals.scenario$region)
    
    # add underscores to regions, to make them match output df
    gridAnimals.scenario$region <- gsub(" ", "_", gridAnimals.scenario$region)
    # and make London the same as South East
    gridAnimals.scenario$region <- ifelse(gridAnimals.scenario$region == "London", "South_East"
                                          , gridAnimals.scenario$region)
    # convert the humber, to match excrete N df
    gridAnimals.scenario$region <- ifelse(gridAnimals.scenario$region == "Yorkshire_and_The_Humber"
                                          , "Yorkshire_and_the_Humber"
                                          , gridAnimals.scenario$region)
    unique(gridAnimals.scenario$region)
    head(gridAnimals.scenario)
    
    cat("      land grid ...\n")
    gridLand.scenario <- st_read(landArea.curr, quiet = T) %>%
      st_drop_geometry() %>%
      merge(ghgGridGB, .) %>%
      merge(gridAnimals.scenario %>% dplyr::select(region, rcFid_1km) %>% st_drop_geometry(), .)
    head(gridLand.scenario)
    
    #### 13bs - enteric fermentation ####
    cat("         calculating enteric fermentation ...\n")
    # load in enteric fermentation - a per-animal emission value
    entFerm <- fread(file.path("results", "animals", "entericFermCH4CO2e.csv")
                     , showProgress = FALSE)
    head(entFerm)
    
    # run through the transposing multiplication function
    totalEntFerm <- regionsTranspose(griddf = gridAnimals.scenario
                                     , rSel = entFerm
                                     , finalCol = "kgCO2e_entferm"
                                     , colofInterest = "kgCO2.50pcScen")
    head(totalEntFerm)
    class(totalEntFerm)  
    
    # add to final dataframe 
    ghgGridResults <- totalEntFerm %>% dplyr::select(rcFid_1km, kgCO2e_entferm) %>%
      merge(ghgGridGB, .
            , by = "rcFid_1km")
    head(ghgGridResults)
    plot(ghgGridResults[1])
    
    #### 13cs - manure management ####
    cat("            calculating manure emissions ...\n")
    # load in manure management - a per-animal emission value - it is a value for 6 months (indoors)
    manMgmt <- fread(file.path("results", "animals", "sixMonthMMCO2e.csv")
                     , showProgress = FALSE)
    head(manMgmt)
    
    # run through the transposing multiplication function
    totalmanMgmt <- regionsTranspose(griddf = gridAnimals.scenario
                                     , rSel = manMgmt
                                     , finalCol = "kgCO2e_manMgmt"
                                     , colofInterest = "MMkgCO2e6MnIn")
    head(totalmanMgmt)
    
    # add to final dataframe 
    ghgGridResults <- totalmanMgmt %>% dplyr::select(rcFid_1km, kgCO2e_manMgmt) %>%
      merge(ghgGridResults, .
            , by = "rcFid_1km")
    head(ghgGridResults)
    
    #### 13ds - animal excretions ####
    cat("               calculating animal excretions ...\n")
    # load in animal excretions - a per-animal emission value
    excrete <- fread(file.path("results", "animals", "annual_grazing_emissions.csv")
                     , showProgress = FALSE)
    head(excrete)
    
    # run through the transposing multiplication function
    totalexcrete <- regionsTranspose(griddf = gridAnimals.scenario
                                     , rSel = excrete
                                     , finalCol = "kgCO2e_excrete"
                                     , colofInterest = "grazEmisYr_kgCO2e")
    head(totalexcrete)
    
    # add to final dataframe 
    ghgGridResults <- totalexcrete %>% dplyr::select(rcFid_1km, kgCO2e_excrete) %>%
      merge(ghgGridResults, .
            , by = "rcFid_1km")
    head(ghgGridResults) 
    
    #### 13es - fertiliser use ####
    cat("                  calculating fertiliser use ...\n")
    # list fertiliser emissions, and reduce to current scenario
    fertEmis <- list.files(file.path("scenario", "results")
                           , pattern = paste0("fertiliser_emissions", sn, ".gpkg")
                           , full.names = T)
    stopifnot(length(fertEmis) == 1)
    
    # load in fertiliser - an amount per km2
    fertUse <- st_read(fertEmis, quiet = T)
    head(fertUse)
    
    # add to final dataframe 
    ghgGridResults <- fertUse %>% dplyr::select(fertliser_CO2e_kgkm2, rcFid_1km) %>%
      st_drop_geometry() %>%
      merge(ghgGridResults, .
            , by = "rcFid_1km"
            , all = T) %>%
      replace(is.na(.), 0)
    head(ghgGridResults)
    
    #### 13fs - land use ####
    cat("                     calculating land use emissions ...\n")
    
    # list land use emissions, and reduce to current scenario
    luEmis <- list.files(file.path("scenario", "results")
                         , pattern = paste0("land_use_emissions", sn, ".gpkg")
                         , full.names = T)
    stopifnot(length(luEmis) == 1)
    
    # load in land use emissions - an amount per km2 for non-ag areas
    landUse <- st_read(luEmis, quiet = T) %>%
      # get a kg co2e value
      mutate(landCoverTotalEmis_kgco2 = landCoverTotalEmis_tco2 * 1000)
    head(landUse)
    
    # add to final dataframe 
    ghgGridResults <- landUse %>% dplyr::select(landCoverTotalEmis_kgco2, rcFid_1km) %>%
      st_drop_geometry() %>%
      merge(ghgGridResults, .
            , by = "rcFid_1km"
            , all = T) %>%
      replace(is.na(.), 0)
    head(ghgGridResults)
    
    #### 13gs - residue management ####
    cat("                        calculating residue emissions ...\n")
    # list residue emissions, and reduce to current scenario
    resEmis <- list.files(file.path("scenario", "results")
                          , pattern = paste0("residEmisComb_kgCO2e", sn, ".gpkg")
                          , full.names = T)
    stopifnot(length(resEmis) == 1)
    
    # load in residue management emissions - an amount per km2 
    residues <- st_read(resEmis, quiet = T) %>%
      # assume all burned
      mutate(totalResid_kgKm2 = Burned..100..)
    head(residues)
    
    # add to final dataframe 
    ghgGridResults <- residues %>% dplyr::select(totalResid_kgKm2, rcFid_1km) %>%
      st_drop_geometry() %>%
      merge(ghgGridResults, .
            , by = "rcFid_1km"
            , all = T) %>%
      replace(is.na(.), 0)
    head(ghgGridResults)
    
    #### 13hs - on-farm energy ####
    cat("                           calculating energy emissions ...\n")
    # list on-farm emissions, and reduce to current scenario
    fuelEmis <- list.files(file.path("scenario", "results")
                           , pattern = paste0("fuel_emissions", sn, ".gpkg")
                           , full.names = T)
    stopifnot(length(fuelEmis) == 1)
    
    # load in on-farm energy emissions - an amount per km2 
    onFarmEnergy <- st_read(fuelEmis, quiet = T) %>%
      rename(fuel_kgco2e = totalkgCO2e)
    head(onFarmEnergy)
    
    # add to final dataframe 
    ghgGridResults <- onFarmEnergy %>% dplyr::select(fuel_kgco2e, rcFid_1km) %>%
      st_drop_geometry() %>%
      merge(ghgGridResults, .
            , by = "rcFid_1km"
            , all = T) %>%
      replace(is.na(.), 0)
    head(ghgGridResults)
    # str(ghgGridResults)  
    
    #### 13is - sum all together ####
    cat("                              summing all emissions ...\n")
    ghgGridResultsSum <- ghgGridResults %>% st_drop_geometry() %>%
      mutate(total_kgco2e = rowSums(dplyr::select(.[,,], c(kgCO2e_entferm:fuel_kgco2e)))) %>%
      # convert to tonnes
      mutate(total_Tco2e = total_kgco2e / 1000) 
    head(ghgGridResultsSum)
    # str(ghgGridResultsSum)
    # max(ghgGridResultsSum$total_Tco2e, na.rm = T)
    
    attr(ghgGridResultsSum$total_kgco2e, "ATT") <- NULL
    attr(ghgGridResultsSum$total_Tco2e, "ATT") <- NULL
    # str(ghgGridResultsSum)
    
    #### 13js - inverse ####
    cat("                                 summing inverse emissions ...\n")
    ghgGridResultsSum <- ghgGridResultsSum %>%
      mutate(total_Tco2e_inv = total_Tco2e * -1)
    head(ghgGridResultsSum)
    # str(ghgGridResultsSum)
    
    # make spatial... again
    ghgGridResultsSpat <- ghgGridResultsSum %>%
      merge(ghgGridGB, .)
    head(ghgGridResultsSpat)
    # str(ghgGridResultsSpat)
    
    #### 13ks - save final output ####
    cat("                                    saving emissions ...\n")
    st_write(ghgGridResultsSpat
             , file.path("scenario", "final_results"
                         , paste0(sn, "_final_emissions_CO2e.gpkg")), append = F)
    fwrite(st_drop_geometry(ghgGridResultsSpat)
           , file.path("scenario", "final_results"
                       , paste0(sn, "_final_emissions_CO2e.csv"))
           , row.names = F)
    
    which(ghgGridResultsSpat$total_Tco2e == max(ghgGridResultsSpat$total_Tco2e, na.rm = T))
    ghgGridResultsSpat[197730, ]
    
    cat("                                       creating maps ...\n")
    # make spatial (raster data)
    vpRast <- st_rasterize(ghgGridResultsSpat %>% dplyr::select(total_Tco2e_inv, geometry))
    plot(vpRast)
    empty_raster <- raster(extent(ghgGridResultsSpat), res = c(1000, 1000))
    vpRast <- fasterize::fasterize(ghgGridResultsSpat, empty_raster, field = "total_Tco2e_inv")
    crs(vpRast) <- "epsg:27700"
    # save raster
    writeRaster(vpRast
                , file.path("scenario", "final_results"
                            , paste0(sn, "_fin_emissions_CO2e_inv.tif"))
                , overwrite = T)
    
    rasterDf <- st_rasterize(ghgGridResultsSpat %>% dplyr::select(total_Tco2e, geometry)
                             , dx = 1000, dy = 1000
                             , crs = 27700)
    plot(rasterDf)
    write_stars(rasterDf
                , file.path("scenario", "final_results"
                            , paste0(sn, "_ghg_final_emissions_CO2e.tif"))
    )
    
    #### 13ls - plot final output ####
    # set colour ramp
    rasterDf2 <- rasterDf %>%
      as.data.frame(xy = T)
    head(rasterDf2)
    
    breaks <- c(-1000, 0, 1000, 2000)  # Define your custom breaks
    colors <- c("#40B0A6", "white", "orange1", "#E66100")
    
    gg <- ggplot() + 
      geom_stars(data = rasterDf, aes(x = x, y = y, fill = total_Tco2e)) +
      scale_fill_gradientn(
        colors = c("#40B0A6","white","#E66100"),
        values = scales::rescale(c(-500, 0 , 500)),
        limits = c(min(rasterDf2$total_Tco2e, na.rm = T)
                   ,  max(rasterDf2$total_Tco2e, na.rm = T))
        , na.value="white"
      ) +
      theme_void() +
      coord_equal() +
      labs(fill = expression(paste("GHG emissions (t ", CO[2]*"e)"))) +
      theme(legend.text = element_text(size = 8))
    
    png(paste0("images/ghg_balance_", sn, ".png")
        , res = 200
        , width = 1000, height = 1200)
    print(gg)
    dev.off()
    
    #### 13ms - convert to MtCO2e ####
    ghgResultsMtCO2e <- ghgGridResultsSpat %>%
      mutate(total_MtCO2e = total_Tco2e)
    
    #### 13ms - replace entFerm emissions with UK-specific equations
    # load in enteric fermentation - a per-animal emission value (uk-specific equation)
    entFermUK <- fread(file.path(resultsPath, "animals", "entericFerm_uk_CH4CO2e.csv"))
    head(entFermUK)
    
    # run through the transposing multiplication function
    totalEntFermUK <- regionsTranspose(griddf = gridAnimals.scenario
                                       , rSel = entFermUK
                                       , finalCol = "kgCO2e_uk_entferm"
                                       , colofInterest = "co2e_kg_year")
    head(totalEntFermUK)
    class(totalEntFermUK)
    
    head(ghgGridResultsSum)
    ghgGridResultsSum.uk <- ghgGridResultsSum %>%
      dplyr::select(-c(kgCO2e_entferm, total_kgco2e: total_Tco2e_inv)) %>%
      # add uk Enteric fermentation in
      merge(totalEntFermUK %>%
              dplyr::select(rcFid_1km, kgCO2e_uk_entferm))
    head(ghgGridResultsSum.uk)
    ## total calculation
    ghgGridResultsSum.uk <- ghgGridResultsSum.uk %>%
      mutate(total_kgco2e = rowSums(dplyr::select(.[,,], c(kgCO2e_manMgmt:kgCO2e_uk_entferm)))) %>%
      # convert to tonnes
      mutate(total_Tco2e = total_kgco2e / 1000)
    head(ghgGridResultsSum.uk)
    # str(ghgGridResultsSum.uk)
    # max(ghgGridResultsSum.uk$total_Tco2e, na.rm = T)
    ghgGridResultsSum.uk <- ghgGridResultsSum.uk %>%
      merge(ghgGridGB %>% dplyr::select(rcFid_1km), .)
    head(ghgGridResultsSum.uk)
    # str(ghgGridResultsSum.uk)
    
    ### create the inverse
    ghgGridResultsSum.uk <- ghgGridResultsSum.uk %>%
      mutate(total_Tco2e_inv = total_Tco2e * -1) %>%
      relocate(rcFid_1km:total_Tco2e, total_Tco2e_inv)
    head(ghgGridResultsSum.uk)
    str(ghgGridResultsSum.uk)
    
    ### write
    st_write(ghgGridResultsSum.uk
             , file.path("scenario", "final_results"
                         , paste0(sn, "_final_emissions_CO2e_uk.gpkg"))
             , append = F)
    fwrite(st_drop_geometry(ghgGridResultsSpat)
           , file.path("scenario", "final_results"
                       , paste0(sn, "_final_emissions_CO2e_uk.csv"))
           , row.names = F)
    
    which(ghgGridResultsSum.uk$total_Tco2e == max(ghgGridResultsSum.uk$total_Tco2e, na.rm = T))
    ghgGridResultsSum.uk[197730, ]
    
    # make spatial (raster data)
    vpRast <- st_rasterize(ghgGridResultsSum.uk %>% dplyr::select(total_Tco2e_inv, geometry))
    plot(vpRast)
    empty_raster <- raster(extent(ghgGridResultsSpat), res = c(1000, 1000))
    vpRast <- fasterize::fasterize(ghgGridResultsSpat, empty_raster, field = "total_Tco2e_inv")
    crs(vpRast) <- "epsg:27700"
    # save raster
    writeRaster(vpRast
                , file.path("scenario", "final_results"
                            , paste0(sn, "_fin_emisUK_CO2e_inv.tif"))
                , overwrite = T)
    
    rasterDf <- st_rasterize(ghgGridResultsSpat %>% dplyr::select(total_Tco2e, geometry)
                             , dx = 1000, dy = 1000
                             , crs = 27700)
    plot(rasterDf)
    write_stars(rasterDf
                , file.path("scenario", "final_results"
                            , paste0(sn, "_ghgfin_emisUK_CO2e.tif"))
    )
    
    ### plot results
    # set colour ramp
    rasterDf2 <- rasterDf %>%
      as.data.frame(xy = T)
    head(rasterDf2)
    
    max(rasterDf2$total_Tco2e, na.rm = T)
    
    breaks <- seq(-1000, 1000, 250)
    colors <- c("#40B0A6", "white", "orange1", "#E66100")
    
    gg <- ggplot() +
      geom_stars(data = rasterDf, aes(x = x, y = y, fill = total_Tco2e)) +
      scale_fill_gradientn(
        colors = c("#40B0A6","white","#E66100"),
        values = scales::rescale(c(-500, 0 , 500)),
        limits = c(min(rasterDf2$total_Tco2e, na.rm = T)
                   ,  max(rasterDf2$total_Tco2e, na.rm = T))
        , na.value="white"
      ) +
      labs(fill = expression(paste("GHG emissions (t ", CO[2]*"e)"))) +
      theme(legend.text = element_text(size = 8))
    
    png(paste0("images/ghg_balance_uk_", sn, ".png")
        , res = 200
        , width = 1000, height = 1200)
    print(gg)
    dev.off()
    
    cat("##############################################\n\n")
    # stop("summed")
  }
} # end of 'runAllScenarios'
