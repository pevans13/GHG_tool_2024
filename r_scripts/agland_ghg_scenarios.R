## ---------------------------
##
## Script name: agland_ghg_scenarios.R
##
## Purpose of script: run the GHG model on different scenarios of land cover.
##                    The land covers were derived from previously-published
##                    values (Redhead et al. 2020)
##                    
## ------------ Notes --------------  ##
## Reference for original scenarios: Redhead, J. W., Powney, G. D.,
## Woodcock, B. A., & Pywell, R. F. (2020). 
## Effects of future agricultural change scenarios on beneficial insects. 
## Journal of Environmental Management, 265, 110550.
## https://doi.org/10.1016/j.jenvman.2020.110550
## ------------ ----- --------------  ##
##
## Run after: total_ghg_script.R
##
## Specific numbered tasks:
## 1 - check different between land cover maps 2007 and 2015
## 2 - create maps for the new scenarios
## 3 - run the GHG model on the new scenarios
## 4 - produce different files maps and graphs based on scenario comparisons
##
## list of final outputs:
##    baseline2007.gpkg 
##    
##
## Author: Dr. Paul M. Evans
##
## Date Created: 2024-01-25
##
## Copyright (c) Paul M. Evans, 2024
## Email: paueva@ceh.ac.uk
## Git: https://github.com/pevans13

## ---------------------------
options(scipen = 6, digits = 4) # for non-scientific notation
## ---------------------------

# remove objects
rm(list = ls())

#### 0 - load libraries ####
## automatic install of packages if they are not installed already
list.of.packages <- c(
  "terra", "stars", "raster"
  , "dplyr", "sf", "tidyr"
  , "tictoc", "purrr"
  , "data.table"
  , "pbapply"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}
# loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(library(package.i, character.only = TRUE)
  )
}
# tidy
rm(list.of.packages, new.packages, package.i)

#### 0 - paths ####


#### 0 - run which parts? ####
runScenMapCreate <- T

#### 0 - load functions ####
# create a tranposed multiply function
## run it through a function, extracting one region at a time
regionsTranspose <- function(griddf = "x"
                             , rSel = "x"
                             , finalCol = "y"
                             , colofInterest = "z"){
  
  ## get all different possible regions
  ukRegionsPossible <- unlist(unique(griddf %>% st_drop_geometry() %>%
                                       dplyr::select(region)) %>% as.vector()) %>% as.character()
  # print(ukRegionsPossible)
  
  x <- ukRegionsPossible[[1]]
  
  # split by regions current region
  regSplit <- pblapply(ukRegionsPossible, function(x){
    cat(x, "| ")
    currRegion <- x
    # use region to extract from over all grid
    regionGrid <- griddf %>% st_drop_geometry() %>%
      filter(region == currRegion)
    head(regionGrid)
    
    # use the same information for the N excrete grid, but add 'All' as well
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
    
    head(regionGrid %>%
           dplyr::select(poultry_per_km), n = 100)
    
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

#### 1 - check different between land cover maps 2007 and 2015 ####
if(runScenMapCreate){
  if(file.exists(file.path("scenario", "scen_maps"
                           , "baseline2007.gpkg"))){
    df2007.final <- st_read(file.path("scenario", "scen_maps"
                                      , "baseline2007.gpkg"))
    df2007.raster <- rast(file.path("scenario", "scen_maps"
                                    , "baseline2007.tif"))
  } else {
    ## ------------ Notes --------------  ##
    ## load in both maps for 2007 and 2015, at the 1 km2 resolution. Due to the
    ## example map having only the 10 broad categories, the 2015 should have too
    ## ------------ ----- --------------  ##
    # 2007 map
    lcm2007 <- rast("ASSET v2/Landcover/baselineluNa_0_afNa_0.tif")
    # 2015 map (with 10 broad categories)
    lcm2015 <- rast("ASSET v2/Landcover/lcm2015_10cat.tif")
    ## create into dataframe
    df2007 <- terra::as.data.frame(lcm2007, xy = TRUE, na.rm = FALSE) %>%
      filter(!is.na(baselineluNa_0_afNa_0)) %>%
      mutate(x = round(x, 0)
             , y = round(y, 0))
    df2015 <- terra::as.data.frame(lcm2015, xy = TRUE, na.rm = FALSE) %>%
      filter(lcm2015_10cat > 0) %>%
      mutate(x = round(x, 0)
             , y = round(y, 0))
    head(df2007)
    head(df2015)
    ### merge both
    lcmDiff <- df2007 %>%
      rename(lcm2007 = baselineluNa_0_afNa_0) %>%
      merge(.
            , df2015 %>%
              rename(lcm2015 = lcm2015_10cat)
            , by = c("x", "y")
            , all = T)
    
    # determine which pixels have changed
    lcmDiff.change <- lcmDiff %>%
      mutate(different = ifelse(lcm2007 != lcm2015, 1, 0))
    table(lcmDiff.change$different)
    51133/190604
    ## convert to point data
    lcmDiff.point <- lcmDiff.change %>%
      st_as_sf(coords = c("x", "y"),
               crs = st_crs(lcm2015))
    lcmDiff.point <- lcmDiff.point %>%
      st_transform(., 27700)
    
    # determine which region each belongs to
    regions <- st_read("N:/Data/UK/Boundary/uk/eng_wales_scot.gpkg")
    
    ## extract region to point data
    lcmDiff.region <- lcmDiff.point %>%
      ## intersect polygons with points, keeping the information from both
      st_intersection(regions, .)
    head(lcmDiff.region)
    st_write(lcmDiff.region, "lcmDiff_region.gpkg", append = F)
    # lcmDiff.region <- st_read("lcmDiff_region.gpkg")
    
    ## ------------ Notes --------------  ##
    ## below, the finals of all animals, crops, etc. from the data inputs
    ## calculated in 'total_ghg_script.R' will be loaded in. Then, where the 
    ## dominant land type in 2007 does not match that in 2015, averages can 
    ## be used, and alterations made
    ## ------------ ----- --------------  ##
    
    ### load animals
    tghg.animals <- st_read("data_in/animals/all_animals_1km.gpkg")
    t(table(tghg.animals$rcFid_1km)) %>% as.data.frame() %>%
      filter(Freq > 1)
    ### load crops (and all land covers)
    tghg.crops <- st_read("data_in/land_cover/land_cover_table.gpkg")
    #### combine, using rcFid_1km
    tghg.ac <- tghg.animals %>%
      merge(., tghg.crops %>% st_drop_geometry()
            , by = "rcFid_1km") 
    
    tic("extraction to centroids")
    # get centroid, in order to be able to extract finer region
    tghg.centroid <- st_centroid(tghg.ac)
    names(tghg.centroid)
    ## extract region to point data
    tghg.region <- tghg.centroid %>%
      ## intersect polygons with points, keeping the information from both
      st_intersection(regions, .)
    toc()
    head(tghg.region)
    
    # save them all
    st_write(tghg.animals, "scenario/tghg_animals.gpkg", append = F)
    st_write(tghg.ac, "scenario/tghg_ac.gpkg", append = F)
    st_write(tghg.crops, "scenario/tghg_crops.gpkg", append = F)
    st_write(tghg.region, "scenario/tghg_region.gpkg", append = F)
    # tghg.region <- st_read("scenario/tghg_region.gpkg")
    # tghg.ac <- st_read("scenario/tghg_ac.gpkg")
    
    ### merge regions with animal and crop numbers
    tghg.ac2 <- tghg.ac %>% st_drop_geometry() %>%
      merge(., tghg.region %>% dplyr::select(Name, Area_Description, rcFid_1km) %>% st_drop_geometry()
            , by = "rcFid_1km"
            , all = F) %>%
      relocate(rcFid_1km, Name, Area_Description) %>%
      # remove duplicate rows
      distinct()
    st_write(tghg.ac2, "ac2.gpkg", append = T)
    # tghg.ac2 <- st_read("ac2.gpkg")
    ### make spatial
    tghg.ac3 <- tghg.ac2 %>%
      merge(tghg.ac %>%
              dplyr::select(rcFid_1km), .
            , by = "rcFid_1km")
    names(tghg.ac3)
    
    # merge with the ones that are the same between 2007 and 2015 and those that are not
    tic("extraction to centroids - region and change")
    ## extract region to point data
    tghg.ac4 <- lcmDiff.region %>% dplyr::select(-c(Name, Area_Description)) %>%
      ## intersect polygons with points, keeping the information from both
      st_intersection(tghg.ac3, .) %>%
      relocate(c(lcm2007, lcm2015, different), .after = rcFid_1km) %>%
      dplyr::select(-Area_Description)
    toc()
    names(tghg.ac4)
    
    # get averages of all animals and crops per region per land type
    tghg.summarise <- tghg.ac4 %>%
      # remove unneeded columns for summary
      dplyr::select(-c(rcFid_1km, lcm2007, different)) %>%
      st_drop_geometry() %>%
      group_by(lcm2015, Name) %>%
      mutate(n = n()) %>%
      relocate(lcm2015, Name, n) %>%
      summarise(across(everything(), ~mean(.)))
    head(tghg.summarise)
    fwrite(tghg.summarise, "scenario/tghg_summarise.csv", row.names = F)
    tghg.summarise <- fread("scenario/tghg_summarise.csv")
    names(tghg.summarise)
    
    # create the 2007 dataframe based on averages for the region / land cover type
    
    ## use county, and 2015 averages 
    df2007.changed <- tghg.ac4 %>%
      st_drop_geometry() %>%
      distinct() %>%
      # keep necessary columns
      dplyr::select(c(rcFid_1km, lcm2007, Name)) %>%
      # merge with summary, to get average animals and crops (from 2015 map)
      merge(., tghg.summarise %>% ungroup()
            , by.x = c("Name", "lcm2007")
            , by.y = c("Name", "lcm2015")
            , all = T) %>%
      filter(!is.na(rcFid_1km)) %>%
      mutate(across(everything(), ~tidyr::replace_na(., 0))) %>%
      distinct()
    fwrite(df2007.changed, "scenario/df2007_changed.csv", row.names = F)
    head(df2007.changed)
    class(df2007.changed)
    # make spatial 
    df2007.spat <- df2007.changed %>%
      # get x and y
      separate(rcFid_1km, c("x", "y"), sep = "_", remove = F) %>%
      # create as points
      st_as_sf(coords = c("x", "y")
               , crs = st_crs(27700))
    
    ### save 
    fwrite(df2007.changed
           , file.path("scenario", "scen_maps"
                       , "baseline2007.csv"), row.names = F)
    ### save as spatial
    st_write(df2007.spat
             , file.path("scenario", "scen_maps"
                         , "baseline2007.gpkg"), append = F)
    
    # write readme
    sink(file = NULL)
    sink(file = file.path("scenario", "scen_maps", "readme_baseline2007.md"))
    cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2024-01-30
Last update:",  format(Sys.Date()), "
Produced by 'agland_ghg_scenarios.R'

Description of 'baseline2007.csv':
  Contains the region, land cover classification, animal numbers, and land cover areas for each 1 km2. The values were derived by:
    the average 2015 values for that region (county), land cover classification combination was used.

Columns of 'baseline2007.csv':
  'rcFid_1km' = centroid point of 1 km reference id (EPSG: 27700) [x and y coordinates seperated by '_']
  'Name' = region of the the UK that the polygon centroid was in.  
  'lcm2007' = the dominant land cover classification that that 1 km2 polygon was assigned to in the land cover map 2007
  '[animal_name]' = how many of the animal type in each column name was found in that polygon. See main project for descriptions
  '[land_cover]_ha' = area, in hectares, of that land cover found in that 1 km2 pixel, as determined (in the original project) by using the 25 m2 land cover raster
  
Spatial resolution: 1 km2
Spatial extent: GB
Temporal coverage: 2007 (although data from animals and crops derived from 2020, and 2015, respectively)
Native projection: 27700")
    
    sink(file = NULL)
  }
  
  #### 2 - create maps for the new scenarios ####
  ## ------------ Notes --------------  ##
  ## The new maps produced here will be based on four of the scenarios from Redhead et al. (2020)
  ## The scenarios are a 15 and 30 per cent increase in agricultural expansion, and 15 and 30 per cent grassland restoration
  ## Here, first, the land covers for the different scenarios will be derived
  ## second, the land cover will be used to determine the animal and crop numbers per 1 km2
  ## ------------ ----- --------------  ##
  
  # read in scenarios maps
  scenMaps <- lapply(c(
    "ASSET v2/Landcover/scenarioluAg_15_afNa_0.tif"
    ,"ASSET v2/Landcover/scenarioluAg_30_afNa_0.tif"
    ,"ASSET v2/Landcover/scenarioluGr_15_afNa_0.tif"
    ,"ASSET v2/Landcover/scenarioluGr_30_afNa_0.tif"
  )
  , rast) %>%
    # stack them
    rast() %>%
    # create as a dataframe
    as.data.frame(., xy = T) %>%
    # create id
    mutate(rcFid_1km = paste(round(x), round(y), sep = "_")) %>%
    relocate(rcFid_1km) %>% dplyr::select(-c(x, y))
  names(scenMaps)
  ## keep only important colnames info
  names(scenMaps) <- gsub(".*riolu(.+)_af.*", "\\1", names(scenMaps))
  names(scenMaps)
  str(scenMaps)
  head(scenMaps)
  
  ## ------------ Notes --------------  ##
  ## for each scenario, merge with the county-averaged data from 2015 map
  
  ## requirements: df with rcFid_1km, region
  ## ------------ ----- --------------  ##
  
  # load in summarised info
  tghg.summarise <- fread("scenario/tghg_summarise.csv")
  head(tghg.summarise)
  sort(names(tghg.summarise))
  
  # load in 2007 basemap
  base2007 <- fread(file.path("scenario", "scen_maps"
                              , "baseline2007.csv")) %>%
    dplyr::select(Name, rcFid_1km)
  head(base2007)
  
  for(i in 2:ncol(scenMaps)){
    
    # isolate column for current scenario
    curr.scen <- scenMaps %>%
      dplyr::select(rcFid_1km, all_of(i)) %>%
      # add 'Name' (county name) via merge
      merge(., base2007)
    ## get scen
    cs <- names(curr.scen)[2]
    ## relpace with generic name
    names(curr.scen)[2] <- "lcm"
    names(curr.scen)
    
    cat("working on ...", cs, "...\n")
    
    ## merge with summary data from 2015
    curr.merge <- curr.scen %>%
      merge(., tghg.summarise
            , by.x = c("Name", "lcm")
            , by.y = c("Name", "lcm2015")) %>%
      # get x and y
      separate(rcFid_1km, c("x", "y"), sep = "_", remove = F) %>%
      # create as points
      st_as_sf(coords = c("x", "y")
               , crs = st_crs(27700))
    head(curr.merge)
    st_crs(curr.merge)
    
    # save
    ## as vector (point)
    cat("saving ...", cs, "...\n")
    st_write(curr.merge, file.path("scenario", "scen_maps"
                                   , paste0(cs, ".gpkg")), append = F)
    
    ## as raster
    curr.rast <- curr.merge %>%
      dplyr::select(lcm) %>%
      st_rasterize(dx = 1000, dy = 1000)
    curr.rast
    write_stars(curr.rast, file.path("scenario", "scen_maps"
                                     , paste0(cs, ".tif")))
  }
  
  # list all the output scenario rasters
  dfs2007 <- list.files(file.path("scenario", "scen_maps")
                        , pattern = ".tif$"
                        , full.names = T)
  cat(dfs2007, sep = "\n")
  stopifnot(length(dfs2007) == 4)
  
  # for each scenario, determine which pixels are different and the same
  for(i in dfs2007){
    print(i)
    
    i.rast <- rast(i)
    
    cat(table(i.rast %>% as.data.frame()), "\n")
  }
  
  # list all the output scenario vectors
  dfs2007 <- list.files(file.path("scenario", "scen_maps")
                        , pattern = ".gpkg$"
                        , full.names = T)
  cat(dfs2007, sep = "\n")
  stopifnot(length(dfs2007) == 5)
  
  # for each scenario, determine totals of a few variables (animals and land cover)
  dfs2007.amounts <- pblapply(dfs2007, function(i){
  
    print(i)
    
    # get name
    i.name <- gsub("^(.+)\\.gp.*", "\\1", basename(i))
    
    i.vect <- st_read(i, quiet = T)
    i.apply <- apply(i.vect %>% as.data.frame() %>%
                       dplyr::select(5:(last_col()-1)) %>%
                       mutate_all(~ ifelse(is.na(.), 0, .))
                       , 2, mean) %>%
      as.data.frame() 
    i.apply <- tibble::rownames_to_column(i.apply, "Name")
    names(i.apply)[2] <- i.name
    return(i.apply)
  }) %>% purrr::reduce(left_join, by = "Name")
  sort(dfs2007.amounts$Name)

  # select only a few rows for illustration
  dfs2007.select <- dfs2007.amounts %>%
    filter(Name %in% c("Medium.dairy_Dairy.cows", "poultry_per_km", "Continental_60...240.months_Beef.cows"
                       , "conifer_ha", "broadleaf_ha", "calc_grass_ha"
                       , "improved_grass_ha", "winterwheat_ha")) %>%
    relocate(Name, Ag_30)
  
  # keep those that are the same
  df2007.noChange <- df2007.scen %>%
    filter(change == 0) %>% dplyr::select(rcFid_1km)
  ## extract those values from the 2007 map
  df2007.noChange <- df2007.final %>%
    filter(rcFid_1km %in% df2007.noChange$rcFid_1km)
  
  # change the other to 2007 averages for region and land cover classification
  df2007.aver <- df2007.scen %>%
    filter(change == 1) %>%
    dplyr::select(rcFid_1km, Name, lc) %>%
    # merge with averages
    merge(., df2007.summarise %>% ungroup() %>% rename(lc = lcm2007)
          , all = T) %>%
    dplyr::select(-lc)
  
  # bind back together
  df2007Scen.bound <- df2007.noChange %>%
    bind_rows(., df2007.aver)
  # stopifnot(identical(nrow(df2007Scen.bound), nrow(df2007.changes)))
  cat(df2007.changes$rcFid_1km[which(!df2007.changes$rcFid_1km %in% df2007Scen.bound$rcFid_1km)]
      , "\n")
  
  ### save 
  st_write(df2007Scen.bound
           , file.path("scenario", "scen_maps"
                       , paste0(i, ".gpkg")), append = T)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path("scenario", "scen_maps", "readme_scenarios.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2024-01-30
Last update:",  format(Sys.Date()), "
Produced by 'agland_ghg_scenarios.R'

Description of scenario maps:
  Contains the region, land cover classification, animal numbers, and land cover areas for each 1 km2. The values were derived in one of two ways:
    1. If the land cover classification was the same as the 2007 baseline, the values were taken from the 2007 dataset (see full project)
    2. If the land cover classification was the different to the 2007 baseline, the average 2007 values for that region, land cover classification combination was used.

Naming of scenario maps:
  Naming follows '[scenario].gpkg', where:
    'scenario' defines the scenario '[ex]_[%]', which consists of what is being expanded (ex) and the per cent (%) amount.
      For 'ex':
        - 'Gr' relates to grassland restoration
        - 'Ag' relates to agricultural expansion

Columns of scenario maps:
  'rcFid_1km' = centroid point of 1 km reference id (EPSG: 27700) [x and y coordinates seperated by '_']
  'Name' = region of the the UK that the polygon centroid was in.  
  'lcm2007' = the dominant land cover classification that that 1 km2 polygon was assigned based on the scenario
  '[animal_name]' = how many of the animal type in each column name was found in that polygon. See main project for descriptions
  '[land_cover]_ha' = area, in hectares, of that land cover found in that 1 km2 pixel, as determined (in the original project) by using the 25 m2 land cover raster
  
Spatial resolution: 1 km2
Spatial extent: GB
Temporal coverage: 2007 (although data from animals and crops derived from 2020, and 2015, respectively)
Native projection: 27700")
  
  sink(file = NULL)
} # end 'runScenMapCreate'
rm(runScenMapCreate)

#### 3 - run the GHG model on the new scenarios (and baseline for 2007) ####
## ------------ Notes --------------  ##
## below follows the exact calculations from the main ghg script 'total_ghg_script.R',
## with only some alterations for reading in different sources (i.e. scenarios)
## It starts from part 4 of that script
## ------------ ----- --------------  ##

# read in initial grid
rcFid.grid <- st_read("data_in/ghgGridGB.gpkg")
rcFid.grid

# list all final files that will be assessed
mapsToEval.list <- list.files(file.path("scenario", "scen_maps")
                              , full.names = T
                              , pattern = ".gpkg")

# get rcFids to match regions (from original analysis)
## load in original regions
regions.broad <- fread("data_in/regions_assigned.csv") %>% as.data.frame() %>%
  distinct()

##### 3.0 - set-ups for analysis #####
###### 3.0a - set-up for animal-based emissions ######
# set-up for part 4 in original model - enteric fermentation
## read in MJ per day table
GEtableOut <- fread(file.path("results", "animals", "mjDayGE.csv"))
# load in table with Net Energy amounts per animal category
netEnergyTable <- fread(file.path("data_in", "animals", "animal_coefs_net_energy.csv")) %>%
  as.data.frame()
head(netEnergyTable)

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

# use the scenario of DEinside = 85 DE and DEoutside = 65 DE to give a per-animal enteric fermentation emission
CH4EntFerm <- DEtable2 %>%
  dplyr::select(animal, category, UK.Region, kgCH4GE.MJdayGEsml_65, kgCH4GE.MJdayGEhoused_85) %>%
  rowwise() %>%
  mutate(kgCH4yr.50pcScen = ifelse(grepl("lamb", category), kgCH4GE.MJdayGEhoused_85
                                   , (kgCH4GE.MJdayGEsml_65/2) + (kgCH4GE.MJdayGEhoused_85/2))) %>%
  # convert to CO2e
  mutate(kgCO2.50pcScen = kgCH4yr.50pcScen * 25)
head(CH4EntFerm)

###### 3.0b - set-up for fertiliser emissions ######
# set-up for part 7 in original model - fertiliser use
# load in soils data
soils <- st_read(file.path("data_in", "crops", "soils.gpkg"))
head(soils)

## ------------ Notes --------------  ##
## From RB209: 'For each field, the amount of lime to apply will depend on the current soil pH,
## soil texture, soil organic matter and the target pH, which should be 0.2 pH points above optimum'
## Therefore, Table 1.2 in RB209 was used to determine how much lime would be applied
## ------------ ----- --------------  ##

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

# create a rainfall grid from 2012 - 2017 precip data. Use the 'all_mean' data
rainExtract <- raster("N:/Data/UK/chessmet/precipitation/year_sum_averageikm.tif") %>%
  # convert to df
  as.data.frame(xy=T) %>%
  # using the info on Step 3 from RB209 (pg 9), put rain into three categories
  mutate(rainCat = if_else(layer < 600, "low"
                           , if_else(layer <= 700, "moderate"
                                     , "high")))

# separate centroid points of the soil data
separatedSoils <- soils %>%
  mutate(X = as.integer(unlist(map(soils$geom, 1))),
         Y = as.integer(unlist(map(soils$geom, 2)))) %>%
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

# read in and then replicate all cereals for the cereal types
snsParameters <- fread(file.path("data_in", "crops", "snsTableCompleted.csv"))
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
# read in table
snsCropParameters <- fread(file.path("data_in", "crops", "snsCropTableCompleted.csv"))
head(snsCropParameters)

# make tables wide
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


# for(i in 1:length(mapsToEval.list)){
for(i in 1){
  if(i == 1){
    i.list <- list()
  }
  
  cat("loading, and working with", mapsToEval.list[[i]], "\n")
  
  # read in polygons with animal and crop numbers
  mapsToEval.in <- st_read(mapsToEval.list[[i]]) %>%
    dplyr::select(-c(lcm2007, Name))
  # get regions
  length(which(mapsToEval.in$rcFid_1km %in% regions.broad$rcFid_1km))
  length(which(!mapsToEval.in$rcFid_1km %in% regions.broad$rcFid_1km))
  mapsToEval.in2 <- merge(mapsToEval.in, regions.broad, by = "rcFid_1km") %>%
    mutate(region_broad = ifelse(grepl("Scotland|Wales|Highlands|Glasgow|Lothian", region), "other", region)) %>%
    dplyr::select(-region)
  
  # add underscores to regions, to make them match output df
  mapsToEval.in2$region_broad <- gsub(" ", "_", mapsToEval.in2$region_broad)
  # and make London the same as South East
  mapsToEval.in2$region_broad <- ifelse(mapsToEval.in2$region_broad == "London", "South_East"
                                        , mapsToEval.in2$region_broad)
  # convert the humber, to match 
  mapsToEval.in2$region_broad <- ifelse(mapsToEval.in2$region_broad == "Yorkshire_and_The_Humber"
                                        , "Yorkshire_and_the_Humber"
                                        , mapsToEval.in2$region_broad)
  
  ## get all differ## get all different possible regions
  ukRegionsPossible <- unique(mapsToEval.in2$region_broad)
  
  ##### 3a - enteric fermentation #####
  # split by regions current region
  x <- ukRegionsPossible[[2]]
  regSplit <- pblapply(ukRegionsPossible, function(x){
    cat(x, "| ")
    currRegion <- x
    # use region to extract from over all grid
    regionGrid <- mapsToEval.in2 %>% st_drop_geometry() %>%
      # remove land covers
      dplyr::select(-contains("_ha")) %>%
      filter(region_broad == currRegion)
    head(regionGrid)
    
    # use the same information for the animal attributes grid, but add 'All' as well
    regionSelect <- CH4EntFerm %>% 
      filter(UK.Region %in% c(currRegion, "All"))
    head(regionSelect)
    
    # ensure lengths match
    stopifnot(identical(
      ncol(regionGrid %>% dplyr::select(-c(region_broad, rcFid_1km)))
      , nrow(regionSelect)
    ))
    
    # ensure colnames of the spatial grid and row names of the amounts match
    ## ensure that across and down are in the same positions,
    ## otherwise change their positions
    if(identical(
      c(names(regionGrid %>% dplyr::select(-c(region_broad, rcFid_1km))))
      , c(regionSelect$animal)
    ))
    {
      regionalMix <- regionSelect
      cat("all names match\n")
      
    } else {
      
      # see if any names are not present
      wx <- colnames(regionGrid)[!which(colnames(regionGrid) %in% regionSelect$animal)]
      stopifnot(length(wx) == 0)
      
      regionalMix <- regionSelect[match(names(regionGrid %>% dplyr::select(-c(region_broad, rcFid_1km)))
                                        , regionSelect$animal), ]   
      
      # recheck - identical?
      stopifnot(identical(
        c(names(regionGrid %>% dplyr::select(-c(region_broad, rcFid_1km))))
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
      dplyr::select(-c(rcFid_1km, region_broad)) %>%
      st_drop_geometry()
    names(anr)
    
    stopifnot(length(regionalMix0$animal[which(!regionalMix0$animal %in% names(anr))]) == 0)
    
    # get index of column of interest
    coltoAnalyse <- which(names(regionalMix0) == "kgCO2.50pcScen")
    # multiply for result
    trya <- as.matrix(anr)
    dim(trya)
    dim(regionalMix0 %>%
          dplyr::select(all_of(coltoAnalyse)) %>%
          as.vector() %>%
          unlist() %>%
          as.data.frame())
    finalCol <- "total" # set row total name
    anrTimes <- as.matrix(anr) %*% diag(regionalMix0 %>%
                                          dplyr::select(all_of(coltoAnalyse)) %>%
                                          as.vector() %>%
                                          unlist()) %>%
      as.data.frame() %>%
      mutate(!!finalCol := rowSums(.[], na.rm = T)) 
    anrTimes[22005:22010, ]
    anrTimes[91:95, ]
    
    # assign names back
    names(anrTimes)[1:(ncol(anrTimes)-1)] <- c(names(regionGrid  %>%
                                                       dplyr::select(-c(rcFid_1km, region_broad))))
    head(anrTimes)
    
    # run some checks
    stopifnot(anr[20, 20] * regionalMix0[20, "kgCO2.50pcScen"] == anrTimes[20, 20])
    stopifnot(anr[100, 100] * regionalMix0[100, "kgCO2.50pcScen"] == anrTimes[100, 100])
    stopifnot(anr[200, 200] * regionalMix0[200, "kgCO2.50pcScen"] == anrTimes[200, 200])
    
    anrTimes <- anrTimes %>%
      # include id col back in
      bind_cols(regionGrid %>% dplyr::select(rcFid_1km), .)
    head(anrTimes)
    return(anrTimes)
  })
  
  # combine back together
  anrTimesbind <- bind_rows(regSplit) %>%
    merge(rcFid.grid, .
          , by = "rcFid_1km")
  str(anrTimesbind)
  class(anrTimesbind)
  
  i.list[[i]] <- anrTimesbind
  
  st_write(anrTimesbind, file.path("scenario", "results"
                                   , paste0("EntFer_", basename(mapsToEval.list[[i]])))
           , append=TRUE)
  
  ##### 3b - manure management #####
  # load in manure produced per animal per region
  sixMonthMMCO2e <- fread(file.path("results", "animals", "sixMonthMMCO2e.csv")) %>%
    as.data.frame()
  analyseColumn <- "MMkgCO2e6MnIn"
  
  # split by regions current region
  x <- ukRegionsPossible[[2]]
  regSplit <- pblapply(ukRegionsPossible, function(x){
    cat(x, "| ")
    currRegion <- x
    # use region to extract from over all grid
    regionGrid <- mapsToEval.in2 %>% st_drop_geometry() %>%
      # remove land covers
      dplyr::select(-contains("_ha")) %>%
      filter(region_broad == currRegion)
    head(regionGrid)
    table(colSums(regionGrid %>% dplyr::select(where(is.numeric)))) %>% as.data.frame()
    
    # use the same information for the animal attributes grid, but add 'All' as well
    regionSelect <- sixMonthMMCO2e %>% 
      filter(UK.Region %in% c(currRegion, "All"))
    head(regionSelect)
    
    # ensure lengths match
    stopifnot(identical(
      ncol(regionGrid %>% dplyr::select(-c(region_broad, rcFid_1km)))
      , nrow(regionSelect)
    ))
    
    # ensure colnames of the spatial grid and row names of the amounts match
    ## ensure that across and down are in the same positions,
    ## otherwise change their positions
    if(identical(
      c(names(regionGrid %>% dplyr::select(-c(region_broad, rcFid_1km))))
      , c(regionSelect$animal)
    ))
    {
      regionalMix <- regionSelect
      cat("all names match\n")
      
    } else {
      
      # see if any names are not present
      wx <- colnames(regionGrid)[!which(colnames(regionGrid) %in% regionSelect$animal)]
      stopifnot(length(wx) == 0)
      
      regionalMix <- regionSelect[match(names(regionGrid %>% dplyr::select(-c(region_broad, rcFid_1km)))
                                        , regionSelect$animal), ]   
      
      # recheck - identical?
      stopifnot(identical(
        c(names(regionGrid %>% dplyr::select(-c(region_broad, rcFid_1km))))
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
      dplyr::select(-c(rcFid_1km, region_broad)) %>%
      st_drop_geometry()
    names(anr)
    
    stopifnot(length(regionalMix0$animal[which(!regionalMix0$animal %in% names(anr))]) == 0)
    
    # get index of column of interest
    coltoAnalyse <- which(names(regionalMix0) == analyseColumn)
    # multiply for result
    trya <- as.matrix(anr)
    dim(trya)
    dim(regionalMix0 %>%
          dplyr::select(all_of(coltoAnalyse)) %>%
          as.vector() %>%
          unlist() %>%
          as.data.frame())
    finalCol <- "total" # set row total name
    anrTimes <- as.matrix(anr) %*% diag(regionalMix0 %>%
                                          dplyr::select(all_of(coltoAnalyse)) %>%
                                          as.vector() %>%
                                          unlist()) %>%
      as.data.frame() %>%
      mutate(!!finalCol := rowSums(.[], na.rm = T)) 
    anr[22005:22010, ]
    anrTimes[22005:22010, ]
    anrTimes[91:95, ]
    
    # assign names back
    names(anrTimes)[1:(ncol(anrTimes)-1)] <- c(names(regionGrid  %>%
                                                       dplyr::select(-c(rcFid_1km, region_broad))))
    head(anrTimes)
    
    # run some checks
    stopifnot(anr[20, 20] * regionalMix0[20, analyseColumn] == anrTimes[20, 20])
    stopifnot(anr[100, 100] * regionalMix0[100, analyseColumn] == anrTimes[100, 100])
    stopifnot(anr[200, 200] * regionalMix0[200, analyseColumn] == anrTimes[200, 200])
    
    anrTimes <- anrTimes %>%
      # include id col back in
      bind_cols(regionGrid %>% dplyr::select(rcFid_1km), .)
    head(anrTimes)
    return(anrTimes)
  })
  
  # combine back together
  anrTimesbind <- bind_rows(regSplit)
  
  st_write(anrTimesbind, file.path("scenario", "results"
                                   , paste0("ManMgmt_", basename(mapsToEval.list[[i]])))
           , append=TRUE)
  
  ##### 3c - animal excretions #####
  # load in animal excretions - a per-animal emission value
  excrete <- fread(file.path("results", "animals", "annual_grazing_emissions.csv"))
  head(excrete)
  analyseColumn <- "grazEmisYr_kgCO2e"
  dfIn.section <- excrete
  
  # split by regions current region
  x <- ukRegionsPossible[[2]]
  regSplit <- pblapply(ukRegionsPossible, function(x){
    cat(x, "| ")
    currRegion <- x
    # use region to extract from over all grid
    regionGrid <- mapsToEval.in2 %>% st_drop_geometry() %>%
      # remove land covers
      dplyr::select(-contains("_ha")) %>%
      filter(region_broad == currRegion)
    head(regionGrid)
    table(colSums(regionGrid %>% dplyr::select(where(is.numeric)))) %>% as.data.frame()
    
    # use the same information for the animal attributes grid, but add 'All' as well
    regionSelect <- dfIn.section %>% 
      filter(UK.Region %in% c(currRegion, "All"))
    head(regionSelect)
    
    # ensure lengths match
    stopifnot(identical(
      ncol(regionGrid %>% dplyr::select(-c(region_broad, rcFid_1km)))
      , nrow(regionSelect)
    ))
    
    # ensure colnames of the spatial grid and row names of the amounts match
    ## ensure that across and down are in the same positions,
    ## otherwise change their positions
    if(identical(
      c(names(regionGrid %>% dplyr::select(-c(region_broad, rcFid_1km))))
      , c(regionSelect$animal)
    ))
    {
      regionalMix <- regionSelect
      cat("all names match\n")
      
    } else {
      
      # see if any names are not present
      wx <- colnames(regionGrid)[!which(colnames(regionGrid) %in% regionSelect$animal)]
      stopifnot(length(wx) == 0)
      
      regionalMix <- regionSelect[match(names(regionGrid %>% dplyr::select(-c(region_broad, rcFid_1km)))
                                        , regionSelect$animal), ]   
      
      # recheck - identical?
      stopifnot(identical(
        c(names(regionGrid %>% dplyr::select(-c(region_broad, rcFid_1km))))
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
      dplyr::select(-c(rcFid_1km, region_broad)) %>%
      st_drop_geometry()
    names(anr)
    
    stopifnot(length(regionalMix0$animal[which(!regionalMix0$animal %in% names(anr))]) == 0)
    
    # get index of column of interest
    coltoAnalyse <- which(names(regionalMix0) == analyseColumn)
    # multiply for result
    trya <- as.matrix(anr)
    dim(trya)
    dim(regionalMix0 %>%
          dplyr::select(all_of(coltoAnalyse)) %>%
          as.vector() %>%
          unlist() %>%
          as.data.frame())
    finalCol <- "total" # set row total name
    anrTimes <- as.matrix(anr) %*% diag(regionalMix0 %>%
                                          dplyr::select(all_of(coltoAnalyse)) %>%
                                          as.vector() %>%
                                          unlist()) %>%
      as.data.frame() %>%
      mutate(!!finalCol := rowSums(.[], na.rm = T)) 
    anr[22005:22010, ]
    anrTimes[22005:22010, ]
    anrTimes[91:95, ]
    
    # assign names back
    names(anrTimes)[1:(ncol(anrTimes)-1)] <- c(names(regionGrid  %>%
                                                       dplyr::select(-c(rcFid_1km, region_broad))))
    head(anrTimes)
    
    # run some checks
    stopifnot(anr[20, 20] * regionalMix0[20, analyseColumn] == anrTimes[20, 20])
    stopifnot(anr[100, 100] * regionalMix0[100, analyseColumn] == anrTimes[100, 100])
    stopifnot(anr[2000, 2000] * regionalMix0[2000, analyseColumn] == anrTimes[2000, 2000])
    
    anrTimes <- anrTimes %>%
      # include id col back in
      bind_cols(regionGrid %>% dplyr::select(rcFid_1km), .)
    head(anrTimes)
    return(anrTimes)
  })
  
  # combine back together
  anrTimesbind <- bind_rows(regSplit)
  
  st_write(anrTimesbind, file.path("scenario", "results"
                                   , paste0("excrete_", basename(mapsToEval.list[[i]])))
           , append=TRUE)
} # end of 'i'

##### 3d - fertiliser use #####
mapsToEval.crops <- mapsToEval.in2 %>%
  # keep land covers
  dplyr::select(rcFid_1km, region_broad, contains("_ha")) 
names(mapsToEval.crops)

# limeRegs are in t/ha, so determine amount for km2 depending on arable and grass areas
limeReqsAmounts <- mapsToEval.crops %>% dplyr::select(rcFid_1km, winterwheat_ha:sugarbeet_ha, improved_grass_ha) %>%
  st_drop_geometry() %>%
  # get total arable area
  mutate(total_arab_ha = rowSums(dplyr::select(., c(winterwheat_ha:sugarbeet_ha)), na.rm = T)) %>%
  mutate(total_arab_ha = ifelse(is.na(total_arab_ha), 0, total_arab_ha)) %>%
  rename(total_grass_ha = improved_grass_ha) %>%
  dplyr::select(rcFid_1km, total_arab_ha, total_grass_ha) %>%
  # make spatial again
  merge(rcFid.grid %>% dplyr::select(rcFid_1km), .) %>%
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
# https://naei.beis.gov.uk/data/ef-all-results?q=184177 states that the emission factor for limestone is 0.12 t CO2e t−1
limeReqsAmounts$limeTco2e <- limeReqsAmounts$limeAmount_tonnes * 0.12
head(limeReqsAmounts)

# adjust crop data
cropSoils <- mapsToEval.crops %>%
  st_centroid()
cropSoils <- cropSoils %>%
  mutate(X = as.integer(unlist(map(cropSoils$geom, 1))),
         Y = as.integer(unlist(map(cropSoils$geom, 2))))
head(cropSoils)

# merge with crop data
cropSoilRain <- rootDepthUK %>%
  merge(., cropSoils %>% st_drop_geometry() %>% dplyr::select(rcFid_1km:improved_grass_ha
                                                              , sugarbeet_ha, X, Y)
        , by = c("X", "Y"), all = T) %>%
  filter(!is.na(rcFid_1km)) %>%
  filter(!is.na(soilTyp))

# concatenate the three categories
cropSoilRain$snsCats <- paste(cropSoilRain$rainCat
                              , cropSoilRain$soilTyp
                              , cropSoilRain$depthRB209
                              , sep = "_") 
head(cropSoilRain)

# merge with SNS parameter data to make spatial
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

# use matrix calculations - just areas
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

# loop through sns categories
m3 <- list()
for(i in 0:3){
  if(i == 0){
    v = 1
  } else {
    v = v+1
  }
  
  print(i)
  nm <- paste0("sns", i, "_")
  
  # new instance
  snsNew <- snsMat2
  print(head(snsNew))
  # make 1 or 0
  M <- ifelse(snsNew==i,1,0)
  print(head(M))
  
  # multiply by the area
  print(head(snsMat1))
  m2 <- M * snsMat1
  print(head(m2))
  
  # sum the rows
  m3[[v]] <- m2 %>%
    as.data.frame() 
  colnames(m3[[v]]) <- paste0(nm, gsub("sns_", "", colnames(m3[[v]])))
  print(head(m3[[v]]))
}

# bind back together
snsParasWide3 <- do.call(cbind, list(snsParasWide2, m3)) %>% as.data.frame() %>%
  # merge back to spatial
  merge(rcFid.grid %>% dplyr::select(rcFid_1km), .)
head(snsParasWide3)

# and with the derived SNS values
snsParas <- merge(snsParasWide3, snsCropParametersWide 
                  , by = "threeParas"
                  , all = T) %>%
  dplyr::select(-c(contains(c("sns_"))))
head(snsParas)

# get all unique names
startCol <- which(colnames(snsParas) == "sns0_forage")
endCol <- which(colnames(snsParas) == "sns3_uncropped")
xCurr <- list()
stop("xxx")
for(i in startCol:endCol){
  
  cat(i, nm <- colnames(snsParas)[i], "\n")
  nmKg <- paste0(nm, "_kgN")
  
  # select only columns that match
  # multiple the amount of N suggested to be applied with the amount of hectares of
  # that crop type
  xCurr[[i]] <- snsParas %>%
    st_drop_geometry() %>%
    dplyr::select(rcFid_1km, contains(nm)) %>%
    mutate(total = .[, 2] * .[, 3]) %>%
    dplyr::select(rcFid_1km, total) %>%
    rename(!!nmKg := total)
}
xCurr <- Filter(length, xCurr) # only get filled objects

for(jx in 1:length(xCurr)){
  
  cat(jx, "\n")
  x1 <- xCurr[[jx]] %>% unnest(cols = c()) %>% as.data.frame() %>% filter(.[[2]] > 0)
  
  if(jx == 1){
    outout <- x1
  } else {
    outout <- merge(outout, x1, by = "rcFid_1km", all = T)
  }
  
}

xx <- merge(xCurr[[1]] %>% unnest(cols = c()) %>% as.data.frame()
            , xCurr[[2]] %>% unnest(cols = c()) %>% as.data.frame()
            , by = "rcFid_1km")
x1 <- xCurr[[1]] %>% unnest(cols = c()) %>% as.data.frame() %>% filter(.[[2]] > 0)
x2 <- xCurr[[2]] %>% unnest(cols = c()) %>% as.data.frame() %>% filter(.[[2]] > 0)

x1t <- table(x1$rcFid_1km) %>% as.data.frame()
x2t <- table(x2$rcFid_1km) %>% as.data.frame()
max(x1t$Freq)
max(x2t$Freq)

# merge all together
snsParasKgN <- Reduce(function(x, y) merge(unlist(x), unlist(y), by = "rcFid_1km", all = TRUE),
                      xCurr[1:2]) %>%
  # get total
  mutate(totalNperKm = rowSums(select(., 2:last_col()), na.rm = T)) %>%
  # merge back to spatial
  merge(cropSoilRain %>% dplyr::select(X, Y, rcFid_1km), .)
head(snsParasKgN)

# } # end of 'i'
