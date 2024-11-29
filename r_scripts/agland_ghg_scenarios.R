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
    
    ## convert to point data
    lcmDiff.point <- lcmDiff.change %>%
      st_as_sf(coords = c("x", "y"),
               crs = st_crs(lcm2015))
    lcmDiff.point <- lcmDiff.point %>%
      st_transform(., 27700)
    
    # determine which region each belongs to
    regions <- st_read("N:/Data/UK/Boundary/uk/eng_wales_scot.gpkg")
    
    if(file.exists("lcmDiff_region.gpkg")){
      lcmDiff.region <- st_read("lcmDiff_region.gpkg")
    } else {
      ## extract region to point data
      lcmDiff.region <- lcmDiff.point %>%
        ## intersect polygons with points, keeping the information from both
        st_intersection(regions, .)
      st_write(lcmDiff.region, "lcmDiff_region.gpkg", append = F)
      # lcmDiff.region <- st_read("lcmDiff_region.gpkg")
    }
    head(lcmDiff.region)
    
    ## ------------ Notes --------------  ##
    ## below, the finals of all animals, crops, etc. from the data inputs
    ## calculated in 'total_ghg_script.R' will be loaded in. Then, regional averages can 
    ## be used, and alterations made. 
    ## ------------ ----- --------------  ##
    
    # Load animals and crops data
    tghg.animals <- st_read("data_in/animals/all_animals_1km.gpkg")
    tghg.crops <- st_read("data_in/land_cover/land_cover_table.gpkg")
    
    # Combine animals and crops data using rcFid_1km
    tghg.ac <- tghg.animals %>%
      left_join(st_drop_geometry(tghg.crops), by = "rcFid_1km")
    
    # Load or create region data
    tghg.region <- if(file.exists("scenario/tghg_region.gpkg")) {
      st_read("scenario/tghg_region.gpkg")
    } else {
      tic("extraction to centroids")
      tghg.region <- tghg.ac %>%
        st_centroid() %>%
        st_intersection(regions, .)
      toc()
      st_write(tghg.region, "scenario/tghg_region.gpkg", append = FALSE)
      tghg.region
    }
    
    # Save data to disk
    st_write(tghg.animals, "scenario/tghg_animals.gpkg", append = FALSE)
    st_write(tghg.ac, "scenario/tghg_ac.gpkg", append = FALSE)
    st_write(tghg.crops, "scenario/tghg_crops.gpkg", append = FALSE)
    st_write(tghg.region, "scenario/tghg_region.gpkg", append = FALSE)
    
    # Merge regions with animal and crop numbers
    tghg.ac2 <- tghg.ac %>%
      left_join(st_drop_geometry(tghg.region) %>%
                  select(Name, Area_Description, rcFid_1km),
                by = "rcFid_1km") %>%
      relocate(rcFid_1km, Name, Area_Description) %>%
      distinct()
    
    # Load or create spatial data
    tghg.ac3 <- if(file.exists("ac3.gpkg")) {
      st_read("ac3.gpkg")
    } else {
      tghg.ac3 <- tghg.ac %>%
        select(rcFid_1km) %>%
        left_join(tghg.ac2, by = "rcFid_1km")
      st_write(tghg.ac3, "ac3.gpkg", append = FALSE)
      tghg.ac3
    }
    
    # Extract region to point data
    tic("extraction to centroids - region and change")
    tghg.ac4 <- st_intersection(lcmDiff.region %>% dplyr::select(-c(Name, Area_Description)),
                                tghg.ac3) %>%
      relocate(c(lcm2007, lcm2015, different), .after = rcFid_1km) %>%
      dplyr::select(-Area_Description)
    toc()
    
    # Summarize animals and crops by region and land type
    tghg.summarise <- tghg.ac4 %>%
      dplyr::select(-c(rcFid_1km, lcm2007, different)) %>%
      st_drop_geometry() %>%
      group_by(lcm2015, Name) %>%
      summarise(across(everything(), mean)) %>%
      mutate(n = n()) %>%
      relocate(lcm2015, Name, n)
    head(tghg.summarise)
    
    # Save summary data
    fwrite(tghg.summarise, "scenario/tghg_summarise.csv", row.names = FALSE)
    
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
  
  ##### calculate differences between 2015 and 2007 #####
  # read in 2007
  check2007 <- st_read(file.path("scenario", "scen_maps"
                         , "baseline2007.gpkg"))
  ## check land cover equal 100 ha (i.e., 100% area) [Sum columns ending with '_ha']
  landCheck <- check2007 %>%
    mutate(sum_ha = rowSums(across(ends_with('_ha')), na.rm = TRUE))
  hist(landCheck$sum_ha, breaks = 100)
  
  check2015 <- st_read("scenario/tghg_ac.gpkg")
  landCheck <- check2015 %>%
    mutate(sum_ha = rowSums(across(ends_with('_ha')), na.rm = TRUE))
  hist(landCheck$sum_ha, breaks = 100)
  
  #### 2 - create maps for the new scenarios ####
  ## ------------ Notes --------------  ##
  ## The new maps produced here will be based on four of the scenarios from Redhead et al. (2020)
  ## The scenarios are a 15 and 30 per cent increase in agricultural expansion, 
  ## and 15 and 30 per cent grassland restoration.
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
  # add the original 2015 map
  dfs2007 <- c(dfs2007, "scenario/tghg_ac.gpkg")
  
  cat(dfs2007, sep = "\n")
  stopifnot(length(dfs2007) == 6)
  
  # for each scenario, determine totals of a few variables (animals and land cover)
  dfs2007.amounts <- pblapply(dfs2007, function(i){
    
    print(i)
    
    # get name
    i.name <- gsub("^(.+)\\.gp.*", "\\1", basename(i))
    i.vect <- st_read(i, quiet = T)
    ## check land cover equal 100 ha (i.e., 100% area) [Sum columns ending with '_ha']
    landCheck <- i.vect %>%
      mutate(sum_ha = rowSums(across(ends_with('_ha')), na.rm = TRUE))
    print(hist(landCheck$sum_ha, breaks = 100))
    
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
    filter(Name %in% c("Medium.dairy_Dairy.cows"
                       , "poultry_per_km", "Continental_60...240.months_Beef.cows"
                       , "conifer_ha", "broadleaf_ha", "calc_grass_ha"
                       , "improved_grass_ha", "winterwheat_ha")) %>%
    relocate(Name, Ag_30)
  dfs2007.select
  
  # write readme
  sink(file = NULL)
  sink(file = file.path("scenario", "scen_maps", "readme_scenarios.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2024-01-30
Last update:",  format(Sys.Date()), "
Produced by 'agland_ghg_scenarios.R'

Description of scenario maps:
  Contains the region, land cover classification, animal numbers, and land cover areas for each 1 km2. The values were derived by:
    the average 2015 values for that region (county), land cover classification combination was used.
    
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