## ---------------------------
##
## Script name: agland_ghg_scenarios.R
##
## Purpose of script: run the GHG model on different scenarios of land cover.
##                    The land covers were derived from previously-published
##                    values
##                    
## ------------ Notes --------------  ##
## ??
## ??
## ------------ ----- --------------  ##
##
## Run after: total_ghg_script.R
##
## Run before: ??
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
  "terra", "stars"
  , "dplyr", "sf"
  , "tictoc"
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
## in this section you select which bits of the script you would ##
## like to run this time. Please note: if this is the first time ##
## the script has been run, you will need to do it sequentially ##

#### 0 - load functions ####

#### 1 - check different between land cover maps 2007 and 2015 ####
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
  51133/241797
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
          , by = "rcFid_1km") %>%
    # remove original region
    dplyr::select(-region)
  
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
  ### merge regions with animal and crop numbers
  tghg.ac2 <- tghg.ac %>% st_drop_geometry() %>%
    merge(., tghg.region %>% dplyr::select(Name, Area_Description, rcFid_1km) %>% st_drop_geometry()
          , by = "rcFid_1km"
          , all = F) %>%
    relocate(rcFid_1km, Name, Area_Description) %>%
    # remove duplicate rows
    distinct()
  st_write(tghg.ac2, "ac2.gpkg", append = T)
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
    summarise(across(everything(), ~mean(.)))
  
  # create the 2007 dataframe based on either:
  # averages for the region / land cover type from 2015 where the land cover has changed, or
  # or same as 2015 if the land cover type has not changed. 
  # For this, use the 'different' column
  
  df2007.changed <- tghg.ac4 %>%
    st_drop_geometry() %>%
    filter(different == 1) %>%
    # keep necessary columns
    dplyr::select(c(rcFid_1km, lcm2007, Name)) %>%
    # merge with summary, to get average animals and crops (from 2015 map)
    merge(., tghg.summarise %>% ungroup()
          , by.x = c("Name", "lcm2007")
          , by.y = c("Name", "lcm2015")
          , all = T) %>%
    filter(!is.na(rcFid_1km)) %>%
    mutate(across(everything(), ~replace_na(., 0)))
  
  df2007.noChange2015 <- tghg.ac4 %>%
    st_drop_geometry() %>%
    filter(different == 0)
  
  ## combine back to create 2007 starting dataset
  df2007.final <- df2007.changed %>%
    bind_rows(., df2007.noChange2015 %>% dplyr::select(-c(lcm2015, different))) %>%
    merge(tghg.ac3 %>% dplyr::select(rcFid_1km)
          , .)
  names(df2007.final)
  
  ## convert to raster and save
  df2007.raster <- df2007.final %>%
    dplyr::select(lcm2007) %>%
    st_centroid() %>%
    st_rasterize()
  
  ### save 
  st_write(df2007.final
           , file.path("scenario", "scen_maps"
                       , "baseline2007.gpkg"), append = T)
  write_stars(df2007.raster
              , file.path("scenario", "scen_maps"
                          , "baseline2007.tif")
              , overwrite = T)
  
  # write readme
  sink(file = NULL)
  sink(file = file.path("scenario", "scen_maps", "readme_baseline2007.md"))
  cat("Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2024-01-30
Last update:",  format(Sys.Date()), "
Produced by 'agland_ghg_scenarios.R'

Description of 'baseline2007.gpkg':
  Contains the region, land cover classification, animal numbers, and land cover areas for each 1 km2. The values were derived in one of two ways:
    1. If the land cover classification was the same in 2007 as 2015, the values were taken from the 2015 dataset (see full project)
    2. If the land cover classification was the different in 2007 compared to 2015, the average 2015 values for that region, land cover classification combination was used.

Columns of 'baseline2007.gpkg':
  'rcFid_1km' = centroid point of 1 km reference id (EPSG: 27700) [x and y coordinates seperated by '_']
  'Name' = region of the the UK that the polygon centroid was in.  
  'lcm2007' = the dominant land cover classification that that 1 km2 polygon was assigned to in the land cover map 2007
  '[animal_name]' = how many of the animal type in each column name was found in that polygon. See main project for descriptions
  '[land_cover]_ha' = area, in hectares, of that land cover found in that 1 km2 pixel, as determined (in the original project) by using the 25 m2 land cover raster
  
Description of 'baseline2007.tif':
  Contains the land cover classification of each pixel. 

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
  as.data.frame(., xy = T)
names(scenMaps)
## keep only important colnames info
names(scenMaps) <- gsub(".*riolu(.+)_af.*", "\\1", names(scenMaps))
names(scenMaps)

### determine which change from the baseline 2007 one
df2007 <- df2007.raster %>% as.data.frame(xy = T) %>%
  mutate(x = round(x, 0)
         , y = round(y, 0)) %>%
  merge(., scenMaps %>%
          mutate(x = round(x, 0)
                 , y = round(y, 0))
  ) %>%
  filter(baseline2007 > 0)

# see which pixels changed
df2007.changes <- df2007 %>%
  mutate(across(Ag_15:Gr_30, ~ifelse(. != baseline2007, 1, 0)
                , .names = "{.col}_changes")) %>%
  # make rcFid_1km, to match other dfs
  mutate(rcFid_1km = paste(x, y, sep = "_")) %>%
  ## use id to match up regions
  merge(., df2007.final %>% st_drop_geometry() %>% dplyr::select(rcFid_1km, Name)) %>%
  relocate(rcFid_1km, Name)

# get averages of all animals and crops per region per land type
df2007.summarise <- df2007.final %>%
  # remove unneeded columns for summary
  dplyr::select(-c(rcFid_1km)) %>%
  st_drop_geometry() %>%
  group_by(lcm2007, Name) %>%
  summarise(across(everything(), ~mean(.)))

# for each scenario, determine which pixels are different and the same, 
# and repeat the methodology in section 1 to get full dataset
for(i in names(df2007)[4:7]){
  print(i)
  
  # use 'i' to extract columns for scenario
  df2007.scen <- df2007.changes %>%
    dplyr::select(c(rcFid_1km:baseline2007
                    , contains(i))) %>%
    rename(lc = 6, change = 7)
  
  table(df2007.scen$change, useNA = "always")
  
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
}

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

#### 3 - run the GHG model on the new scenarios (and baseline for 2007) ####
## ------------ Notes --------------  ##
## below follows the exact calculations from the main ghg script 'total_ghg_script.R',
## with only some alterations for reading in different sources (i.e. scenarios)
## It starts from part 4 of that script
## ------------ ----- --------------  ##

# list all final files that will be assessed
mapsToEval.list <- list.files(file.path("scenario", "scen_maps")
                              , full.names = T
                              , pattern = ".gpkg")



