## ---------------------------
##
## Script name: nonBovines.R
##
## Purpose of script: create maps for non-Bovines in the GHG tool, for Scotland
##                    and Wales
##                    
## ------------ Notes --------------  ##
## This uses maps that were georeferenced in ArcGIS Pro.
## It includes the calculation of pigs, sheep, and poultry
## ------------ ----- --------------  ##
##
## Run after: ghgtool_setup.R
##
## Run before: total_ghg_script.R
##
## Specific numbered tasks:
## 1 - 
## 2 - 
## 3 - 
## 4 - 
##
## list of final outputs:
##    
##    
##
## Author: Dr. Paul M. Evans
##
## Date Created: 2024-02-12
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
  "terra", "stars", "mapview"
  , "dplyr", "sf", "ggplot2"
  , "tictoc", "purrr", "tidyr"
  , "parallel", "doParallel", "foreach" # parallel loops
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
# scotland and wales boundaries
swPath <- "C:/Users/paueva/OneDrive - UKCEH/Data/uk/boundary"
# animal path
animalPath <- "C:/Users/paueva/Documents/ArcGIS/Projects/agland2"

#### 0 - load functions ####

#### 1 - dissolve scotland and wales ####
swShape <- st_read(file.path(swPath, "sw.gpkg"))
## split into wales and scotland
unique(swShape$Area_Code)
# the summarize function is based on the one from dplyr 
swSplit <- swShape %>%
  group_split(Area_Code) 
sw <- lapply(1:2, function(x){
  sw1 <- swSplit[[x]] %>%
    st_cast() %>%
    st_union()
  plot(sw1[1])
  return(sw1)
})

# get ymin of scotland
scotMin <- st_bbox(sw[[1]])[[2]]

# get smaller region of...
## wales
walesRegions <- st_read("N:/Data/UK/Boundary/wales/Wales_Ward_Boundaries_hwm/Wales_Ward_Boundaries_hwm.shp") %>%
  st_transform(27700)
head(walesRegions)
### group by 'FILE_NAME'
walesReg2 <- walesRegions %>% group_by(FILE_NAME) %>% summarize() 
plot(walesReg2)
unique(walesReg2$FILE_NAME)
#### determine in terms of regions from https://statswales.gov.wales/Catalogue/Agriculture/Agricultural-Survey/Area-Survey-Results/total-livestock-in-wales-by-area 
walesReg2 <- walesReg2 %>%
  mutate(broad_file_name = stringr::str_to_title(gsub(".*-_(.+)$", "\\1", FILE_NAME))) %>%
  dplyr::select(-FILE_NAME) %>%
  mutate(broad_region = ifelse(grepl("nglesey|wynedd|onwy|enbighshire|lintshire|rexham|wansea|talbot|ridgend|lamorgan|cynon_taf|ydfil|aerphilly|gwent|orfaen|onmouthshire|ewport|ardiff|owys|eredigion|embrokeshire|armarthenshire", broad_file_name)
                               , "broad_region", "???"))

## Scotland
scotRegions <- st_read("N:/Data/UK/Boundary/Counties/Data/Supplementary_Ceremonial/Boundary-line-ceremonial-counties_region.shp")  %>%
  st_transform(27700)
head(scotRegions)
### get only scottish counties
scotReg <- scotRegions %>%
  filter(NAME %in% c("Tweeddale", "Clackmannan", "Midlothian", "Lanarkshire"
                     , "City of Glasgow", "City of Aberdeen", "City of Dundee"                       
                     , "City of Edinburgh", "Roxburgh, Ettrick and Lauderdale"
                     , "Ross and Cromarty", "Fife", "Aberdeenshire", "Angus"                                
                     ,"Argyll and Bute", "Ayrshire and Arran", "Dumfries"                   
                     , "Banffshire", "Berwickshire", "Caithness", "Dunbartonshire"
                     , "East Lothian", "Inverness", "Kincardineshire", "Nairn"
                     , "Moray", "Orkney", "Perth and Kinross"
                     , "Renfrewshire", "Ross and Cromarty", "Roxburgh, Ettrick and Lauderdale"
                     , "Shetland", "Stirling and Falkirk", "Sutherland"
                     , "The Stewartry of Kirkcudbright", "West Lothian", "Wigtown"))
plot(scotReg[1])

#### 2 - load in animal amounts ####
# pigs
pigs.in <- st_read(file.path(animalPath, "pigs.shp")) %>% st_transform(27700)
## separate between scotland and wales using coordinates
pigs.ymax <- as.data.frame(matrix(unlist(lapply(st_geometry(pigs.in), st_bbox))
                                  , nrow = nrow(pigs.in)
                                  , byrow = T)) %>%
  dplyr::select(4) %>%
  mutate(sw = ifelse(.[[1]] < scotMin, "Wales", "Scot"))
pigs.in <- bind_cols(pigs_per_km = pigs.in %>% dplyr::select(pigs_km)
                     , sw = pigs.ymax$sw) %>% rename(pigs_per_km = 1)

# sheep
sheep.in <- st_read(file.path(animalPath, "sheep.shp")) %>% st_transform(27700)
## separate between scotland and wales using coordinates
sheep.ymax <- as.data.frame(matrix(unlist(lapply(st_geometry(sheep.in), st_bbox))
                                   , nrow = nrow(sheep.in)
                                   , byrow = T)) %>%
  dplyr::select(4) %>%
  mutate(sw = ifelse(.[[1]] < scotMin, "Wales", "Scot"))
sheep.in <- bind_cols(sheep_per_km = sheep.in %>% dplyr::select(sheep_numb)
                      , sw = sheep.ymax$sw) %>% rename(sheep_per_km = 1)

# poultry
poultry.in <- st_read(file.path(animalPath, "poultrySW.shp")) %>% st_transform(27700)
## separate between scotland and wales using coordinates
poultry.ymax <- as.data.frame(matrix(unlist(lapply(st_geometry(poultry.in), st_bbox))
                                     , nrow = nrow(poultry.in)
                                     , byrow = T)) %>%
  dplyr::select(4) %>%
  mutate(sw = ifelse(.[[1]] < scotMin, "Wales", "Scot"))
poultry.in <- bind_cols(poultry_per_km = poultry.in %>% dplyr::select(poul_numb)
                        , sw = poultry.ymax$sw) %>% rename(poultry_per_km = 1)

# tidy
rm(pigs.ymax, sheep.ymax, poultry.ymax, scotMin)

#### 3 - get difference in area between nations and animal shapes ####
# get outlines
scotland.outline <- sw[[1]] %>% st_as_sf()
plot(scotland.outline[1])
wales.outline <- sw[[2]] %>% st_as_sf()
plot(wales.outline[1])

# list animals, by nation
animals.nms <- c("pig", "sheep", "poultry")
animals.list <- list(pigs.in, sheep.in, poultry.in)

## scotland
scotland.animals <- list(pigs.in %>% filter(sw == "Scot")
                         , sheep.in %>% filter(sw == "Scot")
                         , poultry.in %>% filter(sw == "Scot"))
## wales
wales.animals <- list(pigs.in %>% filter(sw == "Wales")
                      , sheep.in %>% filter(sw == "Wales")
                      , poultry.in %>% filter(sw == "Wales"))

### save animals and outlines
for(i in 1:length(animals.list)){
  iFilter <- animals.list[[i]] %>% filter(sw == "Wales")
  st_write(iFilter, file.path(swPath, "pps", paste0(animals.nms[[i]], "_wales.gpkg"))
           , append=F)
  iFilter <- animals.list[[i]] %>% filter(sw == "Scot")
  st_write(iFilter, file.path(swPath, "pps", paste0(animals.nms[[i]], "_scot.gpkg"))
           , append=F)
}
#### save wales
st_write(wales.outline, file.path(swPath, "pps", "wales_outline.gpkg")
         , append=F)
#### save scotland
st_write(scotland.outline, file.path(swPath, "pps", "scot_outline.gpkg")
         , append=F)
## ------------ Notes --------------  ##
## creating differences in area were done in QGIS. It used the 'native:difference'
## function. Example command:  processing.run("native:difference",
## {'INPUT':'.../uk/boundary/pps/wales_outline.gpkg|layername=wales_outline',
## 'OVERLAY':'.../uk/boundary/pps/sheep_wales.gpkg|layername=sheep_wales',
## 'OUTPUT':'sheep_wales_diff.gpkg','GRID_SIZE':None})

## following, union between the differences and original numbers were calculated
## this was also done in QGIS. Example: processing.run("native:union",
## {'INPUT':'.../uk/boundary/pps/pig_scot.gpkg|layername=pig_scot',
## 'OVERLAY':'.../uk/boundary/nonB/pig_scot_diff.gpkg',
## 'OVERLAY_FIELDS_PREFIX':'','OUTPUT':'...uk/boundary/pps/qgis/pig_scot_union.gpkg',
## 'GRID_SIZE':None})

## Below, all the differences from that process are loaded in
## ------------ ----- --------------  ##

#### 4 - get animal numbers ####
## ------------ Pig and Sheep Notes --------------  ##
## livestock data for Scotland were obtained by using 1 km2 published maps (in image form)
## from UK Government data. The original categories were ranges of:
## <=1, 1-5, 5 - 20, 20-80, 80-200, 200-430 for pigs
## <=1, 1-30, 30-60, 60-120, 120-240, 240-460 for sheep
## the lowest number of each bracket was used 
## These data were geoprocessed from the original, polygonised.
## Below, these polygons will rasterised 
## ------------ ----- ----------------------------  ##

## ------------ Poultry Notes --------------  ##
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
## ------------ ------------- --------------  ##

# list the shapefiles that were unioned
union.list <- list.files(file.path(swPath, "pps", "qgis")
                         , pattern = "_union.gpkg"
                         , full.names = T)
union.list
stopifnot(length(union.list) == 6)

## for each, load and assign the correct final value, and tidy up the df
i.list <- list()
for(i in 1:length(union.list)){
  x <- union.list[[i]]
  xAn <- gsub("^(.+)_s.*|_w.*", "\\1", basename(x))
  xAnPer <- paste0(xAn, "_per_km")
  xReg <- ifelse(grepl("wale", x), "wales", "scotland")
  cat(x, paste0("\n", xAn, " | ", xReg), "\n")
  df.in <- st_read(x, quiet = T)
  
  if(grepl("pig_scot_union|pig_wales_union", x)){
    df.in[is.na(df.in$pigs_per_km), "pigs_per_km"] <- 1
  } else if(grepl("poul_scot_union|poul_wales_union", x)){
    df.in[is.na(df.in$poul_numb), "poul_numb"] <- 50
  } else if(grepl("sheep_scot_union|sheep_wales_union", x)){
    df.in[is.na(df.in$sheep_numb), "sheep_numb"] <- 30
  }
  
  # keep only the animals per km column
  df.in2 <- df.in %>%
    dplyr::select(any_of(c("sheep_numb", "pigs_per_km", "poul_numb"))) %>%
    rename(!!xAnPer := 1) %>%
    # convert to 27700
    st_transform(27700) %>%
    mutate(area = st_area(.)) %>%
    mutate(km_area = as.vector(area/1000000)) %>%
    dplyr::select(-area)
  
  # stop if NAs still exist
  if(anyNA(df.in2[, 1])){stop("NAs here")}
  
  # get centroids, to be able to see the specific region of Wales / Scotland
  df.cent <- df.in2 %>%
    st_centroid()
  ## extract the county info
  if(grepl("wales", x)){
    # intersect polygons with points, keeping the information from both
    df.county <- st_intersection(walesReg2, df.in2 %>% st_make_valid()) %>%
      rename(county = broad_file_name) %>%
      dplyr::select(-broad_region) %>%
      mutate(area = st_area(.)) %>%
      mutate(km_area = as.vector(area/1000000)) %>%
      dplyr::select(-area)
    # ggplot() +
    #   geom_sf(data = df.county, col = "blue")
    # ggplot() +
    #   geom_sf(data = walesReg2, col = "red")
  } else if(grepl("scot", x)){
    df.county <- st_intersection(scotReg, df.in2) %>%
      rename(county = NAME) %>%
      dplyr::select(-DESCRIPTIO) %>%
      mutate(area = st_area(.)) %>%
      mutate(km_area = as.vector(area/1000000)) %>%
      dplyr::select(-area)
    # ggplot() +
    #   geom_sf(data = scotReg, col = "red")
  }
  # str(df.county)
  
  # need to determine the most accurate numbers. This was done using
  # published government statistics, and the maps 
  ## ------------ Notes --------------  ##
  ## livestock data for Wales can be obtained in two ways: total livestock
  ## numbers in 2017, by region (https://statswales.gov.wales/Catalogue/Agriculture/Agricultural-Survey/Area-Survey-Results/total-livestock-in-wales-by-area)
  ## and 2020 livestock numbers, by animal type and 'role' (https://www.gov.wales/survey-agriculture-and-horticulture-june-2022)
  ## From the 2017 data, the proportion of animals in different regions of Wales could have been calculated in terms
  ## of percentage. However, to make the methodology consistent, they were calculated the same way as Scotland. See next 'Notes' section
  ## ------------ ----- --------------  ##
  
  # convert top values to middle and lowest values
  if(xAn == "sheep"){
    sheepNumb <- df.numb <- df.county %>%
      mutate(sheep_per_km_high = sheep_per_km
             , sheep_per_km_mid = ifelse(sheep_per_km == 460, 350
                                         , ifelse(sheep_per_km == 240, 180
                                                  , ifelse(sheep_per_km == 120, 90
                                                           , ifelse(sheep_per_km == 60, 45
                                                                    , ifelse(sheep_per_km == 30, 15
                                                                             , ifelse(sheep_per_km == 1, 1, 0
                                                                             ))))))
             , sheep_per_km_low = ifelse(sheep_per_km == 460, 240
                                         , ifelse(sheep_per_km == 240, 120
                                                  , ifelse(sheep_per_km == 120, 60
                                                           , ifelse(sheep_per_km == 60, 30
                                                                    , ifelse(sheep_per_km == 30, 1
                                                                             
                                                                             , ifelse(sheep_per_km == 1, 0, 0
                                                                             ))))))) %>%
      mutate(total_high = sheep_per_km_high * km_area
             , total_mid = sheep_per_km_mid * km_area
             , total_low = sheep_per_km_low * km_area)
    
  } else if(xAn == "poul") {
    
    poultryNumb <- df.numb <- df.county %>%
      mutate(poul_per_km = as.numeric(poul_per_km)) %>%
      mutate(poul_per_km_high = poul_per_km
             , poul_per_km_mid = ifelse(poul_per_km == 39000, 20988
                                        , ifelse(poul_per_km == 2500, 1900
                                                 , ifelse(poul_per_km == 1300, 950
                                                          , ifelse(poul_per_km == 600, 425
                                                                   , ifelse(poul_per_km == 250, 150
                                                                            , ifelse(poul_per_km == 50, 25, 0
                                                                            ))))))
             , poul_per_km_low = ifelse(poul_per_km == 39000, 2501
                                        , ifelse(poul_per_km == 2500, 1301
                                                 , ifelse(poul_per_km == 1300, 601
                                                          , ifelse(poul_per_km == 600, 251
                                                                   , ifelse(poul_per_km == 250, 51
                                                                            
                                                                            , ifelse(poul_per_km == 50, 1, 0
                                                                            ))))))) %>%
      mutate(total_high = poul_per_km_high * km_area
             , total_mid = poul_per_km_mid * km_area
             , total_low = poul_per_km_low * km_area)
    
  } else if(xAn == "pig"){
    
    pigNumb <- df.numb <- df.county %>%
      mutate(pig_per_km_high = pig_per_km
             , pig_per_km_mid = ifelse(pig_per_km == 430, 315
                                       , ifelse(pig_per_km == 200, 140
                                                , ifelse(pig_per_km == 80, 50
                                                         , ifelse(pig_per_km == 5, 3
                                                                  , ifelse(pig_per_km == 250, 150
                                                                           , ifelse(pig_per_km == 1, 1, 0
                                                                           ))))))
             , pig_per_km_low = ifelse(pig_per_km == 430, 200
                                       , ifelse(pig_per_km == 200, 80
                                                , ifelse(pig_per_km == 80, 20
                                                         , ifelse(pig_per_km == 5, 1
                                                                  , ifelse(pig_per_km == 1, 0, 0)))))) %>%
      mutate(total_high = pig_per_km_high * km_area
             , total_mid = pig_per_km_mid * km_area
             , total_low = pig_per_km_low * km_area)
    
  }
  
  # get totals of animals
  xSums <- df.numb %>%
    group_by(county) %>%
    summarise(sumHigh = sum(total_high)
              , sumMid = sum(total_mid)
              , sumLow = sum(total_low)) %>% st_drop_geometry() %>%
    mutate(animal = xAn)
 

  i.list[[i]] <- xSums
  
} # end 'i'

# merge all together
i.bound <- bind_rows(i.list) %>%
  pivot_wider(names_from = animal, values_from = c(sumHigh, sumMid, sumLow))
## save
fwrite(i.bound, file.path("data_in", "animals", "anmials_area.csv")
       , row.names = F)
