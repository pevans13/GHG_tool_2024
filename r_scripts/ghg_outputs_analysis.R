## ---------------------------
##
## Script name: ghg_outputs_analysis.R
##
## Purpose of script: to create the final figures to use in the paper of the model
##                    
## ------------ Notes --------------  ##
## The figures will include the 2015 data, and also the 
## scenarios that were used. See Redhead et al. (2020). 

## The final maps were already produced in the 'total_ghg_script.R' script.
## ------------ ----- --------------  ##
##
## Run after: total_ghg_script.R
##
## Specific figures:
## 1 - total GHG emissions for all scenarios combined 
## 2 - proportion of emissions (based on origin) assigned to each scenario
## 3 - proportions of GB hectads with different levels of percentage change
##
## Author: Dr. Paul M. Evans
##
## Date Created: 2024-07-23
##
## Copyright (c) Paul M. Evans, 2023
## Email: paueva@ceh.ac.uk
## Git: https://github.com/pevans13

## ---------------------------
options(scipen = 6, digits = 4) # for non-scientific notation
## ---------------------------

# remove objects
rm(list = ls())

#### 0 - paths/ names ####
spatialGrid <- "agland_grid_ESW.gpkg"

#### 0 - load libraries ####
source("./r_scripts/proj_packages.R")

#### 0 - load functions ####

#### Table 1 - total land cover ####
# load 2015 land cover map
currPath <- file.path("./data_in", "land_cover")
# load 25 m LCM raster
lcmCrops25 <- rast(file.path(currPath, "LCM2015PlusCrops.tif"))
lcTable <- table(lcmCrops25 %>% as.vector()) %>% as.data.frame()
lcTableArea <- lcTable %>%
  mutate(m2 = 625 * Freq
         , km2 = m2 / 1000000)
sum(lcTableArea$km2)
# calculate percentage cover
lcTableArea$pc_cover <- (lcTableArea$km2 / sum(lcTableArea$km2, na.rm = TRUE) * 100)
# for ag
lcTableArea %>%
  filter(startsWith(as.character(Var1), "3")) %>%
  summarise(total_value = sum(pc_cover))
lcTableArea %>%
  filter(startsWith(as.character(Var1), "3")) %>%
  summarise(total_value = sum(km2))

#### Fig 1 - total GHG emissions for all scenarios combined  ####
##### 1a - load in all final csvs #####
## ------------ Notes --------------  ##
## add a new column to all, signifying the scenario
## ------------ ----- --------------  ##

# original data from 2015
data2015 <- fread(file.path("results", "ghg_final_emissions_CO2e.csv")) %>%
  as.data.frame() %>%
  # add new column
  mutate(data = "d2015")

# data from the scenarios
scenarios <- list.files(file.path("scenario", "final_results")
                        , pattern = "_final_emissions_CO2e(.+).csv"
                        , full.names = T)
dataScenarios <- lapply(scenarios, function(x) { 
  # get name
  xn <- gsub(".*CO2e(.+).csv.*", "\\1", basename(x))
  # read in
  fread(x) %>% as.data.frame() %>%
    # add new column
    mutate(data = xn)})

# for original data and scenarios, combine
allData <- bind_rows(dataScenarios) %>%
  bind_rows(data2015, .)

# tidy
rm(data2015, dataScenarios)

##### 1b - summary data #####
## calculate the final values in MtCO2e
allDataSum <- allData %>%
  rowwise() %>%
  # sum without the minuses
  mutate(pure_emissions_Tco2e = sum(kgCO2e_entferm, kgCO2e_entferm, kgCO2e_manMgmt
                                    , kgCO2e_excrete, fertliser_CO2e_kgkm2, totalResid_kgKm2
                                    , fuel_kgco2e)) %>%
  group_by(data) %>%
  summarise(MtCO2e_balance = (sum(total_Tco2e)/1000000)
            , MtCO2e_emissions = ((sum(pure_emissions_Tco2e)/1000)/1000000)
            , n = n())

## calculate all in MtCO2e and sum in long form
allDataLong <- allData %>%
  mutate(
    across(
      contains("kg"),  # Select columns where the name contains 'kg'
      ~ . / 1000       # divide values by 1000 to get tonnes
      / 1000000        # divide values by 1000000 to get million tonnes
    )) %>%
  rename_with(
    ~ gsub("kg", "Mt", .),  # Replace 'kg' with 'Mt' in column names
    contains("kg")         # Apply only to columns with 'kg' in the name
  ) %>%
  dplyr::select(-rcFid_1km) %>%
  pivot_longer(cols = !c(data, total_Mtco2e, total_Tco2e, total_Tco2e_inv)) %>%
  group_by(data, name) %>%
  summarise(summed_MtCO2e = sum(value))
## save, for exploration
fwrite(allDataLong, "results/allDataLong.csv", row.names = F)

## sum all data based on scenarios
allDataVarSum <- allData %>%
  mutate(
    across(
      contains("kg"),  # Select columns where the name contains 'kg'
      ~ . / 1000       # divide values by 1000 to get tonnes
      / 1000000        # divide values by 1000000 to get million tonnes
    )) %>%
  rename_with(
    ~ gsub("kg", "Mt", .),  # Replace 'kg' with 'Mt' in column names
    contains("kg")         # Apply only to columns with 'kg' in the name
  ) %>%
  dplyr::select(-c(matches("Tco", ignore.case = F))) %>%
  group_by(data) %>%
  # sum all columns
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

# create a long er version, for the purposes of plotting
allDataVarLong <- allDataVarSum %>%
  dplyr::select(-total_Mtco2e) %>%
  pivot_longer(!data) %>%
  mutate(dataFac = as.factor(data))
## make wide, to get amount per category
allDataVarWide <- allDataVarLong %>%
  pivot_wider(names_from = "name", values_from = "value")


## add extra 'blank' bar to the long version
add2 <- allDataVarLong %>%
  filter(data == "d2015") %>%
  mutate(across(everything(), ~ NA)) %>%
  mutate(dataFac = "blank")

# Combine data and assign factor levels
allDataVarLongEx <- allDataVarLong %>%
  bind_rows(add2) %>%
  mutate(dataFac = factor(
    dataFac,
    # convert data factors, to separate 2015 and 2007 data
    levels = c("d2015", "blank", "Gr_30", "Gr_15", "baseline2007", "Ag_15", "Ag_30")
  ))

barOut <- ggplot(allDataVarLongEx, aes(x = dataFac, y = value, fill = name)) +
  # geom_col(colour = "black", position = "fill") +  scale_y_continuous(labels = scales::percent) +
  geom_col(colour = "black", , position = "stack") +  scale_y_continuous() +
  scale_fill_brewer(palette = "Pastel1"
                    , na.translate = FALSE
                    , labels = c("MtCO2e_entferm" = "Enteric"
                                 , "MtCO2e_manMgmt" = "Manure"
                                 , "MtCO2e_excrete" = "Excretion"
                                 , "fertliser_CO2e_Mtkm2" = "Fertilier app."
                                 , "landCoverTotalEmis_Mtco2" = "Land use (non-ag)"
                                 , "totalResid_MtKm2" = "Residues"
                                 , "fuel_Mtco2e" = "On-farm energy")) +
  annotate("text", x = 3.5, y = -27, label = "Grassland expansion", size = 3.5, hjust = 0.5) +
  annotate("text", x = 6.5, y = -27, label = "Agricultural expansion", size = 3.5, hjust = 0.5) +
  annotate("text", x = 1, y = -27, label = "2015", size = 3.5, hjust = 0.5) +
  annotate("text", x = 5, y = -27, label = "2007", size = 3.5, hjust = 0.5) +
  labs(y = expression("Mt CO"[2] * "e") # CO2e with subscript
       , fill = "Source") +
  # Rename ticks
  scale_x_discrete(labels = c("d2015" = "Baseline", "Gr_30" = "30%", "Gr_15" = "15%"
                              , "baseline2007" = "Baseline", "Ag_30" = "30%", "Ag_15" = "15%"
                              , "blank" = "")) +
  scale_y_continuous(breaks = seq(-25, 100, by = 5)) +
  theme_light() + theme(axis.title.x = element_blank()
                        , panel.grid.major.x = element_blank()
                        , panel.grid.minor = element_blank()
                        # remove outside line
                        , panel.border = element_blank()
                        # remove x ticks
                        , axis.ticks.x = element_blank()
                        , axis.text.x = element_text(vjust = 5))

png(file.path("scenario", "final_results"
              , paste0("stack_scen_compare"
                       , Sys.Date(), ".png"))
    , width = 1500, height = 900, units = "px"
    , res = 200
    )
plot(barOut)
dev.off()


ggplot(allDataVarLongEx, aes(x = dataFac, y = value, fill = name)) +
  geom_col(colour = "black", position = "fill") +  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Pastel1")  +
  theme_minimal() +
  annotate(
    "text", x = 3.5, y = -1.1, label = "Grassland expansion", size = 3.5, hjust = 0.5
  ) +
  annotate(
    "text", x = 6.5, y = -1.1, label = "Agricultural expansion", size = 3.5, hjust = 0.5
  ) +
  labs(y = expression("Mt CO"[2] * "e") # CO2e with subscript
  ) +
  # Rename ticks
  scale_x_discrete(labels = c("d2015" = "2015", "Gr_30" = "30%", "Gr_15" = "15%"
                              , "baseline2007" = "Baseline", "Ag_30" = "30%", "Ag_15" = "15%"
                              , "blank" = "")) +
  theme(axis.title.x = element_blank())

##### 1c - compare animal and land covers #####
# original data
animals2015 <- st_read("data_in/animals/all_animals_1km.gpkg") %>% st_drop_geometry() %>%
  # add new column
  mutate(data = "d2015")
land2015 <- st_read(file.path("data_in", "land_cover", "land_cover_area.gpkg")) %>% st_drop_geometry() %>%
  # add new column
  mutate(data = "d2015")

# data from the scenarios
## animals
animalScenarios <- list.files(file.path("scenario", "scen_maps", "finer_detail")
                              , pattern = "_animal"
                              , full.names = T)
animalScenarios <- lapply(animalScenarios, function(x) { 
  # get name
  xn <- gsub("^(.+)_animal.*", "\\1", basename(x))
  # read in
  st_read(x) %>% st_drop_geometry() %>%
    # add new column
    mutate(data = xn)}) %>% bind_rows()

## land
landScenarios <- list.files(file.path("scenario", "scen_maps", "finer_detail")
                            , pattern = "_land"
                            , full.names = T)
landScenarios <- lapply(landScenarios, function(x) { 
  # get name
  xn <- gsub("^(.+)_land.*", "\\1", basename(x))
  # read in
  st_read(x) %>% st_drop_geometry() %>%
    # add new column
    mutate(data = xn)}) %>% bind_rows() %>%
  # get all arable
  mutate(Arable_ha = rowSums(across(winterwheat_ha:sugarbeet_ha)))

# for original data and scenarios, combine
allDataLand <- bind_rows(land2015, landScenarios)
allDataAnimals <- bind_rows(animals2015, animalScenarios)
# tidy
rm(land2015, landScenarios, animals2015, animalScenarios)
gc()

# summarise to see comparisons
allDataLandLong <- allDataLand %>%
  dplyr::select(-rcFid_1km) %>%
  pivot_longer(cols = !data) %>%
  group_by(data, name) %>%
  summarise(area_ha = sum(value, na.rm = T))
## and again
land2 <- allDataLandLong %>%
  group_by(data) %>%
  summarise(area_ha_total = sum(area_ha, na.rm = T))
## save, for exploration
fwrite(allDataLandLong, "results/trials/allDataLandLong.csv", row.names = F)

# summarise to see comparisons
allDataAnimalLong <- allDataAnimals %>%
  dplyr::select(-rcFid_1km) %>%
  pivot_longer(cols = !data) %>%
  group_by(data, name) %>%
  summarise(animal_number = sum(value, na.rm = T))
head(allDataAnimalLong)
## and again
anim2 <- allDataAnimalLong %>%
  group_by(data) %>%
  summarise(animal_number_total = sum(animal_number, na.rm = T))
## save, for exploration
fwrite(allDataAnimalLong, "results/trials/allDataAnimalLong.csv", row.names = F)
gc()

# for each category, see the difference in distirbution
w1 <- which(names(allDataAnimals) == "Continental_0...3.months_Breeding.Heifers")
w2 <- which(names(allDataAnimals) == "rams.for.service..c41.")

pdf(file.path("investigate_differences.pdf"))
for(ww in w1:w2){
  
  # get name
  print(nowName <- names(allDataAnimals)[ww])
  
  # select just that name
  allDataAnimalsCut <- allDataAnimals %>%
    dplyr::select(c(all_of(nowName), data)) %>%
    rename(evalCol = 1)
  ## create boxplot
  dfbox <- ggplot() + 
    geom_boxplot(data = allDataAnimalsCut
                 , aes(x=as.factor(data), y=evalCol)
                 , position=position_dodge(1)) +
    ggtitle(nowName) + theme_bw()
  print(dfbox)
  Sys.sleep(1)
}
dev.off()

allDataLand2 <- allDataLand %>%
  relocate(rcFid_1km, data)
# for each category, see the difference in distirbution
w1 <- which(names(allDataLand2) == "broadleaf_ha")
w2 <- which(names(allDataLand2) == "sugarbeet_ha")

pdf(file.path("investigate_differences_land.pdf"))
for(ww in w1:w2){
  
  # get name
  print(nowName <- names(allDataLand2)[ww])
  
  # select just that name
  allDataLand2Cut <- allDataLand2 %>%
    dplyr::select(c(all_of(nowName), data)) %>%
    rename(evalCol = 1)
  ## create boxplot
  dfbox <- ggplot() + 
    geom_boxplot(data = allDataLand2Cut
                 , aes(x=as.factor(data), y=evalCol)
                 , position=position_dodge(1)) +
    ggtitle(nowName) + theme_bw()
  print(dfbox)
  Sys.sleep(1)
}
dev.off()


ggplot(data = allDataLong %>% filter(name == "total_Tco2e")
       , aes(x = data, y = summed_Tco2e)) +
  geom_bar(stat="identity") +
  theme_bw()

ad2 <- allDataLong
ad2$data <- as.factor(ad2$data)
ad2$data <- factor(ad2$data, levels = c("Ag_30", "Ag_15", "baseline2007", "Gr_15", "Gr_30", "d2015"))
ggplot(data = ad2 %>% filter(name == "total_Tco2e" & data != "d2015")
       , aes(x = data, y = summed_Tco2e)) +
  geom_bar(stat="identity") +
  theme_bw()

#### Fig 2 - proportion of emissions (based on origin) assigned to each scenario ####
# keep only required rows
allDataLong2 <- allDataLong %>%
  filter(name %in% c("fertliser_CO2e_kgkm2", "fuel_kgco2e", "kgCO2e_entferm", "kgCO2e_excrete"
                     , "kgCO2e_manMgmt", "landCoverTotalEmis_kgco2"))
allDataLong2$data <- as.factor(allDataLong2$data)
allDataLong2$data <- factor(allDataLong2$data, levels = c("Ag_30", "Ag_15", "baseline2007", "Gr_15", "Gr_30", "d2015"))

ggplot(allDataLong2, aes(x = data, y = summed_Tco2e, fill = name)) +
  geom_col(colour = "black", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Pastel1")  +
  theme_bw()

ggplot(allDataLong2, aes(x = data, y = summed_Tco2e, fill = name)) +
  geom_col(colour = "black", position = "stack") +
  # scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Pastel1")  +
  theme_bw()
