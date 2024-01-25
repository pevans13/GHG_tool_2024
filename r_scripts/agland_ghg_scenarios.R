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

#### 1 - task1 ####

#### 2 - task2 ####

#### 3 - task3 ####