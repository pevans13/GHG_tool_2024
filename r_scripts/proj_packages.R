#### 0 - load libraries for the GHGtool project ####
## automatic install of packages if they are not installed already
list.of.packages <- c(
'raster'
,'stars'
,'sf'
,'dplyr'
,'data.table'
,'tictoc'
,'terra'
,'tidyr'
,'nngeo'
,'pbapply'
,'stringr'
,'readxl'
,'tidyverse'
,'xml2'
,'rvest'
,'spatialEco'
,'parallel'
,'doParallel'
,'foreach'
,'janitor'
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
