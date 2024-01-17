Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2023-27-01
Last update: 2024-01-17 

The files in this directory contain the data used to create the outputs of the GHG tool. There are three directories:
  1. 'land_cover' = data input relating to land covers, including hedgerows and field margins
  2. 'crops' = data relating to crops
  3. 'animals' = data relating to animals

The other files in this directory are 'agland_grid_ESW.gpkg', and 'ghgGridGB.gpkg' which contains that initial 1 km2 grid that will be used in all future analysis for this GHGtool
The only column is 'rcFid_1km', which is the centroid point of 1 km reference id (EPSG: 27700) [x and y coordinates are seperated by '_']. This column and the values in it will be used to enusre all data is in the same km
Data:
    Spatial resolution: 1 km2
    Spatial extent: GB
    Temporal coverage: N/A
    Native projection: 27700