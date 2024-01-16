Creator: Dr. Paul M. Evans (paueva@ceh.ac.uk) | Git: https://github.com/pevans13
Date created: 2024-01-10
Last update: 2024-01-15 

-----------------------
Overview
-----------------------

This is the readme for the 'GHG_tool_2024' R project (https://github.com/pevans13/GHG_tool_2024.git)

The aim of this project is to combine all of the different elements of greenhouse gas (GHG) emissions calculations to create a final layer showing the GHG emissions or sequestration of 1 km2 cells in Great Britain (GB). 

Important note: this is not an all-encompassing model of GHG emissions, but rather focuses on the GHG activity of agricultural landscapes, adding only supplementary emissions or sequestration values for non-agricultural land uses.  

Specifically, this project aims to model GHG activity from:
	- farmed animals (cows, sheep, pigs, horses, and poultry)
	- farm management (fertiliser use, manure and residue management, and energy use)
	- non-agricultural land uses

-----------------------
Content
-----------------------

File Structure:
    The output is organized into a structured format for ease of use and understanding:
		'r_scripts' directory stores the post-setup scripts
		'data_in' and 'results' store data input, and results, respectively
		'images' stored any figures

Coordinate Reference System (CRS):
	All data were reprojected to EPSG:27700 (OSGB36 / British National Grid; https://epsg.io/27700) before analysis, if required.
	The output of all data is EPSG:27700. 
	
-----------------------
Metadata of input data
-----------------------