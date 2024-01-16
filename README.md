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

This readme contains the metadata of all the input data, how to get started with the data, the licence information, contact details, and a description of the 'ghgtool_setup.R' script.

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

Land cover maps (LCM)
	LCM 2015 at 1 km2
		Contains the land cover classes as of 2015
		Filename: 'lcm2015_gb_1km_dominant_target_class.img'
		Original data:
			Source: UKCEH (https://doi.org/10.5285/c4035f3d-d93e-4d63-a8f3-b00096f597f5)
			Units: land cover class
			Spatial resolution: 1 km2
			Spatial extent: GB
			Data collected: 2015
			Date produced/ published: 2017
			Native projection: EPSG:27700
			Citations: Rowland, C.S.; Morton, R.D.; Carrasco, L.; McShane, G.; O'Neil, A.W.; Wood, C.M. (2017). Land Cover Map 2015 (1km dominant target class, GB). NERC Environmental Information Data Centre. https://doi.org/10.5285/c4035f3d-d93e-4d63-a8f3-b00096f597f5
	Method:
		- Downloaded data
	
	LCM 2015 at 25 m2
		Contains the land cover classes as of 2015
		Filename: 'LCM2015_GB.tif'
		Original data:
			Source: UKCEH
			Units: land cover class
			Spatial resolution: 25 m2
			Spatial extent: GB
			Data collected: 2015
			Date produced/ published: 2017
			Native projection: EPSG:27700
			Citations: Rowland, C.S.; Morton, R.D.; Carrasco, L.; McShane, G.; O'Neil, A.W.; Wood, C.M. (2017). Land Cover Map 2015 (25m raster, GB). NERC Environmental Information Data Centre. https://doi.org/10.5285/bb15e200-9349-403c-bda9-b430093807c7
	Method:
		- Downloaded data	
		
	LCM (including specific crops) 2015 at 25 m2
		Contains the land cover classes (including specific crops) as of 2015
		Filename: 'LCM2015PlusCrops.tif'
		Original data:
			Source: UKCEH
			Units: land cover class
			Spatial resolution: 25 m2
			Spatial extent: GB
			Data collected: 2015 / 2016
			Date produced/ published: 2017
			Native projection: EPSG:27700
			Citations: 
				Rowland, C.S.; Morton, R.D.; Carrasco, L.; McShane, G.; O'Neil, A.W.; Wood, C.M. (2017). Land Cover Map 2015 (25m raster, GB). NERC Environmental Information Data Centre. https://doi.org/10.5285/bb15e200-9349-403c-bda9-b430093807c7
	Method:
		- Downloaded data	
		- converted arable and improved grassland to specific crops

	LCM 2019 at 25 m2
		Contains the land cover classes as of 2019
		Filename: 'gb2019lcm25m.tif'
		Original data:
			Source: UKCEH
			Units: land cover class
			Spatial resolution: 25 m2
			Spatial extent: GB
			Data collected: 2019
			Date produced/ published: 2020
			Native projection: EPSG:27700
			Citations: Morton, R. D., Marston, C. G., O’Neil, A. W., & Rowland, C. S. (2020). Land Cover Map 2019 (25m rasterised land parcels, GB) [Data set]. NERC Environmental Information Data Centre. https://doi.org/10.5285/F15289DA-6424-4A5E-BD92-48C4D9C830CC
	Method:
		- Downloaded data	
	
Animal data
	AgCensus data (Agricultural Census for England)
	Filenames: 'agcensus_5km.shp'
			   'England_2010_5k.csv'
	Original data:
		Source: UKCEH
		Units: number of animals / crops per cell
		Spatial resolution: 5 km2
		Spatial extent: England
		Data collected: 2010
		Date produced/ published: 2015
		Native projection: EPSG:27700
		Citations: ???

	Beef cattle
	Filename: 'BEEF_2020.csv'
	Original data:
	
	Dairy cattle
	Filename: 'DAIRY_2020.csv'
	Original data:
	
-----------------------
'ghgtool_setup.R' script
-----------------------


-----------------------
Getting Started
-----------------------
Download:
	These files can be downloaded and accessed via https://github.com/pevans13/GHG_tool_2024.git
	
-----------------------
Contact
-----------------------
For questions, concerns, or additional information, please contact paueva@ceh.