# System Filestructure

The system filestructure is broken up into the following directories:

**NOTE**: *Directories which are italicized are not located in the repository (or some larger data files are left out)*

* Archive
* *ATL06*
* *BackgroundData*
* *BAK*
* Bounds
* Documentation
* *Flow*
* ICESat2GroundingLineMigration
* *Line*
* *Saves*
* src

## Archive
Contains old source code and notebooks. Also some scratch work. Not necessary for execution of any programs.

## ATL06
Necessary directory for the `Database` class to work in `ICESat2GroundingLineMigration.ICESatHDF5Unpacker.py` to work properly.
Designed to be used to locally host data in the form of .h5 files for processing.

## BackgroundData
Hosts map background data. 

## BAK
Temporary backups for big changes in code

## Bounds
Holds GIS files for study regions

## Documentation
Documents

## Flow
Directory containing the antarctic ice flow database. Used for computing direction of flow for when calculating flowslope
<br>
Additionally the following 2 data sets are contained in this directory, each of which are used for processing:<br>
* MODIS Mosaic of Antarctica 2008-2009 (MOA2009) Image Map, Version 2<sup>[1][2]</sup>
* MEaSUREs Antarctic Grounding Line from Differential Satellite Radar Interferometry, Version 2<sup>[3]</sup>

## ICESat2GroundingLineMigration
Codebase notebooks and source code in `src` rely on. Used to make working with the ICESat2 data much easier by splitting data into multiple levels including `Database`, `Granule`, `Laser` and `Basemap`.
Additionally contains many utilities for working with the flow data, as well as grounding line data, and visualizing.

## Line
Contains grounding line models

## Saves
Contains geojson files as exported by `S3-InitProcess.ipynb`. These are files which contain all of the ICESat2-ATL06 data with some additional parameters. Mainly flowslope (slope in the direction of ice flow).
Within this folder data is divided into folders by RGT, then by beam. Within each beam folder are the datafiles named as such: `{REGION}-{CYCLE}.json`

## src
Contains source code for later processing of the actual break in slope. Functions in this folder are very heavily used after the geojson files are exported by `S3-InitProcess.ipynb`.

# References

1. Haran, T., J. Bohlander, T. Scambos, T. Painter, and M. Fahnestock. 2014, updated 2019. MODIS <br>
   Mosaic of Antarctica 2008-2009 (MOA2009) Image Map, Version 1. Ice Sheet Grounding Lines. Boulder, <br>
   Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center. <br>
   https://doi.org/10.5067/4ZL43A4619AF. 06-27-2023. <br>

2. Scambos, T., T. Haran, M. Fahnestock, T. Painter, and J. Bohlander. 2007. MODIS-based Mosaic of <br>
   Antarctica (MOA) data sets: Continent-wide surface morphology and snow grain size, Remote <br>
   Sensing of Environment. 111. 242-257. https://doi.org/10.1016/j.rse.2006.12.020 <br>
   
3. Rignot, E., J. Mouginot, and B. Scheuchl. 2016. MEaSUREs Antarctic Grounding Line from <br>
   Differential Satellite Radar Interferometry, Version 2. Ice Sheet Grounding Lines. Boulder, Colorado USA. <br>
   NASA National Snow and Ice Data Center Distributed Active Archive Center. <br>
   https://doi.org/10.5067/IKBWW4RYHF1Q. 02-12-2023. <br>