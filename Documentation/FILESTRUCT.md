# System Filestructure

The system filestructure is broken up into the following directories:

**NOTE**: *Directories which are italicized are not located in the repository*

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
