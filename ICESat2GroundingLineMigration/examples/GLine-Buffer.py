# Takes in a antartica shapefile and a grounding line shapefile
# Plots the grounding line on the antarctica basemap with the 
# buffer shown.

import GLineUnpacker as gline
import Visualizations as visualize
import geopandas as gpd

path = "PATH TO BASEMAP"

# Import and generate buffers for the grounding line shapefile
glinedata = gline.GLineData("GROUNDING LINE SHAPEFILE PATH")

# Visualize data
visualize.Polygons.line_buffer(gpd.read_file(path), glinedata.linedata, glinedata.linepolygon)
