import ICESat2GroundingLineMigration.IceSatHDF5Unpacker as unpack
import ICESat2GroundingLineMigration.GLineUnpacker as gline
import ICESat2GroundingLineMigration.FlowUnpacker as flow
import ICESat2GroundingLineMigration.Visualizations as visualize
import numpy as np
import geopandas as gpd
import pandas as pd
import shapely
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

import icepyx as ipx
import fiona

# Specifying the necessary icepyx parameters
short_name = 'ATL06'
spatial_extent = 'doc.kml' # KML polygon centered on Sermeq Kujalleq
date_range = ['2019-05-01', '2019-06-30']

# Open KML file for use
fiona.drvsupport.supported_drivers['LIBKML'] = 'rw' # enable KML support which is disabled by default
jk = gpd.read_file(spatial_extent)

# Setup the Query object
region = ipx.Query(short_name, spatial_extent, date_range)

plotspat_ext = list(region.spatial.extent_as_gdf.geometry.unary_union.minimum_rotated_rectangle.bounds)
plotreg = ipx.Query(short_name, plotspat_ext, date_range)

# Grabbing granule s3 urls
gran_ids = region.avail_granules(ids=True, cloud=True)
print(f"Found {len(gran_ids[0])} granules")

# Establish credentials
EARTHDATA_USERNAME = "byrne"
EARTHDATA_EMAIL = "byrne@mines.edu"

region.earthdata_login(EARTHDATA_USERNAME, EARTHDATA_EMAIL, s3token=True)

credentials = region._s3login_credentials

import s3fs

s3 = s3fs.S3FileSystem(key=credentials['accessKeyId'],
                       secret=credentials['secretAccessKey'],
                       token=credentials['sessionToken'])

proj4_crs = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
dataset = unpack.Dataset("ATL06") # ICESat-2 ATL06

# grounding line data with a buffer of 50000m
groundline = gline.GLineData("Line/InSAR_GL_Antarctica_v02.shp", 50000)

# import flow database
flowdatabase = flow.Database("Flow/antarctic_ice_vel_phase_map_v01.nc")
flowdatabase.compute_all()

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return int(idx)

def generateGrid(xmin, ymin, xmax, ymax, ncellsy, flowx, flowy, flowangles, flowerr):
    
    cell_size = (ymax-ymin)/ncellsy
    print(f"cell_size = {cell_size}")
    
    # create the cells in a loop
    grid_cells = []
    angles = []
    err = []
    for x0 in np.arange(xmin, xmax+cell_size, cell_size):
        for y0 in np.arange(ymin, ymax+cell_size, cell_size):
            # bounds
            x1 = x0-cell_size
            y1 = y0+cell_size
            ix = find_nearest(flowx, x0)
            iy = find_nearest(flowy, y0)
            grid_cells.append(shapely.geometry.box(x0, y0, x1, y1))
            angles.append(flowangles[int(iy),int(ix)])
            err.append(flowerr[int(iy), int(ix)])
    cells = gpd.GeoDataFrame({"geometry":grid_cells, "angle":angles, "angle_err":err}, geometry = "geometry", crs=proj4_crs)
    
    return cells

xlim, ylim = (-400000,200000), (-750000,-300000)
#xlim, ylim = (-3e6,3e6), (-3e6,3e6)
cs = 450 # cell size in m (must be a mulitple of actual grid size (450*n))

# create a mask
mask = shapely.geometry.box(xlim[0], ylim[0], xlim[1], ylim[1])

# get nearest values to x and y lim which are on grid
fxmin, fxmax = xlim[0] - xlim[0] % cs, xlim[1] - xlim[1] % cs
fymin, fymax = ylim[0] - ylim[0] % cs, ylim[1] - ylim[1] % cs
print(f"Bounds: {fxmin},{fymin},{fxmax},{fymax}")

# compute cells on the y side
cellsy = (fymax-fymin)/cs

#flowdatabase.x = flowdatabase.x[::fill]
#flowdatabase.y = flowdatabase.y[::fill]
#flowdatabase.angle = flowdatabase.angle[::fill, ::fill]
#flowdatabase.angle_error = flowdatabase.angle_error[::fill, ::fill]

# generate grid of flowdatabase
cell = generateGrid(fxmin, fymin, fxmax, fymax, cellsy, 
                    flowdatabase.x, flowdatabase.y,
                    flowdatabase.angle, flowdatabase.angle_error)

cell = cell.dropna()

def digitalFilter(lasergdf, fill, window, fullgdf):
    if len(lasergdf["flowslopes"]) >= fill*window:

        # only use every fill(th) point
        lasergdf.iloc[::fill]

        # apply a rolling average
        lasergdf["flowslopes"] = lasergdf["flowslopes"].rolling(window).mean()

        # remove edges from trackx and tracky
        lasergdf.dropna()

        # remove points outside of 2*std deviation
        stddev = lasergdf["flowslopes"].std()
        average = lasergdf["flowslopes"].mean()
        lasergdf[lasergdf["flowslopes"] <= (average+2*stddev)]
        lasergdf[lasergdf["flowslopes"] >= (average-2*stddev)]
        
        # make geoseries to contact with full gdf
        tempgdf = gpd.GeoDataFrame({"geometry":lasergdf["geometry"], "slope":lasergdf["flowslopes"]})
        fullgdf = gpd.GeoDataFrame(pd.concat([fullgdf, tempgdf], ignore_index=True), crs=proj4_crs)
        return fullgdf
    else:
        return None
    
from time import time as Time

# Define constants
filecount = 75
filecountmax = 150
fill = 10
window = 10
folder = "ATL06"

fullgdf = gpd.GeoDataFrame(columns=['geometry', 'slope'], geometry='geometry', crs=proj4_crs)

flowdatabase.compute_all() # computes all parameters from the flowdatabase

starttime = Time()

# iterate through each file
for filename in gran_ids[1][filecount:filecountmax+1]:
    
    print(f"Downloading file {filecount}/{filecountmax} | Total Time: {Time() - starttime} seconds | Avg Time: {(Time() - starttime)/filecount} seconds       ", end="\r")
    
    # open granule (this opens from the hard drive. Not from the cloud)
    granule = dataset.opens3(filename, s3)
    
    # iterate through lasers in granule
    for laser in granule.lasers:
        
        print(f"Processing file {filecount}/{filecountmax} - Laser: {laser.name}| Total Time: {Time() - starttime} seconds | Avg Time: {(Time() - starttime)/filecount} seconds       ", end="\r")
        
        lat, lon, time, dh_fit_dx, dh_fit_dx_sigma, metadata, granuledata = laser.getTrackData()
        dh_fit_dy = laser.returnAcrossTrackSlope()
        
        laserdf = pd.DataFrame({'lat': lat, 'lon':lon, "dh_fit_dx":dh_fit_dx, "dh_fit_dy":dh_fit_dy, "dh_fit_dx_sigma":dh_fit_dx_sigma})
        
        # only keep every 75th point
        laserdf = laserdf[laserdf.index % 200 == 0]
        
        lasergdf = gpd.GeoDataFrame(laserdf, geometry=gpd.points_from_xy(laserdf.lon, laserdf.lat), crs="EPSG:4326")
        lasergdf = lasergdf.to_crs(proj4_crs)
        trackxy = np.array([(point.x, point.y) for point in lasergdf.geometry])
        
        # remove points outside of wanted area
        lasergdf.clip(mask)
        
        # assign cell gdf to lasergdf
        lasergdf = lasergdf.sjoin_nearest(cell, how="inner")

        # Get the projected azumith at each x and y point
        azumith_in_xy = unpack.Basemap.angleTransform(trackxy[:,0], trackxy[:,1])
        # Compute slope in the direction of flow
        try:
            lasergdf["flowslopes"] = flowdatabase.get_flow_slopes(lasergdf["dh_fit_dx"], lasergdf["dh_fit_dy"], 
                                                                 lasergdf["dh_fit_dx_sigma"], azumith_in_xy, 
                                                                 lasergdf["angle"], lasergdf["angle_err"])
        except IndexError:
            break

        # Remove None Values
        lasergdf = lasergdf.dropna()
        
        out = digitalFilter(lasergdf, fill, window, fullgdf)
        if type(out) == gpd.geodataframe.GeoDataFrame:
            fullgdf = out
        
    filecount += 1
    if filecount > filecountmax:
        break

cell = cell.sjoin_nearest(fullgdf, how="inner", max_distance = 4500)
print(cell.keys())

cell.to_file(f"Saves/{filecount-filecountmax}-{filecount}.json", driver="GeoJSON")