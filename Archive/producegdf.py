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
import fiona
from math import sqrt

import icepyx as ipx
import s3fs
import earthaccess

# Specifying the necessary icepyx parameters
short_name = 'ATL06'
spatial_extent = 'Bounds/bungen.gpkg' # KML polygon
date_range = ['2021-01-01', '2021-03-31']

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

# grab creds
def grab_creds():
    
    region.earthdata_login(EARTHDATA_USERNAME, EARTHDATA_EMAIL, s3token=True)
    
    #credentials = region._s3login_credentials
    
    # If the above cell results in an error (specifically a KeyError: 'accessKeyId'), please try this approach instead
    credentials = region._session.get("https://data.nsidc.earthdatacloud.nasa.gov/s3credentials").json()
    
    s3 = s3fs.S3FileSystem(key=credentials['accessKeyId'], secret=credentials['secretAccessKey'], token=credentials['sessionToken'])
    
    return s3

s3 = grab_creds()

dataset = unpack.Dataset("ATL06") # ICESat-2 ATL06

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

proj4_crs = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# gen x and y lim from bounds.kml
# Open KML file for use
boundsgdf = gpd.read_file(spatial_extent)
boundsgdf = boundsgdf.to_crs(proj4_crs)
mask = boundsgdf["geometry"].iloc[0]
ylim = (boundsgdf["geometry"].total_bounds[1], boundsgdf["geometry"].total_bounds[3])
xlim = (boundsgdf["geometry"].total_bounds[0], boundsgdf["geometry"].total_bounds[2])
print(f"xlim: {xlim}")
print(f"ylim: {ylim}")

cs = 450 # cell size in m (must be a mulitple of actual grid size (450*n))

# get nearest values to x and y lim which are on grid
fxmin, fxmax = xlim[0] - xlim[0] % cs, xlim[1] - xlim[1] % cs
fymin, fymax = ylim[0] - ylim[0] % cs, ylim[1] - ylim[1] % cs
print(f"Bounds: {fxmin},{fymin},{fxmax},{fymax}")

# compute cells on the y side
cellsy = (fymax-fymin)/cs

# generate grid of flowdatabase
print("Generating grid")
cell = generateGrid(fxmin, fymin, fxmax, fymax, cellsy, flowdatabase.x, flowdatabase.y,
                    flowdatabase.angle, flowdatabase.angle_error)

cell = cell.dropna()

imbuffsc = 0.05

# crop by xlim & ylim
cell = cell.clip(mask)

def calc_direction(lat):
    ascending = [1 if lat[i+1] > lat[i] else -1 for i in range(len(lat)-1)]
    ascending.append(ascending[-1])
    return np.array(ascending)

def account_ascending(ascending, dh_fit_dx):
    dh_fit_dx = dh_fit_dx * ascending
    return dh_fit_dx

from time import time as Time

# Define constants
filecount = 1
filecountmax = 1000
fill = 1
window = 5

fullgdf = gpd.GeoDataFrame(columns=['geometry', 'slope', 'angle', "azumith_in_xy",
                                   "dh_fit_dx", "dh_fit_dy", "date_time", "quality",
                                   "along_track_dist", "name"], geometry='geometry', crs=proj4_crs)

def gen_lasergdf(laser, rgt):
    
    lat, lon, time, dh_fit_dx, dh_fit_dx_sigma, metadata, granuledata = laser.getTrackData()
    datetime = np.array(laser.land_ice_segments["delta_time"][:])
    quality = np.array(laser.land_ice_segments["atl06_quality_summary"][:])
    dh_fit_dy = laser.returnAcrossTrackSlope()

    ascending = calc_direction(lat)
    dh_fit_dx = account_ascending(ascending, dh_fit_dx)
    dh_fit_dy = account_ascending(ascending, dh_fit_dy)
    
    rgt_arr = np.zeros(len(dh_fit_dx))
    rgt_arr = rgt_arr + int(rgt)
    
    name_arr = [metadata["name"]] * len(dh_fit_dx)

    laserdf = pd.DataFrame({'lat': lat, 'lon':lon, "dh_fit_dx":dh_fit_dx, "dh_fit_dy":dh_fit_dy, "dh_fit_dx_sigma":dh_fit_dx_sigma, "date_time":datetime, "quality":quality, "rgt":rgt_arr, "name":name_arr})
    
    laserdf = laserdf[laserdf["dh_fit_dy"] > -1]

    laserdf = laserdf[laserdf.index % fill == 0]

    lasergdf = gpd.GeoDataFrame(laserdf, geometry=gpd.points_from_xy(laserdf.lon, laserdf.lat), crs="EPSG:4326")
    lasergdf = lasergdf.to_crs(proj4_crs)
    
    return lasergdf, ascending

def crop_gdf(lasergdf, xlim, ylim):
    
    trackxy = np.array([(point.x, point.y) for point in lasergdf.geometry])

    lasergdf["x"], lasergdf["y"] = trackxy[:,0], trackxy[:,1]

    # remove points outside of wanted area
    lasergdf = lasergdf[(lasergdf["x"] >= xlim[0]) & (lasergdf["x"] <= xlim[1])]
    lasergdf = lasergdf[(lasergdf["y"] >= ylim[0]) & (lasergdf["y"] <= ylim[1])]
    
    lasergdf["index"] = range(len(lasergdf))
    lasergdf = lasergdf.set_index("index")
    
    # calculate along track distance
    distances = []
    if len(lasergdf["x"]) > 1:
        xdiff = (lasergdf["x"]-lasergdf["x"].shift(periods=1, fill_value=0))**2
        ydiff = (lasergdf["y"]-lasergdf["y"].shift(periods=1, fill_value=0))**2
        xdiff[0], ydiff[0] = 0, 0
        pointdistances = (xdiff + ydiff)**0.5
        pointdistances = pointdistances / 1000 # turn m into km
        for i in range(len(pointdistances)):
            distances.append(pointdistances[:i+1].sum())

        lasergdf["along_track_dist"] = distances
    
    return lasergdf

import scipy.signal as signal

cut_off = 0.6
order = 5
sos = signal.butter(order, cut_off, btype="lowpass", output="sos")

def calc_fs(lasergdf, ascending, cell):
    
    # assign cell gdf to lasergdf
    lasergdf = lasergdf.sjoin_nearest(cell, how="inner")
    
    #if len(lasergdf) >= 20:
        # do butterworth filter
    #    lasergdf["dh_fit_dx"] = signal.sosfilt(sos, np.array(lasergdf["dh_fit_dx"]))

    try:
        # Get the projected azumith at each x and y point
        lasergdf["azumith_in_xy"] = unpack.Basemap.angleTransform(lasergdf["x"][:], lasergdf["y"][:])
    except KeyError:
        # this raises when the length of the values is not greater than 1
        # hence an angle transform cannot be computed
        return "break"

    # Compute slope in the direction of flow
    try:
        lasergdf["flowslopes"] = flowdatabase.get_flow_slopes(lasergdf["dh_fit_dx"], lasergdf["dh_fit_dy"], 
                                                             lasergdf["dh_fit_dx_sigma"], ascending, lasergdf["azumith_in_xy"], 
                                                             lasergdf["angle"], lasergdf["angle_err"])
    except IndexError:
        return "break"

    # Remove None Values
    lasergdf = lasergdf.dropna()
    
    return lasergdf

pd.options.mode.chained_assignment = None

def digitalFilter(lasergdf, window, fullgdf):
    
    if len(lasergdf["flowslopes"]) >= window:
        
        length = len(lasergdf["flowslopes"])
        
        # drop poor quality points
        lasergdf = lasergdf[lasergdf["quality"] == 0]
        
        if len(lasergdf["quality"]) < 0.5 * length:
            return None
        
        # apply a rolling average
        lasergdf["flowslopes"] = lasergdf["flowslopes"].rolling(window).mean()

        # remove edges from trackx and tracky
        lasergdf = lasergdf.dropna()
        
        # make geoseries to contact with full gdf
        tempgdf = gpd.GeoDataFrame({"geometry":lasergdf["geometry"], "slope":lasergdf["flowslopes"], "angle":lasergdf["angle"],
                                   "azumith_in_xy":lasergdf["azumith_in_xy"], "dh_fit_dx":lasergdf["dh_fit_dx"],
                                   "dh_fit_dy":lasergdf["dh_fit_dy"], "date_time":lasergdf["date_time"], "quality":lasergdf["quality"],
                                   "along_track_dist":lasergdf["along_track_dist"], "rgt":lasergdf["rgt"], "name":lasergdf["name"]})
        fullgdf = gpd.GeoDataFrame(pd.concat([fullgdf, tempgdf], ignore_index=True), crs=proj4_crs)
        return fullgdf
    
    else:
        return None
    
def process_file(filename, filecount, filecountmax, gran_ids, starttime, window, fullgdf):

    print(f"Fetching file {filecount}/{len(gran_ids[0])} | Total Time: {Time() - starttime} seconds | Avg Time: {(Time() - starttime)/filecount} seconds       ", end="\r")

    # open granule
    granule = dataset.opens3(filename, s3)
    rgt = granule.start_rgt

    # iterate through lasers in granule
    for laser in granule.lasers:

        print(f"Interpeting file {filecount}/{len(gran_ids[0])} - Laser: {laser.name}| Total Time: {Time() - starttime} seconds | Avg Time: {(Time() - starttime)/filecount} seconds       ", end="\r")

        lasergdf, ascending = gen_lasergdf(laser, rgt)

        lasergdf = crop_gdf(lasergdf, xlim, ylim)

        print(f"Processing file {filecount}/{len(gran_ids[0])} - Laser: {laser.name}| Total Time: {Time() - starttime} seconds | Avg Time: {(Time() - starttime)/filecount} seconds       ", end="\r")

        lasergdf = calc_fs(lasergdf, ascending, cell)


        if type(lasergdf) == str:
            break

        out = digitalFilter(lasergdf, window, fullgdf)
        if type(out) == gpd.geodataframe.GeoDataFrame:
            fullgdf = out
            
    return fullgdf

starttime = Time()

if __name__ == "__main__":

    # iterate through each file
    for filename in gran_ids[1][filecount-1:filecountmax+1]:

        fullgdf = process_file(filename, filecount, filecountmax, gran_ids, starttime, window, fullgdf)

        filecount += 1

        if filecount % 35 == 0 and filecount >= 5:
            print("\nCredentials are about to time out, need to refresh")
            s3 = grab_creds()

        if filecount > filecountmax:
                break
                
fullgdf = fullgdf.clip(mask)

fullgdf.to_file(f"Saves/RAW-{Time()}.json", driver="GeoJSON")