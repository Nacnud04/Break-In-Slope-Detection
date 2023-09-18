import geopandas as gpd
import shapely as sp
import numpy as np
import vaex as vx
import argparse
import netCDF4 as nc

from shapely.geometry import Polygon, mapping, Point
from shapely.prepared import prep
from rasterio.features import geometry_mask
from time import time as Time

import matplotlib.pyplot as plt
import icepyx as ipx
import s3fs
import pyproj
import multiprocessing

from rich import print, pretty
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)


pretty.install()

from datetime import datetime
def dtm():
    return f'[{datetime.now().strftime("%H:%M:%S")}]'

def core_name():
    return int(multiprocessing.current_process().name.split('-')[-1]) % 4

def crs():
    return "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

@vx.register_dataframe_accessor('mytool', override=True)
class mytool:
    def __init__(self, df):
        self.df = df

    def shift(self, column, n, inplace=False):
        # make a copy without column
        df = self.df.copy().drop(column)
        # make a copy with just the colum
        df_column = self.df[[column]]
        # slice off the head and tail
        df_head = df_column[-n:]
        df_tail = df_column[:-n]
        # stitch them together
        df_shifted = df_head.concat(df_tail)
        # and join (based on row number)
        return df.join(df_shifted, inplace=inplace)

def spatial_extent(filename):
    bgdf = gpd.read_file(filename)
    bgdf = bgdf.to_crs(crs())
    mask = bgdf['geometry'].iloc[0]
    lims = bgdf['geometry'].total_bounds[:]
    return mask, (lims[0], lims[2]), (lims[1], lims[3])

def load_flowslope(filepath, mask, xlim, ylim, gpu):

    print(f"{dtm()} - Importing flow vectors... ")
    f = nc.Dataset(filepath)

    x = f.variables['x'][:]
    y = f.variables['y'][:]
    
    xmin, xmax = xlim
    ymin, ymax = ylim

    print(f"{dtm()} - Cropping geometries...")
    # All indices in bounding box:
    where_j = np.where((y >= ymin) & (y <= ymax))[0]
    where_i = np.where((x >= xmin) & (x <= xmax))[0]

    # Start and end+1 indices in each dimension:
    i0 = where_i[0]
    i1 = where_i[-1]+1

    j0 = where_j[0]
    j1 = where_j[-1]+1
    
    print(f"{dtm()} - Computing flow vectors...")
    if not gpu:
        
        VX = f.variables['VX'][j0:j1, i0:i1]
        VY = f.variables['VY'][j0:j1, i0:i1]
        angle = np.arctan(VY / VX)
        
        xarr = f.variables['x'][i0:i1]
        yarr = f.variables['y'][j0:j1]
        xlen, ylen = len(xarr), len(yarr)
        
        mshx, mshy = np.meshgrid(xarr, yarr)
        
        polygons = []
        
        # Assuming ylen and xlen are defined somewhere in your code
        with Progress() as progress:
            task = progress.add_task(f"{dtm()} - [cyan]Mapping flow geometries...", total=(ylen - 1) * (xlen - 1))

            for i in range(ylen - 1):
                for j in range(xlen - 1):
                    polygon = Polygon([(mshx[i, j], mshy[i, j]),
                                      (mshx[i + 1, j], mshy[i + 1, j]),
                                      (mshx[i + 1, j + 1], mshy[i + 1, j + 1]),
                                      (mshx[i, j + 1], mshy[i, j + 1])])
                    polygons.append(polygon)
                    progress.update(task, advance=1)
                
        angle = angle[:-1, :-1]
                
        print(f"{dtm()} - Migrating into geodataframe...")
        gdf = gpd.GeoDataFrame({'geometry': polygons, 'angle': angle.flatten()})
        
        print(f"{dtm()} - Masking geodataframe...")
        gdf = gdf.clip(mask)
        
    else:
        
        print(gpu)
        
    return gdf


def file_query(spt_ext, cyc):
    
    st = Time()
    # set up query object
    region = ipx.Query("ATL06", spt_ext, cycles = cyc)
    # grab s3 links
    gran_ids = region.avail_granules(ids=True, cloud=True)
    links = gran_ids[1]
    print(f"{dtm()} - Found {len(links)} granules in [bright_cyan]{round(Time() - st, 5)}[/bright_cyan]s")
    
    return region, links
    

def gen_s3(region):
    region.earthdata_login(s3token=True)
    credentials = region._session.get("https://data.nsidc.earthdatacloud.nasa.gov/s3credentials").json()
    s3 = s3fs.S3FileSystem(key=credentials['accessKeyId'], secret=credentials['secretAccessKey'], token=credentials['sessionToken'])
    return s3


def find_nearest(arr, val):
    return int((np.abs(arr-val)).argmin())


def make_laser(laser, rgt, region, cycle, name, bounds, GPU):
    
    land_ice_segments = laser["land_ice_segments"]
    
    lat, lon = land_ice_segments["latitude"], land_ice_segments["longitude"]
    
    xlim, ylim = bounds
    
    if not GPU:
        
        # transform crs
        source_proj4 = '+proj=latlong +datum=WGS84'
        target_proj4 = crs()
        transformer = pyproj.Transformer.from_proj(pyproj.Proj(source_proj4), pyproj.Proj(target_proj4), always_xy=True)
        x, y = transformer.transform(lon[:], lat[:])
        # find where track escapes bounding box
        idxmin, idxmax = find_nearest(x, xlim[0]), find_nearest(x, xlim[1])
        idymin, idymax = find_nearest(y, ylim[0]), find_nearest(y, ylim[1])
    
        
    idmin = min(idxmin, idxmax, idymin, idymax)
    idmax = max(idxmin, idxmax, idymin, idymax)

    fit_statistics = land_ice_segments["fit_statistics"]
    
    data = {"x":x[idmin:idmax], "y":y[idmin:idmax], "time":land_ice_segments["delta_time"][idmin:idmax], 
            "h_li":land_ice_segments["h_li"][idmin:idmax], "dh_fit_dx":fit_statistics["dh_fit_dx"][idmin:idmax], 
            "dh_fit_dy":fit_statistics["dh_fit_dy"][idmin:idmax], 
            "azumith":land_ice_segments["ground_track"]["ref_azimuth"][idmin:idmax],
            "quality":land_ice_segments["atl06_quality_summary"][idmin:idmax]}
    
    df = vx.from_dict(data)
    
    return df


def reduce_dfs(lasers, mask, clpbMsk=False):
    
    nlasers = []
    for df in lasers:
        if clpbMsk:
            keep = []
            tot = len(df)
            for i, row in df.iterrows():
                keep.append(mask.contains(Point(row['x'], row['y'])))
                print(f"{round(100*i/tot, 2)}", end="\r")
            df = df[keep]
        df = df[df["quality"] == 0] # remove points of low quality
        # get rid of extreme or obviously incorrect points
        df = df[(df["dh_fit_dy"] > -1) & (df["dh_fit_dy"] < 1)]
        df = df.extract() # extract and filter
        nlasers.append(df)
    return nlasers
    


def along_track_distance(lasers):
    
    out = []
    for df in lasers:
        #print(f"x:{len(df['x'][:])} | y:{len(df['y'][:])}")
        #df['x_shifted'] = df["x"]
        #df = df.mytool.shift('x_shifted', 1)
        #df['y_shifted'] = df["y"]
        #df = df.mytool.shift('y_shifted', 1)
        
        x, y = df["x"].values, df["y"].values
        
        x_s = np.roll(x, 1)
        y_s = np.roll(y, 1)

        # Calculate squared differences
        xd2 = (x - x_s) ** 2
        yd2 = (y - y_s) ** 2
        
        # Calculate the total distance as the square root of the sum of squared differences
        df['distance'] = (xd2 + yd2) ** 0.5
        
        print(f"{x}\n{y}\n{df['distance'].values}")
        out.append(df)
    
    return out