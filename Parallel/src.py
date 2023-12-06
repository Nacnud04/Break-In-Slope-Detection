import geopandas as gpd
import pandas as pd
import shapely as sp
import numpy as np
import vaex as vx
import argparse
import netCDF4 as nc
import math

from cartopy import crs as ccrs
from shapely.geometry import Polygon, mapping, Point, LineString, MultiPoint
from shapely.ops import split, unary_union
from shapely.prepared import prep
from rasterio.features import geometry_mask
from time import time as Time

import matplotlib.pyplot as plt
import icepyx as ipx
import s3fs
import pyproj
import multiprocessing
import earthaccess
import scipy.signal as signal
from scipy.spatial import cKDTree

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

pd.options.mode.chained_assignment = None

from datetime import datetime
def dtm():
    return f'[{datetime.now().strftime("%H:%M:%S")}]'

cores = multiprocessing.cpu_count()
def core_name():
    return int(multiprocessing.current_process().name.split('-')[-1]) % cores

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
        gdf = gpd.GeoDataFrame({'angle': angle.flatten()}, geometry=polygons, crs=crs())
        
        print(f"{dtm()} - Masking geodataframe...")
        gdf = gdf.clip(mask)
        
    else:
        
        print(gpu)
        
    return gdf


def load_gline(path, xlim, ylim):
    xs, ys = [], []
    data = gpd.read_file(path)
    data = data.to_crs(crs())
    line_strings = [LineString(p.exterior.coords) for p in data.iloc[0].geometry.geoms]
    data = gpd.GeoDataFrame({'geometry': line_strings})
    data = data.clip_by_rect(xlim[0], ylim[0], xlim[1], ylim[1])
    return data


def file_query(spt_ext, cyc):
    
    st = Time()
    # set up query object
    region = ipx.Query("ATL06", spt_ext, cycles = cyc)
    # grab s3 links
    gran_ids = region.avail_granules(ids=True, cloud=True)
    links = gran_ids[1]
    print(f"{dtm()} - Found {len(links)} granules in [bright_cyan]{round(Time() - st, 5)}[/bright_cyan]s")
    
    return region, links
    

def gen_s3():
    auth = earthaccess.login()
    credentials = auth.get_s3_credentials(daac="NSIDC")
    credentials['aws_access_key_id'] = credentials['accessKeyId']
    credentials['aws_secret_access_key'] = credentials['secretAccessKey']
    credentials['aws_session_token'] = credentials['sessionToken']
    return credentials


def find_nearest(arr, val):
    return int((np.abs(arr-val)).argmin())


def calc_direction(lat):
    ascending = [1 if lat[i+1] > lat[i] else -1 for i in range(len(lat)-1)]
    ascending.append(ascending[-1])
    return np.array(ascending)


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
    
        
    idxmin, idxmax = min(idxmin, idxmax), max(idxmin, idxmax)
    idymin, idymax = min(idymin, idymax), max(idymin, idymax)
    idmin = max(idxmin, idymin)
    idmax = min(idxmax, idymax)

    fit_statistics = land_ice_segments["fit_statistics"]
    dh_fit_dx = fit_statistics["dh_fit_dx"]
    
    # account for direction of satellite. (ascending/descending)
    direc = calc_direction(lat[:])
    dh_fit_dx = dh_fit_dx * direc
    
    data = {"x":x[idmin:idmax], "y":y[idmin:idmax], "delta_time":land_ice_segments["delta_time"][idmin:idmax], 
            "h_li":land_ice_segments["h_li"][idmin:idmax], "dh_fit_dx":dh_fit_dx[idmin:idmax], 
            "dh_fit_dy":fit_statistics["dh_fit_dy"][idmin:idmax], 
            "azumith":land_ice_segments["ground_track"]["ref_azimuth"][idmin:idmax],
            "ascending":direc[idmin:idmax],
            "quality":land_ice_segments["atl06_quality_summary"][idmin:idmax]}
    
    df = vx.from_dict(data)
    
    return df


def reduce_dfs(lasers, beams, mask):
    minLen = 10
    nlasers = []
    bs = []
    for df, b in zip(lasers, beams):
        df = df[df["quality"] == 0] # remove points of low quality
        # get rid of extreme or obviously incorrect points
        df = df[(df["dh_fit_dy"] > -1) & (df["dh_fit_dy"] < 1)]
        df = df.dropna()
        if len(df) > minLen:
            nlasers.append(df)
            bs.append(b)
    return nlasers, bs
    
def XYvec_to_ang(origin, vector):
    """
    Takes an xy vector and returns the amount of degrees from the north axis it is.

    Parameters
    ----------
    origin : tuple
          (x, y) - Coordinate at which angle is from
    vector : tuple
          (x, y) - Point which forms vector

    Returns
    -------
    angle : float
          Angle from x in radians on the xy plane.
    """
    dx = vector[0] - origin[0]
    if dx == 0:
        return math.pi / 2
    dy = vector[1] - origin[1]
    angle = math.atan(dy/dx)
    return angle


def angleTransform(xvals, yvals):
    """Performs XYvec_to_ang on an array of data"""
    angles = []
    for i in range(len(xvals)):
        if i >= 0 and i < len(xvals) - 1:
            angles.append(XYvec_to_ang((xvals[i], yvals[i]),(xvals[i+1], yvals[i+1])))
        elif i == len(xvals) - 1:
            angles.append(XYvec_to_ang((xvals[i-1], yvals[i-1]),(xvals[i], yvals[i])))
    angles = np.array(angles)
    return angles

def dist_azumith(lasers, beams, GPU=False):
    
    out = []
    bs = []
    for df, b in zip(lasers, beams):
        if not GPU:

            x, y = df["x"].values, df["y"].values

            x_s = np.roll(x, 1)
            y_s = np.roll(y, 1)

            # Calculate differences
            dx = (x - x_s)
            dy = (y - y_s)

            # Calculate the total distance as the square root of the sum of squared differences
            distances = ((dx**2 + dy**2) ** 0.5)/1000
            distances[0] = 0
            df["along_track_distance"] = np.cumsum(distances)
            
            # projected azumith from x and y coordinates
            frac = [y/x if x != 0 else float('inf') for x, y in zip(dx, dy)]
            df["azumith"] = angleTransform(x, y) #np.arctan(frac)

            out.append(df)
            bs.append(b)

    return out, bs

def calc_along_flow_slope(dh_fit_dx, dh_fit_dy, ascending, slope_azumith, flow_ang):
    """
    Takes the slope and across track slope and compute the slope in the direction of flow

    Parameters
    ----------
    dh_fit_dx : int/float
          Slope in the direction of laser track
    dh_fit_dy : int/float
          Across track slope orthogonal to the direction of the laser track
    dh_fit_dx_sigma : int/float
          Slope error in the direction of the laser track
    slope_azumith : int/float
          Azumith data for the laser track
    flow_ang : int/float
          Direction of flow in the xy plane
    flow_ang_error : int/float
          Error in the direction of flow
    sc_orient : int
          Int which says the direction of orbit.

    Returns
    -------
    None : None
          Data at the point did not have any across track slope data
    slopeflowgrade : flow
          Predicted slope in direction of flow
    """
    slopevec = np.array([math.cos(slope_azumith), math.sin(slope_azumith), dh_fit_dx])
    slopevec = slopevec / (slopevec **2).sum()**0.5
    if ascending == 1:# and slope_azumith < 1:
        acrossslopevec = np.array([math.sin(slope_azumith), -1*math.cos(slope_azumith), dh_fit_dy])
        acrossslopevec = acrossslopevec / (acrossslopevec**2).sum()**0.5
    elif ascending == -1:# or slope_azumith > 1:
        #slope_azumith = slope_azumith * -1
        acrossslopevec = np.array([-1*math.sin(slope_azumith), math.cos(slope_azumith), dh_fit_dy])
        acrossslopevec = acrossslopevec / (acrossslopevec**2).sum()**0.5
    else:
        return None
    flowvec = np.array([math.cos(flow_ang), math.sin(flow_ang), 0])
    slopeflowvec = (flowvec.dot(slopevec)/slopevec.dot(slopevec))*slopevec+(flowvec.dot(acrossslopevec)/acrossslopevec.dot(acrossslopevec))*acrossslopevec
    slopeflowgrade = slopeflowvec[2] / math.sqrt(slopeflowvec[0]**2 + slopeflowvec[1]**2)
    if dh_fit_dy >= 1e30:
        return None
    else:
        return slopeflowgrade

def get_flow_slopes(dh_fit_dx, dh_fit_dy, ascending, slope_azumith, flow_ang):
    """Performs calc_along_flow_slope on an array of data"""
    dh_fit_dx, dh_fit_dy, flow_ang = np.array(dh_fit_dx), np.array(dh_fit_dy), np.array(flow_ang)
    slope_azumith  = np.array(slope_azumith)
    along_flow_slope = np.array([calc_along_flow_slope(dh_fit_dx[i], dh_fit_dy[i], ascending[i], slope_azumith[i], flow_ang[i]) for i in range(len(flow_ang))])
    along_flow_slope = along_flow_slope
    return along_flow_slope


def comp_flowslope(lasers, beams, flowgdf, GPU=False):
    out = []
    bs = []
    for df, b in zip(lasers, beams):
        if not GPU:
            
            # assign cell gdf to lasergdf
            x, y = df['x'].values, df['y'].values
            pddf = pd.DataFrame({'x':x,'y':y})
            lasergdf = gpd.GeoDataFrame(pddf, geometry=[Point(x, y) for x, y in zip(x,y)], crs=crs())
            lasergdf = lasergdf.sjoin_nearest(flowgdf, how="inner")
            if len(df) + 1 == len(lasergdf):
                lasergdf.drop(lasergdf.tail(1).index,inplace=True)
            
            slope_azumiths = df['azumith'].values
            ascending = df['ascending'].values
            dh_fit_dxs = df['dh_fit_dx'].values
            
            # gen along track vector, z component is magnitude of dh_fit_dx
            slopevecs = [np.array([np.cos(slope_azumith), np.sin(slope_azumith), dh_fit_dx]) for slope_azumith, dh_fit_dx in zip(slope_azumiths, dh_fit_dxs)]
            slopevecs = [slopevec / (slopevec**2).sum()**0.5 for slopevec in slopevecs]
            
            # gen across track vector, z component is magnitude of dh_fit_dy
            across_slope_vecs = [np.array([math.sin(sa), -1*math.cos(sa), dhdy]) if asc==1 else np.array([-1*math.sin(sa), math.cos(sa), dhdy]) for sa, asc, dhdy in zip(slope_azumiths, ascending, df['dh_fit_dy'].values)]
            across_slope_vecs = [across_slope_vec / (across_slope_vec**2).sum()**0.5 for across_slope_vec in across_slope_vecs]
            
            # gen ice flow vector, and make z component 0
            angles = lasergdf["angle"].values
            flowvecs = [np.array([np.cos(angle), np.sin(angle), 0]) for angle in angles]
                                
            # project the flow vector onto a plane crated from the along and across track slope vectors
            slopeflowvecs = [(flowvec.dot(slopevec) / slopevec.dot(slopevec)) * slopevec + (flowvec.dot(acrossslopevec) / acrossslopevec.dot(acrossslopevec)) * acrossslopevec for slopevec, acrossslopevec, flowvec in zip(slopevecs, across_slope_vecs, flowvecs)]
                                
            # transform the flow vector into a slope
            slopeflowgrade = [slopeflowvec[2] / math.sqrt(slopeflowvec[0]**2 + slopeflowvec[1]**2) for slopeflowvec in slopeflowvecs]
            
            try:
                df['slope'] = get_flow_slopes(dh_fit_dxs, df['dh_fit_dy'], ascending, slope_azumiths, angles) #np.array(slopeflowgrade)
            except IndexError:
                print(f"{dtm()} - CORE: {core_name()} - [bold red]ERROR:[/bold red] [red]Flow slope computation failed[/red]")
            
            # apply a rolling average
            df["slope"] = df["slope"].rolling(5).mean()
            df["angle"] = angles
            
            out.append(df)
            bs.append(b)
            
    return out, bs
            
        
def get_intersections(lines):
    point_intersections = []
    line_intersections = [] #if the lines are equal the intersections is the complete line!
    lines_len = len(lines)
    for i in range(lines_len):
        for j in range(i+1, lines_len): #to avoid computing twice the same intersection we do some index handling
            l1, l2 = lines[i], lines[j]
            if l1.intersects(l2):
                intersection = l1.intersection(l2)
                if isinstance(intersection, LineString):
                    line_intersections.append(intersection)
                elif isinstance(intersection, Point):
                    point_intersections.append(intersection)
                elif isinstance(intersection, MultiPoint):
                    point_intersections.append(list(intersection.geoms))
                else:
                    print(intersection)
                    raise Exception('What happened?')

    return point_intersections, line_intersections
        
        
def get_intersect_dist(df, intersects, index):
    xi, yi = intersects[index].xy    
    # Calculate the distance from the intersection to each point in the list
    row = df.iloc[0] # capture row 1 
    x, y = row['x'], row['y']
    distance = ((x-xi)**2+(y-yi)**2)**0.5
    distance = distance / 1000
    return distance
        
        
def offset_by_intersect(df, intersects):
    
    intersects = np.flip(intersects) 
    for i in range(len(intersects)):
        # the following line is not necessary but it will produce results consistent
        # with that of the non-parallel process
        distance = get_intersect_dist(df, intersects, i)
        if distance >= 15 and distance <= df['along_track_distance'].iloc[-1] - 20:
            break
    # in theory distance should be approximately the distance of the along track distance at the
    # point of intersection. therefore distance can be subtracted from the along track dist
    #print(f"{dtm()} - CORE: {core_name()} - Along track offset: {distance}")
    df['along_track_distance'] = df['along_track_distance'] - distance
    return df


def interpolation(idx, track):
    
    """
    Interpolate the dataframe over indicies based on idx
    
    Paramters
    ---------
    np.ndarray: idx
          array of along_dist to interpolate points at (np.linspace is good for this)
    pd.DataFrame: track
          dataframe of track data
          output from process_and_extract function
          
    Returns
    -------
    pd.DataFrame: df
          interpolated data as a dataframe
    """
    
    along_dist = track["along_track_distance"].values
    slope_raw = np.interp(idx, along_dist, track["slope"].values)
    slope_filt = np.interp(idx, along_dist, track["slope-filt"].values)
    h_li = np.interp(idx, along_dist, track["h_li"].values)
    x = np.interp(idx, along_dist, track['x'].values)
    y = np.interp(idx, along_dist, track['y'].values)
    dtm = np.interp(idx, along_dist, track["delta_time"].values)
    df = pd.DataFrame(index = idx, data = {"along_track_distance":idx, "delta_time":dtm, "slope":slope_raw, "slope-filt":slope_filt, "h_li":h_li, 'x':x, 'y':y})
    
    return df


def clean_by_gaps(track, df, max_dist, count):
    
    """
    Removes large gaps of missing data which were interpolated over
    
    Parameters
    ----------
    pd.DataFrame: track
          dataframe which contains all the non-interpolated track data
    pd.DataFrame: df
          dataframe which contains the interpolated tracl data
    int: max_dist
          maximum distance between two points which isnt considered a gap
          (any larger is considered a gap)
    int: count
          count of all previous gaps removed.
          
    Returns
    -------
    pd.DataFrame: df
          updated output dataframe
    int: count
          updated gap removal count
    """
    
    along_dist = np.array(track["along_track_distance"].values)
    delta_dist = np.diff(along_dist)
    for j, delt in enumerate(delta_dist):
        if delt > (max_dist / 1000) and j < len(delta_dist):
            min_idx, max_idx = along_dist[j], along_dist[j + 1]
            df.loc[(df["along_track_distance"] > min_idx) & (df["along_track_distance"] < max_idx), ["slope", "h_li"]] = None
            #print(f"Removed data from {min_idx}km to {max_idx}km of delta {delt * 1000}m", end = "           \r")
            count += 1
    return df, count
    

def apply_butter(series, order, cut_off):
    b, a = signal.butter(order, cut_off, btype="lowpass")
    series = signal.filtfilt(b, a, series)
    return series
    

def interp_clean_single(track, fidelity):
    
    """
    Cleans and interpolates a 2d array of track data
    
    Parameters
    ----------
    list: cropped
          2d array where each item is a dataframe for a track
    int: fidelity
          count of how many individual points to create in the interpolation
    int: d_min
          minimum along track distance from the grounding line intersection
    int: d_max
          maximum along track distance from the grounding line intersection
          
    Returns
    -------
    list: interped
          similar to input variable cropped.
          contains the interpolated and cleaned data
    np.ndarray: idx
          values which the data is interpolated along
    """
    
    d_min = track["along_track_distance"].min() + 0.5
    d_max = track["along_track_distance"].max() - 0.5
    track = track[(track["along_track_distance"] > d_min) & (track["along_track_distance"] < d_max)]
    
    idx = np.linspace(d_min, d_max, fidelity)

    dx = (d_max - d_min) / fidelity

    #df_full = interpolation(idx, track)
    # remove points which were interpolated over large chunks of missing data
    order = 5
    cutoff = 0.032 # found by https://doi.org/10.5194/tc-14-3629-2020
    try:
        track["slope-filt"] = apply_butter(track["slope"].values, order, cutoff) # uses scikit to apply the filter.
    except ValueError:
        return None
    
    max_dist = 40 # value in meters
    count = 0
    df = interpolation(idx, track)
    
    df, count = clean_by_gaps(track, df, max_dist, count)
    # remove points below a certain standard deviation threshold
    for dist in reversed(range(int(d_min * 1000), int(d_max * 1000), 100)):
        dist /= 1000
        sel = df[df["along_track_distance"] < dist]
        if sel["slope-filt"].std() < 1e-10:
            df.loc[df["along_track_distance"] < dist, ["slope-filt", "h_li"]] = None
            break
            
    #print(f"Deviation of {rgt}-{name}-{cycle}: {df['slope-filt'].std()}", end = "                                                         \r")
    if df["slope-filt"].values.std() < 1e-10:
        return None
    return df
    
    
def comp_deriv(series, dist):
    derivs = []
    for i in range(len(series) - 1):
        deriv = (series.iloc[i + 1] - series.iloc[i]) / (dist.iloc[i + 1] - dist.iloc[i])
        derivs.append(deriv)
    derivs.append(derivs[-1])
    return np.array(derivs)


def deriv_on_gpd(gpd):
    gpd["slope_deriv_1"] = comp_deriv(gpd["slope-filt"], gpd["along_track_distance"])
    order = 5
    cutoff = 0.3
    gpd[f"slope_deriv_1-filt"] = apply_butter(gpd["slope_deriv_1"], order, cutoff)
    return gpd


def find_peaks(track, amp, dist):
    
    """
    Finds peaks in the slope break
    
    Parameters:
    -----------
    pandas.DataFrame/geopandas.GeoDataFrame: track
          Track to scan for peaks in
    float: amp
          Magnitude requirement to count as a peak
    int: dist
          Minimum distance peaks need to be apart to get detected. If multiple are within this the one with the highest magnitude is chosen.
          
    Returns:
    np.ndarray: peak_dists
          List of along track distances where peaks are.
    """
    
    # find negative peaks and put into array
    peak_index = np.array(signal.find_peaks(track["slope_deriv_1"]*-1, height=amp, distance=dist)[0])
    # append positive peaks to array
    peak_index = np.append(peak_index, np.array(signal.find_peaks(track["slope_deriv_1"], height=amp, distance=dist)[0]))
    # convert indicies into along track distance
    peak_dists = (((track["along_track_distance"].max() - track["along_track_distance"].min()) / len(track)) * peak_index) + track["along_track_distance"].min()
    return peak_dists


def nan_test(local, nan_max):
    
    """
    Tests to see if track region contains too many nan values
    
    Parameters:
    -----------
    pandas.DataFrame/geopandas.GeoDataFrame: local
          Selected region of ground track to check for nans
    float: nan_max
          Maximum percentage of nan values needed to return True (any greater will return false)
    
    Returns:
    --------
    bool: output
          True/False depending on if the amount of nans in the selected region exceeded the nan_max threshold.
    """
    
    nan = np.count_nonzero(np.isnan(local["slope-filt"]))
    if nan < len(local) * nan_max:
        return True
    else:
        return False
    
    
def local_flowslope(local, peak):
    
    """
    Finds the flowslope and the deviation of the flowslope local to a given peak
    
    Parameters:
    -----------
    pandas.DataFrame/geopandas.GeoDataFrame: local
          Segment of full ground track local to the peak
    float/int: peak
          Location of the peak in question
          
    Returns:
    --------
    tuple: output
          (flowslope, flowslope_std)
          The flowslope, as well as the standard deviation of the flowslope in the local area.
          All operations are performed on the low pass filtered data
    """
    
    # sort dataframe by proximity to peak
    df_sort = local.iloc[(local['along_track_distance']-peak).abs().argsort()[:]]
    # compute paramters of nearest
    flowslope = df_sort.iloc[0]['slope-filt']
    flowslope_std = np.std(local['slope-filt'])
    return flowslope, flowslope_std


def comp_means(track, peak, avg_radius, avg_excl):
    
    """
    Computes the left and right hand mean of flowslope local to a peak
    
    Parameters:
    -----------
    pandas.DataFrame/geopandas.GeoDataFrame: track
          Ground track data
    int/float: peak
          Location of peak
    int/float: avg_radius
          Outward radius to compute the average under
    int/float: avg_excl
          Zone of exclusion around peak, to not include for the average calculation.
          This is important as it prevents grounding zone features from skewing the results as much.
          
    Returns:
    --------
    tuple: output
          (ahead_avg, behind_avg)
          ahead_avg is the average flowslope along track in the + direction. behind_avg is the same in the - direction.
    """
    
    along_max = peak + avg_radius
    along_min = peak - avg_radius
    slc = track[(track["along_track_distance"] > peak + avg_excl) & (track["along_track_distance"] < along_max)]["slope-filt"]
    ahead_avg = None if slc.empty == True else np.nanmean(slc)
    slc = track[(track["along_track_distance"] < peak - avg_excl) & (track["along_track_distance"] > along_min)]["slope-filt"]
    behind_avg = None if slc.empty == True else np.nanmean(slc)
    return ahead_avg, behind_avg


def comp_devs(track, peak, std_radius, std_excl):
    
    """
    Performs similarly to comp_means but instead with the standard deviation.
    
    Parameters:
    -----------
    pandas.DataFrame/geopandas.GeoDataFrame: track
          Ground track data
    int/float: peak
          Location of peak
    int/float: std_radius
          Outward radius to compute the deviation under
    int/float: std_excl
          Zone of exclusion around peak, to not include for the deviation calculation.
          This is important as it prevents grounding zone features from skewing the results as much.
          
    Returns:
    --------
    tuple: output
          (ahead_std, behind_std)
          ahead_std is the deviation of the flowslope along track in the + direction. behind_std is the same in the - direction.
    """
    
    along_max = peak + std_radius
    along_min = peak - std_radius
    slc = track[(track["along_track_distance"] > peak + std_excl) & (track["along_track_distance"] < along_max)]["slope-filt"]
    ahead_std = None if slc.empty == True else np.nanstd(slc)
    slc = track[(track["along_track_distance"] < peak - std_excl) & (track["along_track_distance"] > along_min)]["slope-filt"]
    behind_std = None if slc.empty == True else np.nanstd(slc)
    return ahead_std, behind_std



def filt_peaks(peak_dists, track,  debug=False):
    
    """
    Computes the quality score and other parameters to filter the peaks on.
    These parameters include the following: left & right standard deviation and average of flowslope, 
    flowslope at peak, standard deviation around peak, and the quality score.
    
    Parameters:
    -----------
    list/np.ndarray: peak_dists
          List of the along track locations where the peaks are located
    pandas.DataFrame/geopandas.GeoDataFrame: track
          Ground track data
    bool: debug
          Optional debug option. Prints helpful information for tweaking accuracy of the detection.
          Additionally can be used to understand why the program chose the peak it did.
          
    Returns:
    --------
    pandas.DataFrame: peakdf
          Pandas dataframe containing all the peaks and various parameters about them.
    """
    
    final_peaks, int_peaks = [], []
    peakdf = pd.DataFrame(columns=["loc", "slp", "std", "rht", "lft", "r_s", "l_s", "r_hs", "l_hs", "qs"])
    
    buffer = 0.5
    nan_max = 0.05
    avg_radius, avg_excl = 5, 0
    std_radius, std_excl = 40, 5
    
    # comp track params
    track_std = np.nanstd(np.array(track["slope"]))
    track_avg = np.nanmean(np.array(track["slope"]))
    
    atd_max, atd_min = np.max(track["along_track_distance"]), np.min(track["along_track_distance"])
    
    # test to see if each peak passes the requirements
    for i, peak in enumerate(peak_dists):
        local = track[(track["along_track_distance"] < peak + buffer) & (track["along_track_distance"] > peak - buffer)]
        
        # only keep peak if there is very few nans nearby.
        if nan_test(local, nan_max) and peak > atd_min + 50 and peak < atd_max - 50:# and rough_test(local, track_avg, track_std):
            # get flowslope and flowslope deviation near peak
            flowslope, flowslope_std = local_flowslope(local, peak)
            # compute mean both ahead & down track
            ahead_avg, behind_avg = comp_means(track, peak, avg_radius, avg_excl)
            # compute standard deviation up & down track
            ahead_std, behind_std = comp_devs(track, peak, std_radius, std_excl)
            hl_std, hr_std = None, None
            
            if type(ahead_std) != type(None) and type(behind_std) != type(None) and ahead_std != 0 and behind_std != 0:
                mag_dif_std = abs(math.log(ahead_std)-math.log(behind_std)) # promote a magnitude difference between standard deviations
                mag_dif_avg = abs(math.log(abs(ahead_avg))-math.log(abs(behind_avg))) # promote difference in mag between averages
                #mag_dif_hstd = abs(log(hr_std) - log(hl_std))
                dist_punish = abs(peak) / 100 # discourage giving points far from gline
                quality_score = flowslope_std * 200 + mag_dif_std + mag_dif_avg - dist_punish #+ mag_dif_hstd
            else:
                quality_score = None
            
            # add new peak to dataframe
            peakdf.loc[-1] = {"loc":peak, "slp":flowslope, "std":flowslope_std, "rht":ahead_avg, "lft":behind_avg, "r_s":ahead_std, "l_s":behind_std, "r_hs":hr_std, "l_hs":hl_std, "qs":quality_score}
            peakdf.index = peakdf.index + 1
            peakdf = peakdf.sort_index()
            if debug == True:
                if i == 0:
                    print("PICKS:")
                print(f"loc:{round(peak, 2)} slp:{round(flowslope, 4)} std:{round(flowslope_std, 6)} rht:{round(ahead_avg, 4)} lft:{round(behind_avg, 4)} r_s:{round(ahead_std, 6) if type(ahead_std) != type(None) else 'N/A'} l_s:{round(behind_std, 6) if type(behind_std) != type(None) else 'N/A'} qs:{round(quality_score, 4) if type(quality_score) != type(None) else 'N/A'}")
                #print({"loc":peak, "slp":flowslope, "std":flowslope_std, "rht":ahead_avg, "lft":behind_avg, "r_s":ahead_std, "l_s":behind_std})
        
        else:
            if debug == True:
                print(f"Track failed nan or roughness check")
                
    return peakdf


def reduce_peakdf(peakdf):
    
    """
    Reduces the data frame of possible break in slope options by restricting certain parameters based on what we know about ice sheets and ice shelves.
    
    Parameters:
    -----------
    pandas.DataFrame: peakdf
          Peak data
    float: d_min
          Minimum along track distance in reference to the grounding line intersection
    float: d_max
          Maximum along track distance in reference to the grounding line intersection
          
    Returns:
    --------
    pandas.DataFrame: reduced
          Reduced version of peakdf
    """
    
    buffer = 1
    
    std_min = 0.001
    slp_min = 0.0004
    
    shlf_min, shlf_max = -0.001, 0.002
    shet_max = -0.002
    
    reduced = peakdf[peakdf["qs"] > 2.5]
    reduced = reduced[(reduced["rht"] <= shet_max) & (reduced["lft"] >= shlf_min) & (reduced["lft"] < shlf_max) | (reduced["rht"] >= shlf_min) & (reduced["rht"] < shlf_max) & (reduced["lft"] <= shet_max)]
    
    return reduced
