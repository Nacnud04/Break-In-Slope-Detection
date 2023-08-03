import geopandas as gpd
import numpy as np
import pandas as pd
import shapely as sp
import matplotlib.pyplot as plt
import scipy.signal as signal
import scipy.interpolate as interp

def gline_to_df(gpd, mask):
    
    """
    Transforms a GeoDataFrame of the grounding line. Specifically it crops, rearranges,
    and only keeps necessary data.
    
    Parameters:
    -----------
    gpd.GeoDataFrame: gpd
          GeoDataFrame containing the grounding line polygon
    gpd.GeoDataFrame: mask
          Mask of study area. Is used to crop the grounding line.
          NEEDS TO BE IN THE SAME CRS
    
    Returns:
    --------
    gpd.GeoDataFrame: gline_xy
          Output dataframe of the grounding line.
          Not georeferenced. Only x & y.
    """
    
    bounds = mask["geometry"].total_bounds
    xlim, ylim = (bounds[0], bounds[2]), (bounds[1], bounds[3])
    
    xsize, ysize = int((xlim[1] - xlim[0]) / 450), int((ylim[1] - ylim[0]) / 450)

    gline_clip = gpd.clip_by_rect(xlim[0], ylim[0], xlim[1], ylim[1])
    gline_poly = gline_clip.iloc[0] # only 1 polygon in the gdf, so extract it from the dataframe.

    x, y = gline_poly.exterior.coords.xy
    x, y = np.array(x.tolist()), np.array(y.tolist())

    x = np.delete(x, [0, -1, -2]) # removes lines on the border of the clipped region
    y = np.delete(y, [0, -1, -2])

    gline_xy = pd.DataFrame({"x":x,"y":y})
    
    return gline_xy


def interpolate(dataframe, column):
    
    """
    Performs a linear interpolation to fill in missing/nan data in a column of a dataframe
    
    Parameters:
    -----------
    pd.DataFrame: dataframe
          Dataframe to perform interpolation on
    str: column
          Key to column in dataframe to operate on
          
    Returns:
    --------
    pd.DataFrame: dataframe
          Dataframe with filled in column
    """
    
    function = interp.interp1d(dataframe.index, dataframe[column], kind = "linear")
    dataframe[column] = function(dataframe.index)
    return dataframe

def apply_butter(series, order, cut_off):
    
    """
    Takes in a series of data and applies a low pass butterworth filter
    
    Parameters:
    -----------
    (g)pd.Series: series
          Series of data for the filter to be applied on
    int: order
          Order of the filter. Can be interpreted as the rate at which the
          filters response falls in the transition band
    float: cutoff
          Cutoff frequency. Is a ratio. f/f0
          
    Returns:
    --------
    pd.Series: series
          Filtered output of series data
    """
    
    b, a = signal.butter(order, cut_off, btype="lowpass")
    series = signal.filtfilt(b, a, np.array(series))
    return series


def extract_data(data, rgt, name, cycle):
    
    """
    Extracts a specific rgt and laser from a dataframe produced by VisualSlopes.ipynb
    
    Parameters:
    -----------
    gpd.GeoDataFrame: data
          Geodataframe containing all data produced by VisualSlopes.ipynb
    int: rgt
          Region ground track number desired
    str: name
          Name of laser to be extracted. Examples include "gt3r"
          
    Returns:
    --------
    gpd.GeoDataFrame: single_beam
          Geodataframe of just the selected beams data
    """
    
    # extract data
    if cycle:
        data = data[data["cycle"] == str(cycle)]
    if rgt:
        data = data[data["rgt"] == rgt]
    if name:
        data = data[data["name"] == name]
    return data
    

def compute_along_track_dist(single_beam, verbose=True):
    
    """
    Computes along track distance for a given beam's geodataframe
    
    Parameters:
    -----------
    gpd.GeoDataFrame: single_beam
          Geodataframe containing data for a single beam outputted by the extract_data method above
    bool: verbose
          Default to true, prints some helpful outputs for debugging
    
    Returns:
    --------
    gpd.GeoDataFrame: single_beam
          This is identical to the input dataframe, except with an output column {along_dist}
    """
    
    # compute along track distance
    ys = []
    xs = []
    for index, row in single_beam.iterrows():
        x, y = row["geometry"].xy
        xs.append(x[0])
        ys.append(y[0])
    single_beam["y"] = ys
    single_beam["x"] = xs
    date_time = np.array(single_beam["date_time"])
    single_beam = single_beam.set_index("date_time")
    single_beam["date_time"] = date_time
    single_beam.sort_index(inplace=True)
    
    track_maxx, track_minx = single_beam["x"].max(), single_beam["x"].min()
    track_maxy, track_miny = single_beam["y"].max(), single_beam["y"].min()
    
    distances = []
    i = 0
    for index, row in single_beam.iterrows():
        if len(distances) == 0:
            distances.append(0)
        else:
            x_1, y_1 = single_beam["geometry"].iloc[i-1].xy
            x_1, y_1 = x_1[0], y_1[0]
            x_2, y_2 = row["geometry"].xy
            x_2, y_2 = x_2[0], y_2[0]
            dist = ((x_2 - x_1)**2 + (y_2 - y_1)**2)**0.5
            distances.append(dist + distances[-1])
        i += 1
    single_beam["along_dist"] = distances
    
    x_1, y_1 = single_beam["geometry"].iloc[0].xy
    x_1, y_1 = x_1[0], y_1[0]
    x_2, y_2 = single_beam["geometry"].iloc[-1].xy
    x_2, y_2 = x_2[0], y_2[0]
    
    if verbose:
        print("-Distance Check-")
        print(f"Estimated Dist: {distances[-1]} m")
        print(f"Rough Dist: {((x_2 - x_1)**2 + (y_2 - y_1)**2)**0.5} m")

    single_beam["along_dist"] = single_beam["along_dist"] / 1000
    
    return single_beam


def find_gline_dist(single_beam, gline):
    
    """
    Finds the distance to the grounding line for every point in single_beam
    
    Parameters:
    -----------
    gpd.GeoDataFrame: single_beam
          This is the dataframe either extracted from extract_data() or returned 
          by compute_along_track_dist() or find_gline_int()
    sp.geometry.LinesString: gline
          Shapely linestring of the polygon as was returned by find_gline_int()
          
    Returns:
    --------
    gpd.GeoDataframe: single_beam
          Same as single_beam input, except with the new {gline_dist} column
    """
    
    gline_dist = []
    
    for index, row in single_beam.iterrows():
        # compute distance to gline for each point
        x, y = row["geometry"].xy
        point = sp.geometry.Point(x[-1], y[-1])
        out = sp.ops.nearest_points(point, gline)
        point, gline_nearest = out
        g_x, g_y = gline_nearest.xy
        p_x, p_y = point.xy
        # uncomment the following if needed.
        #gline_dist.append(point.distance(gline_nearest))
        if p_y[0] < g_y[0]: # optimally this would be done by checking if the point is on the ice sheet. this is not possible however as we are looking at a segment of the grounding line to boost processing speed. so, after this step, numerous errors are checked and corrected for.
            gline_dist.append(-1*point.distance(gline_nearest))
        else:
            gline_dist.append(point.distance(gline_nearest))
        
            
    single_beam["gline_dist"] = np.array(gline_dist) / 1000
    
    return single_beam



def find_gline_int(single_beam, gline_xy, verbose=True):
    
    """
    Locates the intersection of the single_beam geodataframe and the gline_xy dataframe
    
    Parameters:
    -----------
    gdf.GeoDataFrame: single_beam
          This is the dataframe either extracted from extract_data() or returned by compute_along_track_dist()
    pd.DataFrame: gline_xy
          The grounding line x & y values in EPSG:3031 as returned by gline_to_df
    bool: verbose
          Prints some helpful information for debugging
          
    Returns:
    --------
    gdf.GeoDataFrame: single_beam
          Identical to single_beam input except with the "along_dist" column centered around
          the grounding line
    sp.Geometry.LineString: gline
          Shapely linestring of the grounding line.
    """
    
    track_maxx, track_minx = single_beam["x"].max(), single_beam["x"].min()
    track_maxy, track_miny = single_beam["y"].max(), single_beam["y"].min()
    
    # Hunt for intersect w/ gline
    
    # prepare gline
    gline_xy = gline_xy[(gline_xy["x"] < track_maxx) & (gline_xy["x"] > track_minx)]
    gline_xy = gline_xy[(gline_xy["y"] < track_maxy) & (gline_xy["y"] > track_miny)]
    
    gline_x, gline_y = np.array(gline_xy["x"]), np.array(gline_xy["y"])
    gline_delt = (np.diff(gline_x)**2 + np.diff(gline_y)**2)**0.5
    breaks = [i for i, delt in enumerate(gline_delt) if delt > 10000]
    
    xs = np.split(gline_x, breaks)
    ys = np.split(gline_y, breaks)
    
    if len(xs) != 0:
        segments = []
        for i in range(len(xs)): 
            if len(xs[i]) >= 2:
                gline_np = np.vstack((xs[i], ys[i])).T
                segments.append(sp.geometry.LineString(gline_np))
    else:
        if verbose:
            print("Cannot form into line segment. Gline intersection algorithm failure")
        return None
    
    if len(gline_x) < 2:
        return None
    
    full_gline = sp.geometry.LineString(np.vstack((gline_x, gline_y)).T)
    
    track_maxx, track_minx = single_beam["x"].max(), single_beam["x"].min()
    track_maxy, track_miny = single_beam["y"].max(), single_beam["y"].min()
    track_np = np.vstack((single_beam["x"], single_beam["y"]))
    track = sp.geometry.LineString(track_np.T)
    
    for segment in segments:
        intersect = segment.intersection(track)
        if type(intersect) == sp.geometry.Point:
            x_int, y_int = intersect.xy
            x_int, y_int = x_int[0], y_int[0]
            gline = segment
            break
        elif type(intersect) == sp.geometry.MultiPoint:
            x_int, y_int = intersect.geoms[0].xy
            x_int, y_int = x_int[0], y_int[0]
            gline = segment
            break
        elif type(intersect) == sp.geometry.LineString: 
            x_int, y_int = False, False
            if verbose:
                print(f"No intersection found!")
                #print(find_gline_dist(single_beam, full_gline))
            #return None
            
    if len(segments) < 1:
        return None
            
    if x_int == False:
        return None
        
    if verbose:
        print(f"Found intersection at: {x_int},{y_int}")
        
    track_np = track_np.T
    
    direc_x, direc_y = [], []
    
    i = 0
    while i < len(track_np) - 1:
        if track_np[i, 0] < track_np[i+1, 0]:
            direc_x.append(1)
        else:
            direc_x.append(-1)
        if track_np[i, 1] < track_np[i + 1, 1]:
            direc_y.append(1)
        else:
            direc_y.append(-1)                
        i += 1

    direc_x.append(direc_x[-1])
    direc_y.append(direc_y[-1])

    single_beam["direc_x"] = direc_x
    single_beam["direc_y"] = direc_y

    # get closest point
    out = sp.ops.nearest_points(sp.geometry.MultiPoint(track_np), gline)
    point, gline_nearest = out
    row = single_beam[single_beam["x"] == point.x].iloc[0]

    single_beam_clipped = single_beam

    if row["direc_x"] == 1:
        single_beam_clipped = single_beam_clipped[single_beam_clipped["x"] < x_int]
    else:
        single_beam_clipped = single_beam_clipped[single_beam_clipped["x"] > x_int]

    if row["direc_y"] == 1:
        single_beam_clipped = single_beam_clipped[single_beam_clipped["y"] < y_int]
    else:
        single_beam_clipped = single_beam_clipped[single_beam_clipped["y"] > y_int]

    last_point = single_beam_clipped.iloc[-1]

    int_dist = (((last_point["x"] - x_int)**2 + (last_point["y"] - y_int)**2)**0.5) / 1000 + row["along_dist"]
    
    if verbose:
        print(f"Along track dist @ intersection {int_dist}")
    
    # offset the along track distance, to be 0 at the gline.
    single_beam["along_dist"] = single_beam["along_dist"] - int_dist
    
    return single_beam, gline


def find_max_min(list_of_df, column):
    
    """
    Takes a list of dataframes and find the maximum and minimum value in a specified
    column across all of the dataframes.
    
    Parameters:
    -----------
    list: list_of_df
          List of pandas dataframes
    string: column
          Column in dataframe to find the absolute max and min for
          
    Returns:
    --------
    float: maximum
          Maximum value
    float: minimum
          Minimum value
    """
    
    maximum, minimum = 0, 0
    for df in list_of_df:
        if df[column].max() > maximum:
            maximum = df[column].max()
        if df[column].min() < minimum:
            minimum = df[column].min()
    return maximum, minimum

def process_and_extract(data, rgt, name, gline_xy, debug=False, describe = True, verbose=True):
    
    """
    Extracts desired rgt ground track and cycle from a dataframe, then computes the
    along track distance and intersection. Then it offsets so the intersection
    distance is at an along track distance of 0. Lastly it computes the distance to
    the gline at every point along the track.
    
    Parameters
    ----------
    pandas.DataFrame: data
          dataframe containing all of the ground track data
    string/int: rgt
          desired rgt number
    string: name
          name of desired ground track. (eg. gt1l)
    np.ndarray: gline_xy
          array of xy points which correspond to the grounding line in the active
          coordinate system.
          this should be cropped to the study region
    boolean: describe
          turns on/off statements of which step in the process the function is at
    boolean: verbose
          turns on/off additional information
          
    Returns
    -------
    pd.DataFrame: single_beam
          updated dataframe
    """
    
    # extract data
    if describe:
        print(f"Extracting data     ", end="\r")
    single_beam = extract_data(data, rgt, name, cycle=None)
    
    # compute along track distance
    if describe:
        print(f"Comp along track dist    ", end="\r")
    single_beam = compute_along_track_dist(single_beam, verbose)
    
    # Hunt for intersect w/ gline
    if describe:
        print(f"Finding intersection    ", end="\r")
    out = find_gline_int(single_beam, gline_xy, debug)
    
    if type(out) == type(None):
        return None
    
    single_beam, gline = out
    
    # compute distance to the grounding line
    if describe:
        print(f"Computing gline dist    ", end="\r")
    single_beam = find_gline_dist(single_beam, gline)
    
    order = 5
    cutoff = 0.032 # found by https://doi.org/10.5194/tc-14-3629-2020
    try:
        single_beam["slope-filt"] = apply_butter(single_beam["slope"], order, cutoff) # uses scikit to apply the filter.
    except ValueError:
        # The length of the input vector x must be greater than padlen, which is 18.???
        return None
    
    return single_beam

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
    
    slope_raw = np.interp(idx, track["along_dist"], track["slope"])
    slope_filt = np.interp(idx, track["along_dist"], track["slope-filt"])
    h_li = np.interp(idx, track["along_dist"], track["h_li"])
    gline_dist = np.interp(idx, track["along_dist"], track["gline_dist"])
    along_track_dist = np.interp(idx, track["along_dist"], track["along_track_dist"])
    df = pd.DataFrame(index = idx, data = {"along_dist":idx, "along_track_dist":along_track_dist, "slope":slope_raw, "slope-filt":slope_filt, "h_li":h_li, "gline_dist":gline_dist, "cycle":[track.iloc[0]["cycle"]] * len(idx)})
    
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
    
    along_dist = np.array(track["along_dist"])
    delta_dist = np.diff(along_dist)
    for j, delt in enumerate(delta_dist):
        if delt > (max_dist / 1000) and j < len(delta_dist):
            min_idx, max_idx = along_dist[j], along_dist[j + 1]
            df.loc[(df["along_dist"] > min_idx) & (df["along_dist"] < max_idx), ["slope-filt", "h_li"]] = None
            #print(f"Removed data from {min_idx}km to {max_idx}km of delta {delt * 1000}m", end = "           \r")
            count += 1
    return df, count

def interp_clean(cropped, fidelity, d_min, d_max):
    
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
    
    idx = np.linspace(d_min, d_max, fidelity)
    interped = [[],[],[]]
    interped_full = [[],[],[]]

    dx = (d_max - d_min) / fidelity

    count = 0

    for i, crop in enumerate(cropped):
        for track in crop:

            row = track.iloc[0]
            name, rgt, cycle = row["name"], row["rgt"], row["cycle"]

            interped_full[i].append(interpolation(idx, track))

            # remove points which were interpolated over large chunks of missing data
            max_dist = 40 # value in meters

            df = interpolation(idx, track)
            df, count = clean_by_gaps(track, df, max_dist, count)

            # remove points below a certain standard deviation threshold
            for dist in reversed(range(d_min * 1000, d_max * 1000, 100)):
                dist /= 1000
                sel = df[df["along_dist"] < dist]
                if sel["slope-filt"].std() < 1e-10:
                    df.loc[df["along_dist"] < dist, ["slope-filt", "h_li"]] = None
                    break
                
            #print(f"Deviation of {rgt}-{name}-{cycle}: {df['slope-filt'].std()}", end = "                                                         \r")
            if df["slope-filt"].std() > 1e-10:
                interped[i].append(df)

    print(f"\nRemoved {count} gaps in data")
    
    return interped, interped_full, idx

def interp_clean_single(track, fidelity, d_min, d_max):
    
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
    
    idx = np.linspace(d_min, d_max, fidelity)

    dx = (d_max - d_min) / fidelity

    count = 0

    row = track.iloc[0]
    name, rgt, cycle = row["name"], row["rgt"], row["cycle"]

    df_full = interpolation(idx, track)

    # remove points which were interpolated over large chunks of missing data
    max_dist = 40 # value in meters

    df = interpolation(idx, track)
    df, count = clean_by_gaps(track, df, max_dist, count)

    # remove points below a certain standard deviation threshold
    for dist in reversed(range(int(d_min * 1000), int(d_max * 1000), 100)):
        dist /= 1000
        sel = df[df["along_dist"] < dist]
        if sel["slope-filt"].std() < 1e-10:
            df.loc[df["along_dist"] < dist, ["slope-filt", "h_li"]] = None
            break

    #print(f"Deviation of {rgt}-{name}-{cycle}: {df['slope-filt'].std()}", end = "                                                         \r")
    if df["slope-filt"].std() < 1e-10:
        return None
    
    return df, df_full, idx

import scipy.signal as signal
def corr_lag(interped, d_max, d_min, fidelity):
    
    """
    Calculates the correlation lag for all of the different 
    cycles along each ground track
    
    Parameters
    ----------
    list: interped
          2d array containing the interpolated output from the
          interp clean function
    int: d_max
          maximum value of the interpolated data
    int: d_min
          minimum value of the interpolated data
    int: fidelity
          amount of total points in the interpolated data
          
    Returns
    -------
    list: lag_km
          2d array containing the correlation lag in km
    list: final_lag
          2d array containing the correlation lag in indicies
    """
    
    final_lag = [[], [], []]
        
    for i, beam in enumerate(interped):
        baseline = beam[0]["slope-filt"]
        for track in beam:
            dat = np.nan_to_num(track["slope-filt"], nan = np.average(track["slope-filt"]))
            correlation = signal.correlate(baseline, dat, mode='full')
            lags = signal.correlation_lags(len(baseline), len(track), mode="full")
            final_lag[i].append(lags[np.argmax(np.abs(correlation))])

    print(f"Lags in terms of indices: {final_lag}")

    med_delt_dist = (d_max - d_min) / fidelity

    lag_km = []
    for i, lag_arr in enumerate(final_lag):
        lag_km.append(np.array(lag_arr) * med_delt_dist)
    lag_km = np.array(lag_km)

    print(f"Lags in meters: {lag_km * 1000}")
    
    return lag_km, final_lag

def comp_deriv(series, dist):
    derivs = []
    for i in range(len(series) - 1):
        deriv = (series.iloc[i + 1] - series.iloc[i]) / (dist.iloc[i + 1] - dist.iloc[i])
        derivs.append(deriv)
    derivs.append(derivs[-1])
    return np.array(derivs)


def deriv_on_gpd(gpd):
    gpd["slope_deriv_1"] = comp_deriv(gpd["slope-filt"], gpd["along_dist"])
    order = 5
    cutoff = 0.3
    gpd[f"slope_deriv_1-filt"] = apply_butter(gpd["slope_deriv_1"], order, cutoff)
    return gpd