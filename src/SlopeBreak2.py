import pandas as pd
import geopandas as gpd
import numpy as np
import ICESat2GroundingLineMigration.IceSatHDF5Unpacker as unpack
import matplotlib.pyplot as plt
import src.TrackProfileUtil as util
import scipy.signal as signal
import warnings
import os
from pathlib import Path
from math import log
pd.options.mode.chained_assignment = None


def studyArea(path):
    
    """
    Takes in a filepath of a gpkg or similar file and returns a geodata frame as well as the geographical boundaries the file covers in EPSG:3031
    
    Parameters:
    -----------
    str: path
          This is the filepath of the gpkg or similar file
    
    Returns:
    --------
    geopandas.GeoDataFrame: study_area
          Geodataframe containing output data
    tuple: bounds
          Written as (xlim, ylim). Is a tuple of tuples where each xlim and ylim are (x/ymin, x/ymax).
    """
    
    # proj4 of EPSG:3031
    proj4_crs = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    study_area = gpd.read_file(path)
    study_area = study_area.to_crs(proj4_crs)
    # get x and y limits of study area in EPSG:3031
    bounds = study_area["geometry"].total_bounds
    xlim, ylim = (bounds[0], bounds[2]), (bounds[1], bounds[3])
    
    return study_area, (xlim, ylim)


def find_nearest(array, value):
    
    """
    Finds the index of the nearest value to a given value in an array
    
    Parameters:
    -----------
    list/np.ndarray/pandas.Series: array
          Array to search for nearest value in
    float/int: value
          Reference value
    
    Returns:
    --------
    int: idx
          Integer of index containing number nearest to input value
    """
    
    array = np.asarray(array)
    # subtract value from array and find lowest absolute value to get nearest.
    idx = (np.abs(array - value)).argmin()
    return int(idx)


def returnRNC(track):
    
    """
    Returns RGT, name, and cycle from a given track DataFrame/GeoDataFrame
    
    Parameters:
    -----------
    pandas.DataFrame/geopandas.GeoDataFrame: track
          Input track data
    
    Returns:
    --------
    tuple: rnc
          (rgt, name, cycle)
    """
    
    # assumes that the information in the first row of the df is constant throughout
    row = track.iloc[0]
    rgt, name, cycle = row["rgt"], row["name"], row["cycle"]
    return rgt, name, cycle


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


def rough_test(track, avg, std):
    
    """
    Poor test which checks if the percent of the data which is within 2 standard deviations of the mean, is greater than c (0.975).
    """
    
    c = 0.975 # percent within 2 standard deviations
    if c * len(track) <= len(track[(track["slope"] < avg + 2*std) & (track["slope"] > avg - 2*std)]):
        return True
    else:
        return False


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
    peak_dists = (((track["along_dist"].max() - track["along_dist"].min()) / len(track)) * peak_index) + track["along_dist"].min()
    return peak_dists


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
    df_sort = local.iloc[(local['along_dist']-peak).abs().argsort()[:]]
    # compute paramters of nearest
    flowslope = df_sort.iloc[0]['slope-filt']
    flowslope_std = np.std(local['slope-filt'])
    return flowslope, flowslope_std


def comp_hstd(track, radius, peak, rgt, name, cycle, gline_xy, debug=False):
    
    # this is matt's idea and it should be implemented more!!
    # only problem is it uses a shitttonnn of processing time and power, also
    # relies on tides and stuff,
    
    # original idea was to look at the standard deviation of all the data across time on each side of the peak.
    # this is not great as the ice sheet has a high variation in elevations, so a high standrd deviation
    
    # so, we could compare each point in the track with its corresponding point. except this is very computationally
    # expensive
    
    # so instead we take the average on each side at each time and compute the deviation of that relative to the mean
    
    """
    Takes a peak and returns the standard deviation of all heights across multiple cycles on either side of the peak.
    
    Parameters:
    -----------
    pandas.DataFrame/geopandas.GeoDataFrame: track
          Input ground track. (Only for 1 cycle!)
    float/int: radius
          Radius to look on either side of the peak
    float/int: peak
          Location of peak in question.
    str: rgt
          Name of the granule/rgt specifically
    str: name
          Name of the individual beam on ICESat-2 (eg. gt3l)
    str: cycle
          Number of the cycle in question
    np.ndarray: gline_xy
          Grounding line listed as an array of xy pairs.
    bool: debug
          Can turn on/off for more information.
          
    Returns:
    --------
    tuple: output
          (hl_std, hr_std)
          Left and right standard deviations of the elevations over time.
    """
    
    # check to see which cycles data exists for and put in list
    paths = []
    for cyc in range(1, 20):
        path = f"Saves/{rgt}/{name}/Bung-{cyc}.json"
        if os.path.isfile(path) == True and cyc != cycle:
            paths.append(path)
            
    l_means, r_means = [], []
    
    # compute for input array
    hl_arr = list(track[(track["along_dist"] < peak) & (track["along_dist"] > peak - radius)]["h_li"])
    hr_arr = list(track[(track["along_dist"] > peak) & (track["along_dist"] < peak + radius)]["h_li"])
    l_means.append(np.nanmean(hl_arr))
    r_means.append(np.nanmean(hr_arr))
    
    dmin = min(list(track[(track["along_dist"] < peak) & (track["along_dist"] > peak - radius)]["along_track_dist"]))
    dmid = max(list(track[(track["along_dist"] < peak) & (track["along_dist"] > peak - radius)]["along_track_dist"]))
    dmax = max(list(track[(track["along_dist"] > peak) & (track["along_dist"] < peak + radius)]["along_track_dist"]))
    
    # compute parameters for each cycle and append to list
    for path in paths:
        track = gpd.read_file(path)
        l = track[(track["along_track_dist"] < dmid) & (track["along_track_dist"] > dmin)]
        r = track[(track["along_track_dist"] > dmid) & (track["along_track_dist"] < dmax)]
        hl = list(l["h_li"])
        hr = list(r["h_li"])
        l_means.append(np.nanmean(hl))
        r_means.append(np.nanmean(hr))
        #hl_arr.extend(hl)
        #hr_arr.extend(hr)
        
    # perform calculations on conglomerated means and deviations
    hl_std, hr_std = np.nanstd(l_means), np.nanstd(r_means)
    hl_avg, hr_avg = np.nanmean(l_means), np.nanmean(r_means)
    
    std_rl = hl_std / hl_avg
    std_rr = hr_std / hr_avg
    
    if debug:
        print(f"Hl std: {std_rl}  |  Hr std: {std_rr}")
    
    return std_rl, std_rr


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
    slc = track[(track["along_dist"] > peak + avg_excl) & (track["along_dist"] < along_max)]["slope-filt"]
    ahead_avg = None if slc.empty == True else np.nanmean(slc)
    slc = track[(track["along_dist"] < peak - avg_excl) & (track["along_dist"] > along_min)]["slope-filt"]
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
    slc = track[(track["along_dist"] > peak + std_excl) & (track["along_dist"] < along_max)]["slope-filt"]
    ahead_std = None if slc.empty == True else np.nanstd(slc)
    slc = track[(track["along_dist"] < peak - std_excl) & (track["along_dist"] > along_min)]["slope-filt"]
    behind_std = None if slc.empty == True else np.nanstd(slc)
    return ahead_std, behind_std


def filt_peaks(peak_dists, track, metadata, gline_xy, debug=False):
    
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
    tuple: metadata
          Tuple of (rgt, name, cycle)
    np.npdarray: gline_xy
          List of XY pairs which make up the grounding line. This is produced by a function in TrackProfileUtil
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
    
    rgt, name, cycle = metadata
    
    if debug:
        print(f"TRACK PARAMS FOR {rgt}-{name}-{cycle}:")
        print(f"avg: {track_avg}")
        print(f"2nd std dev: {track_avg - 2*track_std} <-> {track_avg + 2*track_std}")
    
    # test to see if each peak passes the requirements
    for i, peak in enumerate(peak_dists):
        local = track[(track["along_dist"] < peak + buffer) & (track["along_dist"] > peak - buffer)]
        
        # only keep peak if there is very few nans nearby.
        if nan_test(local, nan_max):# and rough_test(local, track_avg, track_std):
            # get flowslope and flowslope deviation near peak
            flowslope, flowslope_std = local_flowslope(local, peak)
            # compute mean both ahead & down track
            ahead_avg, behind_avg = comp_means(track, peak, avg_radius, avg_excl)
            # compute standard deviation up & down track
            ahead_std, behind_std = comp_devs(track, peak, std_radius, std_excl)
            # compute deviation for h_li
            #hl_std, hr_std = comp_hstd(track, 30, peak, rgt, name, cycle, gline_xy, debug=debug)
            hl_std, hr_std = None, None
            
            if type(ahead_std) != type(None) and type(behind_std) != type(None) and ahead_std != 0 and behind_std != 0:
                mag_dif_std = abs(log(ahead_std)-log(behind_std)) # promote a magnitude difference between standard deviations
                mag_dif_avg = abs(log(abs(ahead_avg))-log(abs(behind_avg))) # promote difference in mag between averages
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


def ice_plain_check(track, peakdf, debug=False):
    
    """
    Attempts to check if an ice plain is present between two break in slope picks. As of 7/21/23 this needs substantial improvement.
    
    Parameters:
    -----------
    pandas.DataFrame/geopandas.GeoDataFrame: track
          Ground track data
    pandas.DataFrame: peakdf
          Dataframe containing each peak location and each peaks parameters
    bool: debug
          Provides additional debug information if set to true.
    
    Returns:
    --------
    tuple: output
          (loc_1, loc_2)
          Along track distance at the start of the predicted ice plain and the end of the predicted ice plain.
    """
    
    if len(peakdf) < 2:
        return None
    
    # first find the minimum of the heights
    # ice plain appears immediately before ice sheet minimum.
    min_id = np.argmin(track["h_li"])
    min_loc = track.iloc[min_id]["along_dist"]
    
    sort = peakdf.sort_values(by=['qs'], ascending=False) # sort peak picks by quality score
    sort = sort[sort["qs"] > -1.5] # only keep quality scores above a certain value
    
    max_plane_dist = 75
    min_plane_dist = 10
    sort["dist_to_best"] = abs(sort["loc"] - sort.iloc[0]["loc"]) # find the along track distane from each peak to the best peak
    sort = sort.sort_values(by=['dist_to_best'], ascending=True) 
    sort.reset_index(drop=True, inplace=True)
    
    devs = 2
    
    for i in range(len(sort)):
        
        # grab the two different peak locations. this process is done starting with the largest distance btwn loc1 and loc2
        loc1, loc2 = sorted([sort.iloc[0]["loc"], sort.iloc[-1 * i]["loc"]])
        loc1, loc2 = loc1 + 2, loc2 - 2 # implement 2km buffer to limit weird sheet/iceplain coupling from affecting the alg
        
        # sel should always be of length greater than 0. For some reason it isnt. So we skip over when it isnt.
        sel = track[(track["along_dist"] > loc1) & (track["along_dist"] < loc2)]
        if len(sel) == 0:
            continue
            
        sel_mean = np.nanmean(sel["slope-filt"])
        sel_std = np.nanstd(sel["slope-filt"])
        sel_upper, sel_lower = sel_mean + devs * sel_std, sel_mean - devs * sel_std
        if debug:
            print(f"btwn {loc1} and {loc2} upper is {sel_upper} lower is {sel_lower}")
        if sel_lower > -0.003 and sel_mean < 0.005:
            if debug:
                print(f"Ice plain predicted.")
            return (loc1 - 2, loc2 + 2)
    
    return None

def reduce_peakdf(peakdf, d_min, d_max):
    
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
    
    #reduced = peakdf[(abs(peakdf["std"]) > std_min) & (abs(peakdf["slp"]) > slp_min) & (peakdf["loc"] > d_min + buffer) & (peakdf["loc"] < d_max - buffer)] 
    
    reduced = peakdf[peakdf["qs"] > 2.5]
    reduced = reduced[(reduced["rht"] <= shet_max) & (reduced["lft"] >= shlf_min) & (reduced["lft"] < shlf_max) | (reduced["rht"] >= shlf_min) & (reduced["rht"] < shlf_max) & (reduced["lft"] <= shet_max)]
    
    return reduced


def findIb(track, gline_xy, debug=False):
    
    rgt, name, cycle = returnRNC(track)
    metadata = (str(int(rgt)), name, str(cycle))
    
    # ensures removal of weak/bad tracks
    track = util.process_and_extract(track, rgt, name, gline_xy, debug=debug, verbose=False, describe=False)
    if type(track) == type(None):
        return None
    
    # track data to export.
    out_track = track
    
    # test with and without this. There is code preventing any peaks from being in the outer kilometer of the track. So this should not have any effect on the outcome.
    d_min = track["along_dist"].min() + 0.5
    d_max = track["along_dist"].max() - 0.5
    track = track[(track["along_dist"] > d_min) & (track["along_dist"] < d_max)]
    
    # interpolate
    fidelity = 5000
    interped, interped_full, idx = util.interp_clean_single(track, fidelity, d_min, d_max)
    track = interped
    
    # remove track if too many nan
    fill = 0.66
    nonnan = np.count_nonzero(~np.isnan(track["slope-filt"]))
    if nonnan < len(track["slope-filt"]) * fill:
        if debug:
            print(f"Nan count of track too high {nonnan} < {len(track['slope-filt']) * fill}")
        return None
    
    # take deriv
    track = util.deriv_on_gpd(track)
    
    slope_breaks = [0.003, 0.002, 0.001, 0.0005, 0.0001]
    for thresh in slope_breaks:
        # find peaks
        peak_dists = find_peaks(track, thresh, 400)
        if len(peak_dists) > 3:
            break
    
    # filt peaks
    peakdf = filt_peaks(peak_dists, track, metadata, gline_xy, debug=debug)
    
    # attempt to find an ice plain
    try:
        plain = ice_plain_check(track, peakdf, debug=debug)
    except IndexError:
        plain = None
    
    # remove certain peak predictions based on individual values (NOT QUALITY SCORE!)
    final_peaks = reduce_peakdf(peakdf, d_min=d_min, d_max=d_max)
    
    # round peak values for display
    round_peaks = [str(round(peak, 2)) for peak in final_peaks["loc"]]
    if debug:
        print(f"REDUCED PICKS:\n {', '.join(round_peaks)}")
    
    
    if len(final_peaks) > 0:
        
        # choose the peak with the highest quality score and output.
        final_peaks = final_peaks.sort_values(by=['qs'], ascending=False)
        
        peak = final_peaks.iloc[0]["loc"]
        qs = final_peaks.iloc[0]["qs"]
        
        if debug:
            print(f"CHOICE:\n {peak}")
        
        #peak = min(final_peaks, key=abs)
        #direc = direc[find_nearest(final_peaks, peak)]
        
        return out_track, peak, qs, plain
    else:
        return out_track
