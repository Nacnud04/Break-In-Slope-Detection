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
    
    proj4_crs = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    study_area = gpd.read_file(path)
    study_area = study_area.to_crs(proj4_crs)
    
    bounds = study_area["geometry"].total_bounds
    xlim, ylim = (bounds[0], bounds[2]), (bounds[1], bounds[3])
    
    return study_area, (xlim, ylim)


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return int(idx)


def returnRNC(track):
    row = track.iloc[0]
    rgt, name, cycle = row["rgt"], row["name"], row["cycle"]
    return rgt, name, cycle


def nan_test(local, nan_max):
    nan = np.count_nonzero(np.isnan(local["slope-filt"]))
    if nan < len(local) * nan_max:
        return True
    else:
        return False
    
    
def rough_test(track, avg, std):
    c = 0.975
    if c * len(track) <= len(track[(track["slope"] < avg + 2*std) & (track["slope"] > avg - 2*std)]):
        return True
    else:
        return False
    

def find_peaks(track, amp, dist):
    peak_index = np.array(signal.find_peaks(track["slope_deriv_1"]*-1, height=amp, distance=dist)[0])
    peak_index = np.append(peak_index, np.array(signal.find_peaks(track["slope_deriv_1"], height=amp, distance=dist)[0]))
    peak_dists = (((track["along_dist"].max() - track["along_dist"].min()) / len(track)) * peak_index) + track["along_dist"].min()
    return peak_dists
    
    
def local_flowslope(local, peak):
    df_sort = local.iloc[(local['along_dist']-peak).abs().argsort()[:]]
    flowslope = df_sort.iloc[0]['slope-filt']
    flowslope_std = np.std(local['slope-filt'])
    return flowslope, flowslope_std


def comp_hstd(track, radius, peak, rgt, name, cycle, gline_xy, debug=False):
    
    paths = []
    for cyc in range(2, 7):
        path = f"Saves/{rgt}/{name}/Bung-{cyc}.json"
        if os.path.isfile(path) == True:
            paths.append(path)
    
    hl_arr = list(track[(track["along_dist"] < peak) & (track["along_dist"] > peak - radius)]["h_li"])
    hr_arr = list(track[(track["along_dist"] > peak) & (track["along_dist"] < peak + radius)]["h_li"])
    
    dmin = min(list(track[(track["along_dist"] < peak) & (track["along_dist"] > peak - radius)]["along_track_dist"]))
    dmid = max(list(track[(track["along_dist"] < peak) & (track["along_dist"] > peak - radius)]["along_track_dist"]))
    dmax = max(list(track[(track["along_dist"] > peak) & (track["along_dist"] < peak + radius)]["along_track_dist"]))
    
    for path in paths:
        track = gpd.read_file(path)
        l = track[(track["along_track_dist"] < dmid) & (track["along_track_dist"] > dmin)]
        r = track[(track["along_track_dist"] > dmid) & (track["along_track_dist"] < dmax)]
        hl = list(l["h_li"])
        hr = list(r["h_li"])
        hl_arr.extend(hl)
        hr_arr.extend(hr)
        
    hl_std = np.nanstd(np.array(hl_arr))
    hr_std = np.nanstd(np.array(hr_arr))
    
    if debug:
        print(f"Hl std: {hl_std}  |  Hr std: {hr_std}")
    
    return hl_std, hr_std
        

def comp_means(track, peak, avg_radius, d_min, d_max):
    along_max = peak + avg_radius if peak + avg_radius < d_max else d_max
    along_min = peak - avg_radius if peak - avg_radius > d_min else d_min
    # catching warnings to disable notifications about taking mean of empty arrays
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        ahead_avg = np.nanmean(track[(track["along_dist"] > peak) & (track["along_dist"] < along_max)]["slope-filt"])
        behind_avg = np.nanmean(track[(track["along_dist"] < peak) & (track["along_dist"] > along_min)]["slope-filt"])
    return ahead_avg, behind_avg


def comp_devs(track, peak, std_radius, std_excl, d_min, d_max):
    along_max = peak + std_radius
    along_min = peak - std_radius
    slc = track[(track["along_dist"] > peak + std_excl) & (track["along_dist"] < along_max)]["slope-filt"]
    ahead_std = None if slc.empty == True else np.nanstd(slc)
    slc = track[(track["along_dist"] < peak - std_excl) & (track["along_dist"] > along_min)]["slope-filt"]
    behind_std = None if slc.empty == True else np.nanstd(slc)
    return ahead_std, behind_std
    

def filt_peaks(peak_dists, track, metadata, gline_xy, d_min, d_max, debug=False):
    
    final_peaks, int_peaks = [], []
    peakdf = pd.DataFrame(columns=["loc", "slp", "std", "rht", "lft", "r_s", "l_s", "qs"])
    
    buffer = 0.5
    nan_max = 0.05
    avg_radius = 5
    std_radius, std_excl = 40, 3
    
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
            ahead_avg, behind_avg = comp_means(track, peak, avg_radius, d_min, d_max)
            # compute standard deviation up & down track
            ahead_std, behind_std = comp_devs(track, peak, std_radius, std_excl, d_min, d_max)
            # compute deviation for h_li
            #hl_std, hr_std = comp_hstd(track, 30, peak, rgt, name, cycle, gline_xy, debug=debug)
            
            if type(ahead_std) != type(None) and type(behind_std) != type(None):
                mag_dif_std = abs(log(ahead_std)-log(behind_std)) # promote a magnitude difference between standard deviations
                mag_dif_avg = abs(log(abs(ahead_avg))-log(abs(behind_avg))) # promote difference in mag between averages
                dist_punish = abs(peak) / 100 # discourage giving points far from gline
                quality_score = flowslope_std * 200 + mag_dif_std + mag_dif_avg - dist_punish
            else:
                quality_score = None
            
            # add new peak to dataframe
            peakdf.loc[-1] = {"loc":peak, "slp":flowslope, "std":flowslope_std, "rht":ahead_avg, "lft":behind_avg, "r_s":ahead_std, "l_s":behind_std, "qs":quality_score}
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
                 

def reduce_peakdf(peakdf, d_min, d_max):
    
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
    track = util.process_and_extract(track, rgt, name, gline_xy, verbose=False, describe=False)
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
        return None
    
    # take deriv
    track = util.deriv_on_gpd(track)
    
    # find indicies where the track crosses the gline
    ints = [i for i in range(len(track) - 1) if track.iloc[i]["gline_dist"]/abs(track.iloc[i]["gline_dist"]) != track.iloc[i+1]["gline_dist"]/abs(track.iloc[i+1]["gline_dist"])]
    cuts = [int(ints[0]/2)] + [int((ints[i] + ints[i+1])/2) for i in range(len(ints)-1)] + [int((ints[-1] + len(track))/2)]
    cuts_mod = [0] + cuts + [max(cuts)]
    if cuts_mod[-1] == cuts_mod[-2]:
        cuts_mod.pop(-1)
    if debug:
        print(f"Splitting into {len(cuts_mod)-1} segments\nTrack ints @: {ints} & cuts @: {cuts_mod}")
        
    segments = [track.iloc[cuts_mod[n]:cuts_mod[n+1]] for n in range(len(cuts)-1)]
    
    peaks, qss = [], []
    
    for segment in segments:
        
        # find peaks in each segment
        peak_dists = find_peaks(segment, 0.003, 400)

        # filt peaks using data from all segments
        peakdf = filt_peaks(peak_dists, track, metadata, gline_xy, d_min=d_min, d_max=d_max, debug=debug)
        final_peaks = reduce_peakdf(peakdf, d_min=d_min, d_max=d_max)


        round_peaks = [str(round(peak, 2)) for peak in final_peaks["loc"]]
        if debug:
            print(f"REDUCED PICKS:\n {', '.join(round_peaks)}")


        if len(final_peaks) > 0:

            final_peaks = final_peaks.sort_values(by=['qs'], ascending=False)

            peak = final_peaks.iloc[0]["loc"]
            qs = final_peaks.iloc[0]["qs"]

            if debug:
                print(f"CHOICE:\n {peak}")
                
            peaks.append(peak)
            qss.append(qs)
        
        #peak = min(final_peaks, key=abs)
        #direc = direc[find_nearest(final_peaks, peak)]
        
    if len(peaks) > 0:
        return out_track, peaks, qss
    else:
        return out_track