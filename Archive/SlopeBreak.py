import pandas as pd
import geopandas as gpd
import numpy as np
import ICESat2GroundingLineMigration.IceSatHDF5Unpacker as unpack
import matplotlib.pyplot as plt
import src.TrackProfileUtil as util
import scipy.signal as signal
import warnings
pd.options.mode.chained_assignment = None


def studyArea(path):
    
    proj4_crs = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # proj4 string of the coord system to be used
    study_area = gpd.read_file(path)
    study_area = study_area.to_crs(proj4_crs)
    
    bounds = study_area["geometry"].total_bounds
    xlim, ylim = (bounds[0], bounds[2]), (bounds[1], bounds[3])
    
    return study_area, (xlim, ylim)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return int(idx)

def filtPeaks(peak_dists, track, d_min, d_max, debug=False):
    
    final_peaks = []
    int_peaks = []
    
    peakdf = pd.DataFrame(columns=["loc", "slp", "std", "rht", "lft", "r_s", "l_s"])
    slf = []
    avg_radius = 5
    std_radius = 40
    
    for peak in peak_dists:
        sel = track[(track["along_dist"] < peak + 0.5) & (track["along_dist"] > peak - 0.5)]
        nan = np.count_nonzero(np.isnan(sel["slope-filt"]))
        if nan < len(sel) * 0.05:
            df_sort = track.iloc[(track['along_dist']-peak).abs().argsort()[:]]
            flowslope_std = np.std(sel["slope-filt"]) 
            flowslope = df_sort.iloc[0]["slope-filt"]
            
            along_max = peak + avg_radius if peak + avg_radius < d_max else d_max
            along_min = peak - avg_radius if peak - avg_radius > d_min else d_min
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                ahead_avg = np.nanmean(track[(track["along_dist"] > peak) & (track["along_dist"] < along_max)]["slope-filt"])
                behind_avg = np.nanmean(track[(track["along_dist"] < peak) & (track["along_dist"] > along_min)]["slope-filt"])
            
            along_max = peak + std_radius if peak + std_radius < d_max else d_max
            along_min = peak - std_radius if peak - std_radius > d_min else d_min
            slc = track[(track["along_dist"] > peak + 2.5) & (track["along_dist"] < along_max)]["slope-filt"]
            if slc.empty:
                ahead_std = None
            else: 
                ahead_std = np.std(slc)
            slc = track[(track["along_dist"] < peak - 2.5) & (track["along_dist"] > along_min)]["slope-filt"]
            if slc.empty:
                behind_std = None
            else:
                behind_std = np.std(slc)
            
            peakdf.loc[-1] = {"loc":peak, "slp":flowslope, "std":flowslope_std, "rht":ahead_avg, "lft":behind_avg, "r_s":ahead_std, "l_s":behind_std}
            peakdf.index = peakdf.index + 1  # shifting index
            peakdf = peakdf.sort_index()  # sorting by index
            if debug == True:
                try:
                    print(f"loc: {peak} std: {round(flowslope_std, 6)} slp: {round(flowslope, 4)} rht: {round(ahead_avg, 4)} lft: {round(behind_avg, 4)} r_s: {round(ahead_std, 6)} l_s: {round(behind_std, 6)}")
                except TypeError:
                    print(f"loc: {peak} std: {round(flowslope_std, 6)} slp: {round(flowslope, 4)} rht: {round(ahead_avg, 4)} lft: {round(behind_avg, 4)}")
    
    shlf_min, shlf_max = -0.001, 0.002
    shet_max = -0.002
    
    reduced = peakdf[(abs(peakdf["std"]) > 0.001) & (abs(peakdf["slp"]) > 0.0004) & (peakdf["loc"] > d_min + 1) & (peakdf["loc"] < d_max - 1)]
    
    reduced = reduced[(reduced["rht"] <= shet_max) & (reduced["lft"] >= shlf_min) & (reduced["lft"] < shlf_max) | 
                      (reduced["rht"] >= shlf_min) & (reduced["rht"] < shlf_max) & (reduced["lft"] <= shet_max)]
             
    if reduced.empty:
        reduced = peakdf
        reduced = reduced[(reduced["rht"] <= shet_max) & (reduced["lft"] >= shlf_min) & (reduced["lft"] < shlf_max) | 
                          (reduced["rht"] >= shlf_min) & (reduced["rht"] < shlf_max) & (reduced["lft"] <= shet_max)]
        
    for index, row in reduced.iterrows():
        if row["rht"] <= -0.002 and row["lft"] >= -0.0013 and row["lft"] < 0.002:
            slf.append(-1)
        elif row["rht"] >= -0.001 and row["rht"] < 0.002 and row["lft"] < -0.002:
            slf.append(1)
        if debug:
            print(f"loc: {row['loc']} ice shelf side: {'+' if slf[-1] == 1 else '-'}")
    """
    
    reduced = reduced[(reduced["l_s"] < 0.0009) & (reduced["r_s"] > reduced["l_s"] + 0.0001) | (reduced["r_s"] < 0.0009) & (reduced["l_s"] > reduced["r_s"] + 0.0001)]
    if reduced.empty:
        reduced = peakdf
        reduced = reduced[(reduced["l_s"] < 0.0009) & (reduced["r_s"] > reduced["l_s"] + 0.0001) | (reduced["r_s"] < 0.0009) & (reduced["l_s"] > reduced["r_s"] + 0.0001)]
        
    for index, row in reduced.iterrows():
        if row["l_s"] < 0.0009 and row["r_s"] > row["l_s"] + 0.0001:
            slf.append(-1)
        elif row["r_s"] < 0.0009 and row["l_s"] > row["r_s"] + 0.0001:
            slf.append(-1)
        if debug:
            print(f"loc: {row['loc']} ice shelf side: {'+' if slf[-1] == 1 else '-'}")
    """
    
    final_peaks = np.array(reduced["loc"])  
                
    return final_peaks, slf

def findIb(track, gline_xy, debug=False):
    
    row = track.iloc[0]
    rgt, name, cycle = row["rgt"], row["name"], row["cycle"]
    track = util.process_and_extract(track, rgt, name, gline_xy, verbose=False, describe=False)
    if type(track) == type(None):
        return None
    
    out_track = track
    
    d_min = track["along_dist"].min() + 0.5
    d_max = track["along_dist"].max() - 0.5
    track = track[(track["along_dist"] > d_min) & (track["along_dist"] < d_max)]
    
    fidelity = 5000
    interped, interped_full, idx = util.interp_clean_single(track, fidelity, d_min, d_max)
    track = interped
    
    fill = 0.66
    nonnan = np.count_nonzero(~np.isnan(track["slope-filt"]))
    if nonnan < len(track["slope-filt"]) * fill:
        return None
    
    # take deriv
    track = util.deriv_on_gpd(track)
    
    """
    # split track by intersection with grounding line
    intersections = []
    for i in range(len(track) - 1):
        bef, aft = track.iloc[i], track.iloc[i+1]
        dist_bef, dist_aft = bef["gline_dist"], aft["gline_dist"]
        if dist_bef/abs(dist_bef) != dist_aft/abs(dist_aft):
            intersect = (bef["along_dist"] + aft["along_dist"]) / 2
            intersections.append(intersect)  
    splits = [(intersections[i]+intersections[i+1])/2 for i in range(len(intersections)-1)]
    #segments = [track[track["along_dist"] < splt] if i == 0 else track[track["along_dist"] > splt] if i == len(splits)-1 else track[(track["along_dist"] > splt)] for i, splt in enumerate(splits)]
    segments = []
    for i, splt in enumerate(splits):
        if i == 0:
            segments.append(track[track["along_dist"] < splt])
            segments.append(track[(track["along_dist"] > splt) & (track["along_dist"] < splits[i+1])])
        elif i == len(splits) - 1:
            segments.append(track[track["along_dist"] > splt])
        else:
            segments.append(track[(track["along_dist"] > splt) & (track["along_dist"] < splits[i+1])])
    
    if debug:
        colors = ["red", "blue", "green"]
        for i, segment in enumerate(segments):
            plt.plot(segment["along_dist"], segment["slope-filt"], c=colors[i])
        plt.show()
    
    # find peaks
    amp = 0.003
    final_peaks = []
    direcs = []
    for trk in segments:
        peak_index = np.array(signal.find_peaks(trk["slope_deriv_1"] * -1, height=amp, distance = 400)[0])
        peak_index = np.append(peak_index, np.array(signal.find_peaks(trk["slope_deriv_1"], height=amp, distance = 400)[0]))
        peak_dists = (((trk["along_dist"].max() - trk["along_dist"].min()) / len(trk)) * peak_index) + trk["along_dist"].min()
        filt_peaks, direc = filtPeaks(peak_dists, track, d_min=d_min, d_max=d_max, debug=debug)
        if len(filt_peaks) > 0:
            peak = min(filt_peaks, key=abs)
            direc = direc[find_nearest(filt_peaks, peak)]
            final_peaks.append(peak)
            direcs.append(direc)
    """
    
    amp = 0.003
    peak_index = np.array(signal.find_peaks(track["slope_deriv_1"] * -1, height=amp, distance = 400)[0])
    peak_index = np.append(peak_index, np.array(signal.find_peaks(track["slope_deriv_1"], height=amp, distance = 400)[0]))
    peak_dists = (((track["along_dist"].max() - track["along_dist"].min()) / len(track)) * peak_index) + track["along_dist"].min()
    final_peaks, direc = filtPeaks(peak_dists, track, d_min=d_min, d_max=d_max, debug=debug)
    
    if len(final_peaks) > 0:
        peak = min(final_peaks, key=abs)
        direc = direc[find_nearest(final_peaks, peak)]
        #final_peaks = final_peaks[~np.isnan(final_peaks)]
        #direc = direc[~np.isnan(direc)]
        #peak = final_peaks[0]
        #direc = int(direc[0])
        return out_track, peak, direc
    
    else:
        
        return out_track
