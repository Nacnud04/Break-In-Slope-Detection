import IceSatHDF5Unpacker as unpack
import FlowUnpacker as flow
import GLineUnpacker as gline
import Visualizations as visualize
import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.interpolate import griddata

import matplotlib.pyplot as plt

print("Importing dataset")
dataset = unpack.Dataset(r"C:\Users\Nacnu\Documents\IceSat2")

print("Importing grounding line data")
groundline = gline.GLineData("grounding_line_data/InSAR_GL_Antarctica_v02.shp", 50000)

polygons = groundline.linepolygon.set_crs(crs="+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

print("Grabbing filenames near polygon")
files = np.array(dataset.returnFilenamesByPolygon(polygons))
files = np.unique(files)

print(f"Original filecount: {len(dataset.datafiles)}  |  New filecount: {len(files)}")

# Import the ice flow database
print("Importing and computing flow")
flowdatabase = flow.Database("data/antarctic_ice_vel_phase_map_v01.nc")

# Compute all angles and errors from the flow database
flowdatabase.compute_all()

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def plotTrackSlope():
    xlim, ylim = (-400000,200000), (-750000,-300000)
    latlim, lonlim = (-86, -81)
    i = 0
    fill = 40
    window = 10
    allx, ally, allslope = np.array([]), np.array([]), np.array([])
    for filename in files:
        print(f"Processing granule {filename}")
        granule = dataset.openfilename(filename)
        lasercount = 0
        for laser in granule.lasers:
            print(f"laser {lasercount}")
            # grab all necessary values from laser


            lat, lon, time, dh_fit_dx, dh_fit_dx_sigma, metadata, granuledata = laser.getTrackData()
            dh_fit_dy = laser.returnAcrossTrackSlope()
            print(1)
            xyindices, trackx, tracky = unpack.Basemap.latlonindex_to_grid(lat, lon, flowdatabase)
            print(2)

            # transform into a pandas dataframe
            laserdf = pd.DataFrame({"x":trackx, "y":tracky, "dh_fit_dx":dh_fit_dx, "dh_fit_dy":dh_fit_dy, "dh_fit_dx_sigma":dh_fit_dx_sigma, "xyindices":xyindices})

            # remove points outside of wanted area
            laserdf = laserdf[laserdf["x"] > xlim[0]]  # deal with x values
            laserdf = laserdf[laserdf["x"] < xlim[1]]
            laserdf = laserdf[laserdf["y"] > ylim[0]]  # deal with y values
            laserdf = laserdf[laserdf["y"] < ylim[1]]
        
            # Rearrange flow angle and flow angle error in the order of the xy indices
            flow_angle, flow_angle_error = flowdatabase.rearrange_angle_data(laserdf["xyindices"])
            laserdf["flow_angle"], laserdf["flow_angle_error"] = flow_angle, flow_angle_error

            # Get the projected azumith at each x and y point
            azumith_in_xy = unpack.Basemap.angleTransform(trackx, tracky)
            # Compute slope in the direction of flow
            laserdf["flowslopes"] = flowdatabase.get_flow_slopes(laserdf["dh_fit_dx"], laserdf["dh_fit_dy"], laserdf["dh_fit_dx_sigma"], azumith_in_xy, laserdf["flow_angle"], laserdf["flow_angle_error"])
        
            # Remove None Values
            laserdf.dropna()
        
            if len(laserdf["flowslopes"]) >= fill*window:

                # only use every fill(th) point
                laserdf.iloc[::fill]

                # apply a rolling average
                laserdf["flowslopes"] = laserdf["flowslopes"].rolling(window).mean()

                # remove edges from trackx and tracky
                laserdf.dropna()

                # remove points outside of std deviation
                stddev = laserdf["flowslopes"].std()
                average = laserdf["flowslopes"].mean()
                laserdf[laserdf["flowslopes"] <= (average+2*stddev)]
                laserdf[laserdf["flowslopes"] >= (average-2*stddev)]
                allx, ally, allslope = np.hstack((allx, laserdf["x"])), np.hstack((ally, laserdf["y"])), np.hstack((allslope, laserdf["flowslopes"]))
            lasercount += 1
        
        i += 1
        if i >= 1:
            break

    print("preparing to plot")

    xi = np.arange(np.min(allx), np.max(allx), 1500)
    yi = np.arange(np.min(ally), np.max(ally), 1500)
    xi, yi = np.meshgrid(xi, yi)
    zi = griddata((allx, ally), allslope, (xi, yi), method="linear")

    basemap = unpack.Basemap("C:/Users/Nacnu/Documents/IceSat2/backgrounddata/ATA_adm0.shp")
    fig, ax = plt.subplots(1,1,sharex=True,sharey=True,figsize=(11,11))
    basemap.basemap_gpd.plot(ax=ax, color="white", edgecolor="black")
    plt.contourf(xi, yi, zi, levels = 15)
    groundline.linedata.plot(ax = ax, color = "crimson")
    plt.xlim(xlim[0], xlim[1])
    plt.ylim(ylim[0], ylim[1])
    plt.show()


plotTrackSlope()


