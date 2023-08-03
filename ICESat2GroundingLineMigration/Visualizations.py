# Visualizations file
# Used to visualize data
from math import sin, cos, pi
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import geopandas as gpd
import numpy as np

class Angles:

    """
    Allows for trivial visualizations of parameters related to angles
    """

    def __init__(self):
        pass

    @staticmethod
    def visualize_ang_xy(basemap, xs, ys, azumiths, vectorlength:float, fill:int):
        """
        Visualize the azumith of the track at regular point intervals

        Parameters
        ----------
        basemap : GeoDataFrame
              Geodataframe of antarctica basemap.
        xs : List
              List of projected x values of track in polar stereographic projection
        ys : List
              List of projected y values of track in polar stereographic projection
        azumiths : List
              Azumiths in radians from the x axis going counter clockwise.
        vectorlength : Float
              Length of vector to generate which shows azumith direction
        fill : Int
              Number designating how many points should be plotted.
              The number plotted is the inverse of this value.
              So 100 is 1% of points
        """
        crs_proj4 = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

        basemap_gpd = basemap.to_crs(crs_proj4)

        maxx, minx = max(xs)*1.25, min(xs)/1.25
        maxy, miny = max(ys)*1.25, min(ys)/1.25

        ax = basemap_gpd.plot(color="white", edgecolor="black")
        xs = [x for i, x in enumerate(xs) if i % fill == 0]
        ys = [y for i, y in enumerate(ys) if i % fill == 0]
        azumiths = [a for i, a in enumerate(azumiths) if i % fill == 0]
        ax.scatter(xs, ys)
        for i in range(len(xs)):
            ax.plot([xs[i], xs[i]+vectorlength*cos(azumiths[i])], [ys[i], ys[i]+vectorlength*sin(azumiths[i])], color = "red")
        plt.xlim(minx, maxx)
        plt.ylim(miny, maxy)
        plt.show()

    @staticmethod 
    def flowslope(dh_fit_dx, dh_fit_dy, azumith, flowdirec, flowslope):
        """
        Plots the vectors used to calculate the slope in the direction of flow in 3d

        Parameters
        ----------
        dh_fit_dx : Float
              Slope in the track direction
        dh_fit_dy : Float
              Slope in the across track direction
        azumith : Float
              Radians of slope direction from x axis counterclockwise
        flowdirec : Float
              Radians of flow direction from x axis counterclockwise
        flowslope : Float
              Slope in the direction of flow
        """
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        trackslope = np.array([cos(azumith), sin(azumith), dh_fit_dx])
        ts_hat = trackslope / np.linalg.norm(trackslope)
        acrosstrackslope = np.array([sin(azumith), -1*cos(azumith), dh_fit_dy])
        ats_hat = acrosstrackslope / np.linalg.norm(acrosstrackslope)
        flowslope = np.array([cos(flowdirec), sin(flowdirec), flowslope])
        f_hat = flowslope / np.linalg.norm(flowslope)
        ax.plot([0, cos(azumith)], [0, sin(azumith)], [0, 0], color="red")
        ax.plot([0, ts_hat[0]], [0, ts_hat[1]], [0, ts_hat[2]], color="red")
        ax.plot([0, ats_hat[0]], [0, ats_hat[1]], [0, ats_hat[2]], color="green")
        ax.plot([0, cos(flowdirec)], [0, sin(flowdirec)], [0, 0], color="blue")
        ax.plot([0, f_hat[0]], [0, f_hat[1]], [0, f_hat[2]], color="blue")
        plt.show()


class Polygons():

    """
    Visualize polygon data
    """

    def __init__(self):
        pass

    @staticmethod 
    def line_buffer(basemap, line, buff):
        """
        Allows for visualization of the grounding line with a set buffer around it

        Parameters
        ----------
        basemap : GeoDataFrame
              Basemap data of antarctica
        line : GeoSeries
              Grounding line geodataframe
        buff : GeoSeries
              Grounding line buffer geodataframe
        """
        crs_proj4 = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
        basemap_gpd = basemap.to_crs(crs_proj4)
        fig, ax = plt.subplots(1,1,sharex=True,sharey=True,figsize=(11,11))
        basemap_gpd.plot(ax = ax, color="white", edgecolor="black")
        buff.plot(ax = ax, color = "lightblue", edgecolor = "steelblue", alpha = 0.5)
        line.plot(ax = ax, color = "crimson")
        plt.show()