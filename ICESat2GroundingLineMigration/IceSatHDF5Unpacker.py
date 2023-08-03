import os
import h5py
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np
from cartopy import crs as ccrs
import math
from scipy.interpolate import griddata

import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)

colordict = {
        "gt1l":"lightcoral", "gt1r":"red",
        "gt2l":"lightblue", "gt2r":"blue",
        "gt3l":"lightgreen", "gt3r":"green"
        }

path = r'C:\Users\Nacnu\Documents\IceSat2\backgrounddata\ATA_adm0.shp'

def getColorDict():
    return colordict

class Basemap:
    def __init__(self, path:str):
        self.data = gpd.read_file(path)

        self.crs = ccrs.SouthPolarStereo()
        self.crs_proj4 = self.crs.proj4_init
        self.basemap_gpd = self.data.to_crs(self.crs_proj4)

    def getBasemap(self):
        return self.basemap_gpd

    def getCRS(self):
        return self.crs_proj4

def getColorDict():
    return colordict

class Dataset:

    """
    Passes an entire folder of h5 files and interprets it into granules and lasers.
    This is the highest level object in this program
    Allows for much ease and efficency with extracting hundreds of granules.

    Parameters
    ----------
    folder : str
          Path to folder containing dataset

    Variables
    ---------
    folder : str
          Identical to folder parameter
    files : list
          List of all files in folder
    datafiles : list
          List of all h5 files in folder.
    coverage : dict
          This identifies what files cover what segments and what regions. The keys for this dict
          are the region id's covered. Within each region is another dict containging keys which 
          are the RGT id's covered for each of the regions. In each RGT parameter in the dict there
          is a list of filenames which cover that RGT.
    """

    def __init__(self, folder:str):

        self.folder = folder

        files = os.listdir(folder)
        self.files = [f for f in files if os.path.isfile(folder+'/'+f)]

        self.datafiles = [f for f in files if f.split('.')[-1] == "h5"]

        ## group similar/overlapping files
        #self.coverage = {}
        #for datafile in self.datafiles:
        #    granule = Granule(self.folder, datafile)
        #    rgts, (start_region, end_region) = granule.getGranulePath()
        #    for region in range(start_region, end_region+1):
        #        for rgt in [x for i, x in enumerate(rgts) if i == 0 or x != rgts[1]]:
        #            if region not in self.coverage:
        #                self.coverage[region] = {rgt:[granule.filename]}
        #            else:
        #                if rgt not in self.coverage[region]:
        #                    self.coverage[region][rgt] = [granule.filename]
        #                if granule.filename not in self.coverage[region][rgt]:
        #                    self.coverage[region][rgt].append(granule.filename)
        #    del granule

    def coverageReport(self):
        """
        Returns a sorted list of all the regions and RGTs coverd by the various segments in the dataset.
        The number next to each region/RGT number is the count of files which cover that area.
        """
        regioncount = {}
        rgtcounts = {}
        regions = self.coverage.keys()
        for region in regions:
            regioncount[region] = sum([len(self.coverage[region][key]) for key in self.coverage[region].keys()])
            rgtcounts = {**rgtcounts, **{key:len(self.coverage[region][key]) for key in self.coverage[region].keys()}}
        print("""======================\nCOVERAGE REPORT\n=====================""")
        print("Region : Count")
        for key in sorted(regioncount.keys()):
            print(f"{key} : {regioncount[key]}")
        print("RGT : Count")
        for key in sorted(rgtcounts.keys()):
            print(f"{key} : {rgtcounts[key]}")

    def openspecific(self, index):
        """
        Opens and returns a specific granule via index.

        Parameters
        ----------
        index : int
              Index of file in list datafiles

        Returns
        -------
        granule : class Granule
        """
        granule = Granule(self.folder, self.datafiles[index])

        return granule

    def openfilename(self, filename):
        """
        Opens and returns a specific granule via filename

        Parameters
        ----------
        filename : str
              Name of the file

        Returns
        -------
        granule : class Granule
        """
        granule = Granule(self.folder, filename)

        return granule
    
    def opens3(self, link, s3):
        granule = Granule(self.folder, link, s3=s3, linked=True)
        return granule

    # analyzes self.coverage 
    def returnFilenamesByRegion(self, region):
        """
        Takes in a region id and returns a list of filenames in that region

        Parameters
        ----------
        region : int
              Id of region desired

        Returns
        -------
        filenames : list
              List of filenames in desired region
        """
        total = self.coverage[region]
        filenames = tuple([item for sublist in [total[key] for key in total.keys()] for item in sublist])
        return filenames
    
    def returnFilenamesByRGT(self, rgt):
        """
        Takes in a desired RGT and returns the filenames which cover that RGT

        Parameters
        ----------
        rgt : int
              Id of desired rgt

        Returns
        -------
        filenames : list
              List of filenames in desired rgt
        """
        filenames = []
        keyskeys = {highlevelkey:self.coverage[highlevelkey].keys() for highlevelkey in self.coverage.keys()}
        for highkey in list(keyskeys.keys()):
            for lowkey in keyskeys[highkey]:
                if int(lowkey) == int(rgt):
                    filenames.extend(self.coverage[highkey][lowkey])
        return filenames

    def returnFilenamesByPolygon(self, polygons):
        """
        Take a polygon and return all filenames which have data within the desired polygon

        Parameters
        ----------
        polygon : shapely.polygon
              Polygon for which to extract data from

        Returns
        -------
        filenames : list
              List of filenames within the polygon
        """

        #pd.concat([polygon_series.contains(point) for point in point_series], axis=1)

        files = []

        crs_proj4 = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

        for filename in self.datafiles:
            granule = self.openfilename(filename)
            granpoly = granule.captureCoverage().to_crs(crs_proj4)[0]
            for polygon in polygons.geometry:
                if polygon.intersects(granpoly) == True:
                    files.append(filename)
        return files

    def returnTrackDataByRegion(self, region, orientation=None, name=None, num=None, strength=None):
        """
        Returns file track data if in desired region
        Input parameters can be used to filter data

        Parameters
        ----------
        region : int
              Desired region to pull data from
        orientation : string
              Denotes direction of satellite (** NEEDS TO BE FIXED **)
              Can either be forwards or backwards
        name : str
              Laser name
        num : int
              Id of laser number
        strength : bool
              Strength of laser. If set to True the laser is strong.

        Returns
        -------
        alltrackdata : list
              This is a list of tuples containing the output track data
        """
        filenames = self.returnFilenamesByRegion(region)
        alltrackdata = []
        for filename in filenames:
            granule = Granule(self.folder, filename)
            trackdata = granule.getTrackData()
            for lasertrack in trackdata:
                (lat, lon, time, dh_fit_dx, dh_fit_dx_sigma, metadata, granuledata) = lasertrack
                if orientation and orientation == granuledata["orientation"]:
                    alltrackdata.append(lasertrack)
                elif name and name == metadata["name"]:
                    alltrackdata.append(lasertrack)
                elif num and num == metadata["number"]:
                    alltrackdata.append(lasertrack)
                elif strength and strength == metadata["strength"]:
                    alltrackdata.append(lasertrack)
                elif not orientation and not name and not num and not strength:
                    alltrackdata.append(lasertrack)

        return alltrackdata

    def returnTrackDataByRGT(self, RGT, orientation=None, name=None, num=None, strength=None):
        """
        Returns file track data if in desired RGT
        Input parameters can be used to filter data

        Parameters
        ----------
        RGT : int
              Desired RGT to pull data from
        orientation : string
              Denotes direction of satellite (** NEEDS TO BE FIXED **)
              Can either be forwards or backwards
        name : str
              Laser name
        num : int
              Id of laser number
        strength : bool
              Strength of laser. If set to True the laser is strong.

        Returns
        -------
        alltrackdata : list
              This is a list of tuples containing the output track data
        """
        filenames = self.returnFilenamesByRGT(RGT)
        alltrackdata = []
        for filename in filenames:
            granule = Granule(self.folder, filename)
            trackdata = granule.getTrackData()
            for lasertrack in trackdata:
                (lat, lon, time, dh_fit_dx, dh_fit_dx_sigma, metadata, granuledata) = lasertrack
                if orientation and orientation == granuledata["orientation"]:
                    alltrackdata.append(lasertrack)
                elif name and name == metadata["name"]:
                    alltrackdata.append(lasertrack)
                elif num and num == metadata["number"]:
                    alltrackdata.append(lasertrack)
                elif strength and strength == metadata["strength"]:
                    alltrackdata.append(lasertrack)
                elif not orientation and not name and not num and not strength:
                    alltrackdata.append(lasertrack)

        return alltrackdata

    def returnTrackDataByPolygon(self, polygon):
        """
        Returns track data for all tracks within specified polygon.

        Parameters
        ----------
        polygon : class shapely.Polygon
              Polygon for which to extract data from

        Returns
        -------
        trackdata : list
              List of trackdata tuples within the polygon
        """
        alldata = ()
        for filename in self.datafiles:
            granule = self.openfilename(filename)
            granpoly = granule.captureCoverage()
            if polygon.overlaps(granpoly).any() == True:
                alldata = alldata + granule.getTrackData()
        return alldata

    @staticmethod
    def showBackgroundMap():
        """Displays entire background map with no data"""

        basemap = Basemap(path)
        basemap.getBasemap.plot()

        plt.show()

    def showTracks(self, local=False):
        """Shows all tracks on basemap for the passed dataset"""

        alltrackdata = ()

        for i in range(len(self.datafiles)):

            granule = self.openspecific(i)
            print(f"Getting tracks from: {granule.filename}")
            alltrackdata = alltrackdata + granule.getTrackData()
            granule.close()

        basemap = Basemap(path)
        basemap.plotTracks(alltrackdata, datatype=tuple, local=local)

class Granule:

    """
    Object which parses h5py file. Allows for ease of data extraction.

    Parameters
    ----------
    folder : str
          Path to folder containing the granule
    filename : str
          Filename inside of folder containing granule
    
    Variables
    ---------
    folder : str
          Identical to folder parameter
    filename : str
          Identical to filename parameter
    start_rgt | end_rgt : int
          Index of start and end region ground track
    start_region | end_region : int
          Index of start and end regions
    start_cycle | end_cycle : int
          Number of cycles at the start and end of granule
    metadata : dict
          Compaction of previous 6 variables for ease of passing between classes.
          Keys are the names of previous 6 variables
    laser** : class Laser
          6 class Laser variables each labed by laser name (e.g. laser1l)
    lasers : tuple
          contains all Laser objects for easier interpretation
    """

    def __init__(self, folder:str, filename:str, s3=None, linked = False):

        self.folder = folder
        self.filename = filename
        
        if linked == False:
            filedata = h5py.File(f'{self.folder}/{self.filename}','r')
        elif linked == True:
            filedata = h5py.File(s3.open(filename,'rb'),'r')

        # capture data for matching of granules
        self.start_rgt, self.end_rgt = filedata["ancillary_data"]["start_rgt"][0], filedata["ancillary_data"]["end_rgt"][0]
        self.start_region, self.end_region = filedata["ancillary_data"]["start_region"][0], filedata["ancillary_data"]["end_region"][0]
        self.start_cycle, self.end_cycle = filedata["ancillary_data"]["start_cycle"][0], filedata["ancillary_data"]["end_cycle"][0]

        if self.start_region > 4 and self.start_region < 11:
            orientation = "Backward"
        elif self.start_region < 4 and self.start_region > 11:
            orientation = "Forward"
        else:
            orientation = "Forward"

        self.metadata = {"folder":self.folder,"filename":self.filename,"orientation":orientation,"start_rgt":self.start_rgt,"end_rgt":self.end_rgt,
                         "start_region":self.start_region,"end_region":self.end_region,"start_cycle":self.start_cycle,"end_cycle":self.end_cycle}

        self.laser1l, self.laser1r = Laser(filedata["gt1l"],"gt1l", filename, self.metadata), Laser(filedata["gt1r"],"gt1r", filename, self.metadata)
        self.laser2l, self.laser2r = Laser(filedata["gt2l"],"gt2l", filename, self.metadata), Laser(filedata["gt2r"],"gt2r", filename, self.metadata)
        self.laser3l, self.laser3r = Laser(filedata["gt3l"],"gt3l", filename, self.metadata), Laser(filedata["gt3r"],"gt3r", filename, self.metadata)

        self.lasers = (self.laser1l, self.laser1r, self.laser2l, self.laser2r, self.laser3l, self.laser3r)
        
        self.sc_orient = filedata["orbit_info"]["sc_orient"]

    def getGranulePath(self):
        """
        Returns granule location in terms of rgt and region.
        Used for finding overlap between granules more easily.

        Returns
        -------
        pathparams : list
            This is a list where the 1st term is a list of length 2 containing start and end rgt data.
            The second term is a tuple of length 2 containing start and end region data.
        """
        return [[self.start_rgt, self.end_rgt],(self.start_region, self.end_region)]

    def getLasers(self):
        """Returns a list of laser classes for each laser in the granule."""
        return self.laser1l, self.laser1r, self.laser2l, self.laser2r, self.laser3l, self.laser3r

    def getTrackData(self):
        # returns a tuple with track data for each laser (each laser is a tuple index)
        trackdata = ()
        for laser in self.lasers:
            trackdata = trackdata + (laser.getTrackData(),)
        return trackdata

    def getAzimuthData(self):
        """
        Gathers azumith data from all lasers in track
        
        Returns
        -------
        azumiths : tuple
            Is a tuple containing track data
        """
        azumithdata = ()
        for laser in self.lasers:
            azumithdata = azumithdata + (laser.returnAzumith(),)
        return azumithdata

    def showTrack(self, local=False):

        """
        Plots all tracks in the granule

        Parameters
        ----------
        local : bool
              Designates local or entire basemap visualization of the granule tracks.
        """

        trackmap = Basemap(path)
        trackmap.plotTracks(self.lasers, datatype=object, local=local)

    def captureCoverage(self):
        """
        Finds maximum and minimum lat and lon of granule
        
        Returns
        -------
        coveredpolygon : class shapely.Polygon
              Rectangular shapely polygon which all data is covered by
        """
        alllats, alllons = np.array([]), np.array([])
        for laser in self.lasers:
            data = laser.getTrackData()
            lat, lon = data[0], data[1]
            alllats = np.append(alllats, lat)
            alllons = np.append(alllons, lon)
        maxlat, minlat = np.max(alllats), np.min(alllats)
        maxlon, minlon = np.max(alllons), np.min(alllons)
        coveredpolygon = Polygon([(minlon,minlat),(maxlon,minlat),(maxlon,maxlat),(minlon,maxlat),(minlon,minlat)])
        return gpd.GeoSeries([coveredpolygon], crs="EPSG:4326")

    def close(self):
        """Closes the h5 file"""
        del self.lasers

    def returnDates(self):
        """Returns date of granule in start and end times. Formatted as utc"""
        return (self.start_utc, self.end_utc)

class Laser:
    """
    The laser class stores and interprets data for each laser aboard ICESat-2
    This only manages one laser from one granule at a time

    Parameters
    ----------
    group : class (h5py_hl.group.Group)
          The data itself. Aka that lasers group data in the h5 file.
    name : str
          The name of the laser as given by its group name. (can be either gt1l, gt1r, gt2l, etc.)
    filename : str
          The name of the file for which the laser data came from.
    granuledata : dict
          Granule data about the lasers parent granule.
          Contains attributes: folder, filename, orientation, start_rgt, end_rgt, start_region, 
                    end_region, start_cycle, end_cycle.

    Attributes
    ----------
    name : str
          Identical to name parameter
    data : class (h5py_gl.group.Group)
          Identical to group parameter
    parentfile : str
          Identical to filename parameter
    granuledata : dict
          Identical to granuledata parameter
    land_ice_segments : dict
          Contains all land_ice_segments group data. Is a dictionary of dictionaries
    segment_quality : dict
          Contains all segment_quality group data. Is a dictionary of dictionaries
    residual_histogram : dict
          Contains all residual_histogram group data. Is a dictionary of dictionaries
    beamnum : int
          Laser beam number. Estimated via orientation parameter from granule data
    strength : bool
          Strength of laser. True = Strong
    metadata : dict
          Laser metadata. Has following keys: name, number, strength, orientation
    """

    gt2beamnumBackwards = {"gt1l":1,"gt1r":2,"gt2l":3,"gt2r":4,"gt3l":5,"gt3r":6}
    gt2beamnumForwards = {"gt1l":6,"gt1r":5,"gt2l":4,"gt2r":3,"gt3l":2,"gt3r":1}
    strengthDict = {1:True, 2:False, 3:True, 4:False, 5:True, 6:False}

    def __init__(self, group, name:str, filename:str, granuledata:dict):

        self.name = name
        self.data = group
        self.parentfile = filename
        self.granuledata = granuledata

        orientation = granuledata["orientation"]

        try:
            self.land_ice_segments = self.data["land_ice_segments"]
            fit_statistics = self.land_ice_segments["fit_statistics"]
            self.dh_fit_dx = fit_statistics["dh_fit_dx"][()] #Along tack slope from along track segment fit
            self.dh_fit_dx_sigma = fit_statistics["dh_fit_dx_sigma"][()] #Propagated error in the along-track segment slope
        except KeyError:
            raise Exception(f"land_ice_segments in file: {filename} not found")

        #try :
        #    self.segment_quality = self.data["segment_quality"]
        #except KeyError:
        #    warnings.warn(f"segment_quality in file: {filename} not found")

        #try :
        #    self.residual_histogram = self.data["residual_histogram"]
        #except KeyError:
        #    warnings.warn(f"residual_histogram in file: {filename} not found")

        if orientation == "Forward":
            self.beamnum = Laser.gt2beamnumForwards[name]
        elif orientation == "Backward":
            self.beamnum = Laser.gt2beamnumBackwards[name]
        else:
            raise Exception("orientation given which was neither forwards nor backwards")

        self.strength = Laser.strengthDict[self.beamnum]

        self.metadata = {"name":self.name, "number":self.beamnum, "strength":self.strength, "orientation":orientation}
        # metadata is a useful way of compacting small bits of data into a single variable
        # metadata has name, number, strength, and orientation

    def showTrack(self, local=False):

        """
        Shows track data on basemap

        Parameters
        ----------
        local : bool
              If set True only area local to track will be plotted (rest of basemap is cropped out).
              Otherwise entire basemap shows.
        """

        lat = self.land_ice_segments["latitude"][()]
        lon = self.land_ice_segments["longitude"][()]
        time = self.land_ice_segments["delta_time"][()]
        h_li = self.land_ice_segments["h_li"][()]

        trackmap = Basemap(path)
        trackmap.plotTrack(lat, lon, time, h_li, self.metadata, local=local)

    def getTrackData(self):

        """
        Returns laser track data

        Returns
        -------
        trackdata : tuple
            Data about laser track. Formatted as (lat, lon, time, dh_fit_dx, dh_fit_dx_sigma, metadata, granuledata)
        """

        lat = self.land_ice_segments["latitude"][()]
        lon = self.land_ice_segments["longitude"][()]
        time = self.land_ice_segments["delta_time"][()]
        h_li = self.land_ice_segments["h_li"][()]
        return (lat, lon, time, h_li, self.dh_fit_dx, self.dh_fit_dx_sigma, self.metadata, self.granuledata)

    def returnAzimuth(self):
        """
        Returns laser track azumith data.

        Returns
        -------
        azumith : np.array
            Array of azumith points in radians measured eastward from north.
        """
        azumith = self.land_ice_segments["ground_track"]["ref_azimuth"][()]
        return azumith

    def returnAcrossTrackSlope(self):
        """
        Returns across track slope for each datapoint

        Returns
        -------
        dh_fit_dy : np.array
              Array of across track slope data
        """
        dh_fit_dy = self.land_ice_segments["fit_statistics"]["dh_fit_dy"][()]
        return dh_fit_dy

class Basemap:

    """
    Takes in data to be visualized and outputs it onto a map of the local area.

    Parameters
    ----------
    path : str
          Path to shapefile containing basemap shapefile

    Variables
    ---------
    data : class GeoDatabase
          Geodatabase for the basemap
    crs : object
          Cartopy coordinate reference system for projection (south polar stereographic)
    crs_proj4 : object
          Cartopy projection of crs
    """

    def __init__(self, path:str):
        self.data = gpd.read_file(path)

        self.crs = ccrs.SouthPolarStereo()
        self.crs_proj4 = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
        self.basemap_gpd = self.data.to_crs(self.crs_proj4)

    def plotTrack(self, lat, lon, time, dh_fit_dx, dh_fit_dx_sigma, metadata, local = False):

        """
        Plots a single track of data on the basemap

        Parameters
        ----------
        lat : list
              list of latitude points
        lon : list
              list of longitude points
        time : list
              list of timestamps for each datapoint (seconds since reference epoch)
        dh_fit_dx : list
              slope at each datapoint
        dh_fit_dx_sigma : list
              propogated error of the slope at each datapoint
        metadata : dict
              metadata for input track. Has the following parameters: name, number, strength, orientation
        """

        track1l_df = pd.DataFrame({'time': time, 'slope': dh_fit_dx, 'slope_error': dh_fit_dx_sigma, 'lat': lat, 'lon':lon})
        track1l_gdf = gpd.GeoDataFrame(track1l_df, geometry=gpd.points_from_xy(track1l_df.lon, track1l_df.lat), crs="EPSG:4326")

        if local == True:
            maxlat, minlat = max(lat)*1.25, min(lat)/1.25
            maxlon, minlon = max(lon)*1.25, min(lon)/1.25

            mask = Polygon([(minlon,minlat),(maxlon,minlat),(maxlon,maxlat),(minlon,maxlat),(minlon,minlat)])

            track1l_gdf = gpd.clip(track1l_gdf, gpd.GeoSeries([mask], crs="EPSG:4326"))
            basemap_gpd_local = gpd.clip(self.data, gpd.GeoSeries([mask]))

            track1l_gdf = track1l_gdf.to_crs(self.crs_proj4)
            basemap_gpd_local = basemap_gpd_local.to_crs(self.crs_proj4)

            ax = basemap_gpd_local.plot(color="white", edgecolor="black")

            track1l_gdf.plot(ax = ax, color=colordict[metadata["name"]], markersize = 4)

        else:
            ax = self.basemap_gpd.plot(color="white", edgecolor="black")

            track1l_gdf_polar = track1l_gdf.to_crs(self.crs_proj4)

            track1l_gdf_polar.plot(ax = ax, color=colordict[metadata["name"]], markersize = 4)

        plt.show()

    def plotTracks(self, tracks, datatype, local = False):

        """
        Plots all inputted tracks on the basemap
        
        Parameters
        ----------
        tracks : list
              This is a list of either track data as a tuple or laser objects.
        datatype : type
              Can either be the type object or tuple. This specifies the type of data input in tracks.
        local : bool
              Specifies whether to plot locally or not. If set to true it only plots the basemap around the input data.
        """

        alllats, alllons = np.array([]), np.array([])
        
        if local:
            for track in tracks:

                if datatype == object:
                    lat, lon, time, dh_fit_dx, dh_fit_dx_sigma, metadata, granuledata = track.getTrackData() # each item in lasertracks is another tuple which goes as follows (lat, lon, time, height, name)
                elif datatype == tuple:
                    lat, lon, time, dh_fit_dx, dh_fit_dx_sigma, metadata, granuledata = track
                alllats = np.append(alllats, lat)
                alllons = np.append(alllons, lon)

            maxlat, minlat = np.max(alllats) * 1.25, np.min(alllats) / 1.25
            maxlon, minlon = np.max(alllons) * 1.25, np.min(alllons) / 1.25

            mask = Polygon([(minlon,minlat),(maxlon,minlat),(maxlon,maxlat),(minlon,maxlat),(minlon,minlat)])

        if local:
            basemap_gpd = gpd.clip(self.data, gpd.GeoSeries([mask], crs="EPSG:4326"))
        else:
            basemap_gpd = self.data

        basemap_gpd = basemap_gpd.to_crs(self.crs_proj4)
        ax = basemap_gpd.plot(color="white", edgecolor="black")

        for track in tracks:

            if datatype == object:
                lat, lon, time, dh_fit_dx, dh_fit_dx_sigma, metadata, granuledata = track.getTrackData() # each item in lasertracks is another tuple which goes as follows (lat, lon, time, height, name)
            elif datatype == tuple:
                lat, lon, time, dh_fit_dx, dh_fit_dx_sigma, metadata, granuledata = track

            name = metadata["name"]

            track_df = pd.DataFrame(
                {'time': time, 'slope': dh_fit_dx, 'slope_err': dh_fit_dx_sigma, 'lat': lat, 'lon':lon})
            track_gdf = gpd.GeoDataFrame(track_df, geometry=gpd.points_from_xy(track_df.lon, track_df.lat), crs="EPSG:4326")

            if local:
                track_gdf = gpd.clip(track_gdf, gpd.GeoSeries([mask], crs="EPSG:4326"))

            track_gdf = track_gdf.to_crs(self.crs_proj4)

            track_gdf.plot(ax = ax, color=colordict[metadata["name"]], markersize=4)

        plt.show()

    def getBasemap(self):
        """Returns the basemap as a geodatabase"""
        return self.basemap_gpd

    def getCRS(self):
        """Returns the cartopy projection of the crs"""
        return self.crs_proj4

    @staticmethod
    def transformProjection(geodataframe):
        """
        Takes in a geodata frame and transforms the data into south polar sterographic

        Reuturns
        --------
        gpd : GeoDataFrame
              Input geodataframe converted into the south polar stereographic coordinate system
        """
        crs = ccrs.SouthPolarStereo()
        crs_proj4 = crs.proj4_init
        return geodataframe.to_crs(crs_proj4)

    @staticmethod 
    def find_nearest(array, value, side):
        """
        Takes in an array and does a binary search to find the closest value to the desired value.

        Parameters
        ----------
        array : numpy array
              Sorted 1D input array to search through.
        value : float/int
              Desired value to search for
        side : str
              Side of 1D array with the smallest values. Can be either "left" or "right"

        Returns
        -------
        idx : int
              Index of closest value in array, crs
        """
        idx = np.searchsorted(array, value, side=side)
        if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
            return idx
        else:
            return idx

    @staticmethod 
    def find_nearest_grid(arrx, arry, valx, valy):
        """
        Takes in a point and finds its nearest point on an xy grid.

        Parameters
        ----------
        arrx : numpy.array
              Increasing list of x values in grid
        arry : numpy.array
              Decreasing list of y values in grid
        valx : int/float
              X location
        valy : int/float
              Y kocation

        Returns
        -------
        id : tuple
              Contains idx, idy. Is the index of the grid for the nearest point to the input x, y.
        """

        idx = Basemap.find_nearest(arrx, valx, side = "left")
        idy = len(arry) - Basemap.find_nearest(np.flip(arry), valy, side = "left")

        return idx, idy

    @staticmethod 
    def latlonindex_to_grid(lat, lon, grid):
        """
        Takes in a list of longitude and latitude and creates an index of xy points nearest to each lat lon for a given input grid.

        Parameters
        ----------
        lat : numpy.array
              Array of latitude points
        lon : numpy.array
              Array of longitude points
        grid : FlowUnpacker.Dataset
              Object representing the flow data containing the grid information needed

        Returns
        -------
        xyindices : list
              List of indexes for the x and y values of the flow data to assign to the slope data
        xvals : list
              List of x positions of the slope data
        yvals : list
              List of y positions of the slope data
        """
        crs_proj4 = ccrs.SouthPolarStereo().proj4_init
        crs_proj4 = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
        track_df = pd.DataFrame({'lat': lat, 'lon':lon})
        track_gdf = gpd.GeoDataFrame(track_df, geometry=gpd.points_from_xy(track_df.lon, track_df.lat), crs="EPSG:4326")
        track_gdf = track_gdf.to_crs(crs_proj4)
        track = np.array([(point.x, point.y) for point in track_gdf.geometry])

        xyindices = track_gdf.apply(lambda row: Basemap.find_nearest_grid(grid.x, grid.y, row.geometry.x, row.geometry.y), axis=1)
        
        return xyindices, track[:,0], track[:,1]

    @staticmethod
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

    @staticmethod 
    def angleTransform(xvals, yvals):
        """Performs angNorth_to_angXY on an array of data"""
        angles = []
        for i in range(len(xvals)):
            if i >= 0 and i < len(xvals) - 1:
                angles.append(Basemap.XYvec_to_ang((xvals[i], yvals[i]),(xvals[i+1], yvals[i+1])))
            elif i == len(xvals) - 1:
                angles.append(Basemap.XYvec_to_ang((xvals[i-1], yvals[i-1]),(xvals[i], yvals[i])))
        angles = np.array(angles)
        return angles
