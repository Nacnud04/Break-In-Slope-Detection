import geopandas as gpd

class GLineData():

    """
    Takes in a shapefile for grounding line data and makes it usable and easy to manipulate

    Parameters
    ----------
    filepath : str
          Path to shapefile containing grounding line data
    bufferize : int
          Radius in meters around grounding line to generate buffer

    Variables
    ---------
    crs_proj4 : str
          Proj4 string of polar stereographic projection
    linedata : GeoDataFrame
          Geodata frame containing data transformed into polar projection
    linepolygon : GeoDataFrame
          Polygon with 100km radius around grounding line.
    """

    def __init__(self, filepath, buffersize):
        
        data = gpd.read_file(filepath)
        self.crs_proj4 = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
        self.linedata = data.to_crs(self.crs_proj4)
        
        self.linepolygon = self.generateBuffer(buffersize)

    def generateBuffer(self, buffersize):

        """
        Generates buffer around grounding line

        Parameters
        ----------
        buffersize : int
              Radius in m around grounding line to generate polygon

        Returns
        -------
        linepolygon : GeoDataFrame
              Buffer polygon around the grounding line
        """

        linepolygon = self.linedata.buffer(buffersize)
        linepolygon = gpd.GeoDataFrame(
              geometry=[linepolygon.unary_union], crs="+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs").explode(
              index_parts=False).reset_index(
              drop=True)
        return linepolygon