# Takes in a track from laser 1l of a given file and visualizes the azumith
# Additionally this visualizes the slope vectors

import IceSatHDF5Unpacker as unpack
import FlowUnpacker as flow
import Visualizations as visualize
import geopandas as gpd

# Path to antarctica basemap
path = r'C:\Users\Nacnu\Documents\IceSat2\backgrounddata\ATA_adm0.shp'

# Import the ice flow database
flowdatabase = flow.Database("data/antarctic_ice_vel_phase_map_v01.nc")

# Compute all angles and errors from the flow database
flowdatabase.compute_all()

# Unpack the icesat2 dataset
dataset = unpack.Dataset(r"C:\Users\Nacnu\Documents\IceSat2")

# Open a specific granule in the dataset
granule = dataset.openfilename("ATL06_20181114050603_07110110_005_01.h5")

# Get the data from track 1l in the granule
lat, lon, time, dh_fit_dx, dh_fit_dx_sigma, metadata, granuledata = granule.laser1l.getTrackData()

# Return across track slope
dh_fit_dy = granule.laser1l.returnAcrossTrackSlope()

# Transform the track lat lon values to the correct projected x and y values and return xyincices
xyindices, trackx, tracky = unpack.Basemap.latlonindex_to_grid(lat, lon, flowdatabase)

# Rearrange flow angle and flow angle error in the order of the xy indices
flow_angle, flow_angle_error = flowdatabase.rearrange_angle_data(xyindices)

# Get the projected azumith at each x and y point
azumith_in_xy = unpack.Basemap.angleTransform(trackx, tracky)

# Visualize the azumiths
visualize.Angles.visualize_ang_xy(gpd.read_file(path), trackx, tracky, azumith_in_xy, 1000, 100) 

# Compute slope in the direction of flow
flowslopes = flowdatabase.get_flow_slopes(dh_fit_dx, dh_fit_dy, dh_fit_dx_sigma, azumith_in_xy, flow_angle, flow_angle_error)

# Visualize vectors at index 2000.
visualize.Angles.flowslope(dh_fit_dx[2000], dh_fit_dy[2000], azumith_in_xy[2000], flow_angle[2000], flowslopes[2000])

