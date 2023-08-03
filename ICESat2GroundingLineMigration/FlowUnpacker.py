from doctest import ELLIPSIS_MARKER
from netCDF4 import Dataset as CDFData
import matplotlib.pyplot as plt
import numpy as np
import math
import warnings

class Database:
    """
    This class represents the NetCDF Antarctica flow data file.

    Parameters
    ----------
    filepath : str
          The path to the file containing the NetCDF flow data file

    Attributes
    ----------
    ERRX : numpy.array
          Ice flow error in X-direction
    ERRY : numpy.array
          Ice flow error in Y-direction
    lat : numpy.array
          Latitude for each datapoint
    lon : numpy.array
          Longitude for each datapoint
    VX : numpy.array
          Ice velocity in the X-direction
    VY : numpy.array
          Ice velocity in the Y-direction
    x : numpy.array
          X location of each datapoint
    y : numpy.array 
          Y location of each datapoint
    """
    
    def __init__(self, filepath):
        # use .filled to stop masked array and fill mask with nan values
        rootgroup = CDFData(filepath, "r", format="NETCDF4")
        self.ERRX, self.ERRY = rootgroup.variables["ERRX"][()].filled(np.nan), rootgroup.variables["ERRY"][()].filled(np.nan)
        self.lat, self.lon = rootgroup.variables["lat"][()].filled(np.nan), rootgroup.variables["lon"][()].filled(np.nan)
        self.VX, self.VY = rootgroup.variables["VX"][()].filled(np.nan), rootgroup.variables["VY"][()].filled(np.nan)
        self.x, self.y = rootgroup.variables["x"][()].filled(np.nan), rootgroup.variables["y"][()].filled(np.nan)

    def compute_angle(self):
        """Computes the angle of the ice velocity in radians from the X axis"""
        self.angle = np.arctan(self.VY / self.VX)

    def compute_speed(self):
        """Computes the speed of the ice velocity"""
        self.speed = np.sqrt(self.VX**2 + self.VY**2)

    def compute_error(self):
        """Computes the total error of the ice velocity"""
        self.error = np.sqrt(self.ERRX**2 + self.ERRY**2)

    def compute_angle_error(self):
        """Computes the total angle error of the ice velocity"""
        try:
            self.angle_error = self.error / (2*self.speed)
        except:
            raise Exception(f"Both error and speed need to be computed first to calculate angle error")

    def compute_all(self):
        """Computes angle of ice velocity, ice speed, total error, and total angle error"""
        self.compute_angle()
        self.compute_speed()
        self.compute_error()
        self.compute_angle_error()

    def rearrange_flow_data(self, xyindices):
        """
        Rearranges the flow data to the x and y indices

        Parameters
        ----------
        xyindices : numpy.array
              Numpy array containing X and Y indices. For each index, that index of the slope data will have the data at 
              the specified x and y indices of the flow data placed at it's index.

        Returns
        -------
        data : tuple
              (ERR, VEL, ANG, SPD, TOTERR, ANGERR)
              Rearranged flow data in order of slope data.
        """
        ERR = np.array([[self.ERRX[idx, idy], self.ERRY[idx, idy]] for idx, idy in xyindices])
        VEL = np.array([[self.VX[idx, idy], self.VY[idx, idy]] for idx, idy in xyindices])
        if self.angle:
            ANG = np.array([self.angle[idx, idy] for idx, idy in xyindices])
        if self.speed:
            SPD = np.array([self.speed[idx, idy] for idx, idy in xyindices])
        if self.error:
            TOTERR = np.array([self.error[idx, idy] for idx, idy in xyindices])
        if self.angle_error:
            ANGERR = np.array([self.angle_error[idx, idy] for idx, idy in xyindices])
        return ERR, VEL, ANG, SPD, TOTERR, ANGERR

    def rearrange_angle_data(self, xyindices):
        """
        Rearranges the angle data in the order of xy indices

        Parameters
        ----------
        xyindices : numpy.array
              Numpy array containing X and Y indices. For each index, that index of the slope data will have the data at 
              the specified x and y indices of the flow data placed at it's index.

        Returns 
        -------
        data : tuple
              (ANG, ANGERR)
              The rearranged angle and angle error data
        """
        try:
            ANG = np.array([self.angle[idx, idy] for idx, idy in xyindices])
            ANGERR = np.array([self.angle_error[idx, idy] for idx, idy in xyindices])
            return ANG, ANGERR
        except:
            raise Exception("Angle data for self does not exist")

    # known to work perfectly.
    @staticmethod
    def calc_along_flow_slope(dh_fit_dx, dh_fit_dy, dh_fit_dx_sigma, ascending, slope_azumith, flow_ang, flow_ang_err):
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

    def get_flow_slopes(self, dh_fit_dx, dh_fit_dy, dh_fit_dx_sigma, ascending, slope_azumith, flow_ang, flow_ang_err):
        """Performs calc_along_flow_slope on an array of data"""
        dh_fit_dx, dh_fit_dy, dh_fit_dx_sigma, flow_ang, flow_ang_err = np.array(dh_fit_dx), np.array(dh_fit_dy), np.array(dh_fit_dx_sigma), np.array(flow_ang), np.array(flow_ang_err)
        slope_azumith  = np.array(slope_azumith)
        along_flow_slope = np.array([Database.calc_along_flow_slope(dh_fit_dx[i], dh_fit_dy[i], dh_fit_dx_sigma[i], ascending[i], slope_azumith[i], flow_ang[i], flow_ang_err[i]) for i in range(len(flow_ang))])
        self.along_flow_slope = along_flow_slope
        return along_flow_slope

    def rolling_average(self, w):

        self.along_flow_slope_avg = None
        return self.along_flow_slope_avg

    def close(self):
        """Closes the dataset"""
        self.rootgroup.close()

    def plot_angles(self):
        image = plt.imshow(self.angle)
        plt.colorbar()
        plt.show()

