# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 12:55:11 2018

Testing horizon height calculation based on DEM
Samuli Launiainen 11.-12.6.2018
@author: slauniai
"""

import numpy as np
import matplotlib.pyplot as plt

# demfile = r'sve_1_dem_16m_aggr.asc'

def read_AsciiGrid(fname, setnans=True):
    
    """ reads AsciiGrid format in fixed format as below:
    
        ncols         750
        nrows         375
        xllcorner     350000: latitude in ykj
        yllcorner     6696000: longitude in ykj -kaista  3 edestä
        cellsize      16
        NODATA_value  -9999
        -9999 -9999 -9999 -9999 -9999
        -9999 4.694741 5.537514 4.551162
        -9999 4.759177 5.588773 4.767114
    IN:
        fname - filename (incl. path)
    OUT:
        data - 2D numpy array
        info - 6 first lines as list of strings
        (xloc,yloc) - lower left corner coordinates (tuple)
        cellsize - cellsize (in meters?)
        nodata - value of nodata in 'data'
    Samuli Launiainen Luke 7.9.2016
    
    NOTE! in reality at auras SVE-files: yllcorner = yulcorner?
    """
    import numpy as np

    fid = open(fname, 'r')
    info = fid.readlines()[0:6]
    fid.close()

    # print info
    # conversion to float is needed for non-integers read from file...
    ncols = float(info[0].split(' ')[-1])
    nrows = float(info[1].split(' ')[-1])    
    xloc = float(info[2].split(' ')[-1])
    yloc = float(info[3].split(' ')[-1])
    cellsize = float(info[4].split(' ')[-1])
    nodata = float(info[5].split(' ')[-1])

    # coordinate arrays    
    lon = xloc + cellsize * np.arange(ncols)
    lat = yloc + cellsize* (nrows - np.arange(nrows))

    # read rest to 2D numpy array
    data = np.loadtxt(fname, skiprows=6)

    if setnans is True:
        data[data == nodata] = np.NaN
        nodata = np.NaN

    return data, lat, lon, (xloc, yloc)


class Horizon(object):
    """
    data structure for horizon height calcualtions
    """
    def __init__(self, dem, dx = 1, lat=None, lon=None):
        """
        Args:
            dem - [m] elevation array n x m
            dx  - [m] grid size
            lat - latitude (not currently used)
            lon - longiture (not currently used)
        """
        self.Elev = dem # n x m matrix
        self.dx = dx # grid size (m)
        self.shape = np.shape(dem) # n rows, cols
        self.Lat = np.arange(0, self.shape[0]) # lat # row index, rows in n x m matrix
        self.Lon = np.arange(0, self.shape[1]) # lon # column index, columns n x m matrix        
        
    def calc_horizon(self, P0, R=500.0, Az_deg=[0.0], figs=False):
        """
        Computes horizon height (deg) for point P0 at direction Az_deg
        Args:
            P0 - Point - object
            R - viewing distance, float [m]
            Az_deg - Azimuth angles [deg, clockwise from North], array
            figs - True plots figures
        Returns:
            Az_deg - Azimuth angle [deg]
            horizon_angle - [deg relative to zero-plane]
            terrain_profile - dict with keys 'Az' and 'elev'
        """

        R = R / self.dx # in relative coordinates
        Az = np.deg2rad(Az_deg)
        
        # compute tangent z / r between P0 and all dem points
        x, y = np.meshgrid(self.Lon, self.Lat)
        r = np.sqrt(self.dx*((y - P0.lat)**2 + (x - P0.lon)**2))        
        z = (self.Elev - P0.elev)
        # tan_a = z / r
        elev_angle = np.rad2deg(np.arctan(z / r))

        target = np.array([P0.lat, P0.lon])
        dlat = -np.rint(R * np.cos(Az)) # change in lat, index
        dlon = np.rint(R * np.sin(Az)) # change in lon, index
        
        # outputs
        horizon_angle = np.zeros(np.shape(Az))* np.NaN
        terrain_profile = {'Az': Az_deg, 'elev': []}

        if figs:
            plt.figure(999)
            plt.imshow(self.Elev)
            plt.xlabel('Lon ix') # column index
            plt.ylabel('Lat ix') # row index
            cb = plt.colorbar()
            cb.set_label('Elev [m]')
            
        for k in range(len(horizon_angle)):
            # select grid cells self.R distance from P0
            source = np.array([[P0.lat + dlat[k], P0.lon + dlon[k]]])
            
            # call bresenham, return indices of cells ray travels through
            ix = bresenhamline(source, target)
            
            # convert to int and tuple, check that ix is within dem bounds
            ix = ix.astype(int)
            ix = (ix[:,0], ix[:,1])
            
            a = np.where((ix[0] >= 0) & (ix[0] <= self.shape[0]) & \
                         (ix[1] > 0) & (ix[1] < self.shape[1]))
            ix = (ix[0][a], ix[1][a])

            horizon_angle[k] = np.nanmax(elev_angle[ix])
            terrain_profile['elev'].append(self.Elev[ix][::-1])
            del a, ix
            
            if figs:
                plt.figure(999)
                plt.plot(source[0,1], source[0,0], 'ro', target[1], target[0], 'bs')
                plt.xlabel('lon id'); plt.ylabel('lat id')
                # plt.axis('square')
        
        if figs:
            plt.figure(888)
            plt.plot(Az_deg, horizon_angle, 'ro-')
            plt.xlabel('Az (deg)'); plt.ylabel('elev (deg)')
            
        return Az_deg, horizon_angle, terrain_profile
        
class Point(object):
    """
    defines point object; used as target
    """
    def __init__(self, lat, lon, elev):
        self.lat = lat # row index
        self.lon = lon # column index
        self.elev = elev # elevation value


#def tangent(X, Y, Z, p0):
#    # X, Y, Z - n x m matrixes
#    # p0 = (z0, y0, z0) point tuple
#    r = np.sqrt((X - p0[0])**2 + (Y - p0[1])**2)
#    z = (Z - p0[2])
#    
#    t = z / r
#    return t


"""
N-D Bresenham line algorithm (https://github.com/fjug/BobSeg/blob/master/bresenham.py)
"""
def _bresenhamline_nslope(slope):
    """
    Normalize slope for Bresenham's line algorithm.
    >>> s = np.array([[-2, -2, -2, 0]])
    >>> _bresenhamline_nslope(s)
    array([[-1., -1., -1.,  0.]])
    >>> s = np.array([[0, 0, 0, 0]])
    >>> _bresenhamline_nslope(s)
    array([[ 0.,  0.,  0.,  0.]])
    >>> s = np.array([[0, 0, 9, 0]])
    >>> _bresenhamline_nslope(s)
    array([[ 0.,  0.,  1.,  0.]])
    """
    scale = np.amax(np.abs(slope), axis=1).reshape(-1, 1)
    zeroslope = (scale == 0).all(1)
    scale[zeroslope] = np.ones(1)
    normalizedslope = np.array(slope, dtype=np.double) / scale
    normalizedslope[zeroslope] = np.zeros(slope[0].shape)
    return normalizedslope

def _bresenhamlines(start, end, max_iter):
    """
    Returns npts lines of length max_iter each. (npts x max_iter x dimension) 
    >>> s = np.array([[3, 1, 9, 0],[0, 0, 3, 0]])
    >>> _bresenhamlines(s, np.zeros(s.shape[1]), max_iter=-1)
    array([[[ 3,  1,  8,  0],
            [ 2,  1,  7,  0],
            [ 2,  1,  6,  0],
            [ 2,  1,  5,  0],
            [ 1,  0,  4,  0],
            [ 1,  0,  3,  0],
            [ 1,  0,  2,  0],
            [ 0,  0,  1,  0],
            [ 0,  0,  0,  0]],
    <BLANKLINE>
           [[ 0,  0,  2,  0],
            [ 0,  0,  1,  0],
            [ 0,  0,  0,  0],
            [ 0,  0, -1,  0],
            [ 0,  0, -2,  0],
            [ 0,  0, -3,  0],
            [ 0,  0, -4,  0],
            [ 0,  0, -5,  0],
            [ 0,  0, -6,  0]]])
    """
    if max_iter == -1:
        max_iter = np.amax(np.amax(np.abs(end - start), axis=1))
    npts, dim = start.shape
    nslope = _bresenhamline_nslope(end - start)

    # steps to iterate on
    stepseq = np.arange(1, max_iter + 1)
    stepmat = np.tile(stepseq, (dim, 1)).T

    # some hacks for broadcasting properly
    bline = start[:, np.newaxis, :] + nslope[:, np.newaxis, :] * stepmat

    # Approximate to nearest int
    return np.array(np.rint(bline), dtype=start.dtype)

def bresenhamline(start, end, max_iter=-1):
    """
    Returns a list of points from (start, end] by ray tracing a line b/w the
    points.
    Parameters:
        start: An array of start points (number of points x dimension)
        end:   An end points (1 x dimension)
            or An array of end point corresponding to each start point
                (number of points x dimension)
        max_iter: Max points to traverse. if -1, maximum number of required
                  points are traversed
    Returns:
        linevox (n x dimension) A cumulative array of all points traversed by
        all the lines so far.
    >>> s = np.array([[3, 1, 9, 0],[0, 0, 3, 0]])
    >>> bresenhamline(s, np.zeros(s.shape[1]), max_iter=-1)
    array([[ 3,  1,  8,  0],
           [ 2,  1,  7,  0],
           [ 2,  1,  6,  0],
           [ 2,  1,  5,  0],
           [ 1,  0,  4,  0],
           [ 1,  0,  3,  0],
           [ 1,  0,  2,  0],
           [ 0,  0,  1,  0],
           [ 0,  0,  0,  0],
           [ 0,  0,  2,  0],
           [ 0,  0,  1,  0],
           [ 0,  0,  0,  0],
           [ 0,  0, -1,  0],
           [ 0,  0, -2,  0],
           [ 0,  0, -3,  0],
           [ 0,  0, -4,  0],
           [ 0,  0, -5,  0],
           [ 0,  0, -6,  0]])
    """
    # Return the points as a single array
    return _bresenhamlines(start, end, max_iter).reshape(-1, start.shape[-1])