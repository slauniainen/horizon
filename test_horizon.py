# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 16:56:36 2018

@author: slauniai
"""

import numpy as np
import matplotlib.pyplot as plt
import dem_horizon as h

# import dem example
demfile = r'sve_1_dem_16m_aggr.asc'

dem, lat, lon, _ = h.read_AsciiGrid(demfile)

#plt.figure(1)
#plt.imshow(dem)
#plt.xlabel('Lon ix') # column index
#plt.ylabel('Lat ix') # row index
#cb = plt.colorbar()
#cb.set_label('Elev [m]')


# set parameters
dx = 16.0 # dem resolution (m)
view_dist = 500.0 # viewing distance
sectors = 36.0 # sectors = 360 (deg) / azimuth interval (deg)

# create Horizon object
H = h.Horizon(dem, dx=16)
print('rows, cols', H.shape)

# let's compute horizon for point Lat_ix = 40, Lon_ix 50
lat_ix = 140
lon_ix = 60
P0 = h.Point(lat_ix, lon_ix, dem[lat_ix, lon_ix])

# calculate horizon height angle (deg)
Az = np.linspace(0, 360, sectors)

Az, elev_angle, terrain_profile = H.calc_horizon(P0, R=view_dist, Az_deg=Az, figs=True)

# plot elevation profiles for few directions

Az_ix = np.arange(0, len(Az), 5)

plt.figure()

for k in range(len(Az_ix)):
    y = (terrain_profile['elev'][Az_ix[k]])
    x = np.arange(0, len(y)) * dx
    s = '%.1f deg' % terrain_profile['Az'][Az_ix[k]]
    plt.plot(x, y, '-', label=s)

plt.legend()
plt.title('terrain profiles')
plt.ylabel('elev [m]')
plt.xlabel('distance from target [m]')