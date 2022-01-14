# horizon height and spatial radiation fields

Computes horizon height for a point in DEM (Python 3.7)

First, computes elevation angle between target point P0 and all DEM points using trigonometry. Then, uses Bresenham's line algorithm to find DEM cells that are traversed through when looking into a given direction.

Returns horizon height as function of Azimuth angle (0deg = N, increases clockwise) from P0, and terrain profiles.

Development in progress.
