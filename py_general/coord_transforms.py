#!/usr/bin/env python

"""
this is quick and dirty
generalize w euler matrices or quaternions later
"""

import numpy as np

# assumes coordinates are of the mdtraj.load('something.pdb').xyz variety
# and theta is in radians
def about_z(coords, theta):
    qs = coords[0]
    x = qs[:,0]
    y = qs[:,1]
    xp = x*np.cos(theta) - y*np.sin(theta)
    yp = x*np.sin(theta) + y*np.cos(theta)
    qs[:,0] = xp
    qs[:,1] = yp
    coords[0] = qs
    return coords
