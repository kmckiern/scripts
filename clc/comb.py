#!/usr/bin/env python

import mdtraj
from mdtraj import Trajectory as t
import os, sys
from toolz import *

# combine trajectories
trjs = [f for f in os.listdir('.') if 'equil' in f and 'netcdf' in f] 
trjs = natural_sort(trjs)
ts = {}
for i in trjs:
    ts[i] = t.load(i, top='clc2_eq.pdb', stride=5)
z = ts[trjs[0]][-1]
for i in trjs[1:]:
    print 'z: ', z
    print 'i: ', i
    z = z + ts[i][-1]

# save combined data
z.save_pdb('out.pdb')
