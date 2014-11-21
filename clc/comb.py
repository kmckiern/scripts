#!/usr/bin/env python

import mdtraj
from mdtraj import Trajectory as t
import os, sys
from toolz import *

# combine trajectories
trjs = [f for f in os.listdir('.') if 'equil' in f and 'netcdf' in f] 
ts = {}
tn = []
for i in trjs:
    num = i.split('.')[0].split('j')[-1]
    tn.append(num)
    ts[i] = t.load(i, top='clc2_eq.pdb', stride=5)
tn = natural_sort(tn)
z = ts[tn[0]][-1]
for i in tn[1:]:
    print 'z: ', z
    print 'i: ', i
    z = z + ts[i][-1]

# save combined data
zp.save_pdb('out.pdb')
