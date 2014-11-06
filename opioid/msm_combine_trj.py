#!/usr/bin/env python

import mdtraj
from mdtraj import Trajectory as t
import os

# combine trajectories
trjs = [f for f in os.listdir('.') if 'trj' in f]
ts = {}
tn = []
for i in trjs:
    num = i.split('.')[0].split('j')[-1]
    tn.append(num)
    ts[i] = t.load('./trj%i.h5' % int(num))
tn.sort()
z = ts['trj0.h5']
for i in tn:
    z = z.join(ts['trj%s.h5' % i])

# trim data to have a frame every 1 ns
frames = []
for i in range(len(z)):
    if i % 10 == 0:
        frames.append(i)
zp = z.slice(frames)

# save combined data
zp.save_pdb('ns.pdb')
