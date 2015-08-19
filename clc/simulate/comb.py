#!/usr/bin/env python

import numpy as np
import mdtraj
from mdtraj import Trajectory as t
import os, sys
import argparse
import IPython

parser = argparse.ArgumentParser(description='lig-template pdb generation')
parser.add_argument('--trj_dir', type=str, help='directory of simulation trajectories')
parser.add_argument('--trj_ext', type=str, help='simulation trajectory extension', default='netcdf')
parser.add_argument('--w_ext', type=str, help='written sim trajectory extension', default='dcd')
parser.add_argument('--top', type=str, help='reference pdb topology')
parser.add_argument('--stride1', type=int, help='individual trj subsample rate', default=1)
parser.add_argument('--stride2', type=int, help='combined trj subsample rate', default=1)
parser.add_argument('--cut', type=int, help='frames before this index will be cut', default=0)
parser.add_argument('--vs', action='store_true', help='write trj of voltage sensitive residues', default=False)
parser.add_argument('--sr', type=str, help='script root', default='/home/kmckiern/scripts/')
args = parser.parse_args()

sr = args.sr
sys.path.insert(0, sr + 'py_general/')
from toolz import natural_sort

ext_i = args.trj_ext
ext_o = args.w_ext
td = args.trj_dir

# combine trajectories
trjs = [f for f in os.listdir(td) if ext_i in f] 
trjs = natural_sort(trjs)

if args.vs:
    ts = t.load(trjs[0], top=args.top, stride=args.stride1)
    arg = ts.top.select('resid 31')
    lysglu = ts.top.select('resid 126 to 127')
    vs = np.concatenate([arg, lysglu])
    ts = ts.atom_slice(vs)
    # don't read last trj in case simulations are currently running
    for i in trjs[1:-1]:
        ts += t.load(i, top=args.top, atom_indices=vs, stride=args.stride1)
else:
    ts = t.load(trjs[0], top=args.top, stride=args.stride1)
    # don't read last trj in case simulations are currently running
    for i in trjs[1:-1]:
        ts += t.load(i, top=args.top, stride=args.stride1)


# save combined data
ts[args.cut::args.stride2].save(td + '/out.' + ext_o)
