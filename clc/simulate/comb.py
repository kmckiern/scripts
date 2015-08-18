#!/usr/bin/env python

import mdtraj
from mdtraj import Trajectory as t
import numpy as np
import os, sys
import argparse

parser = argparse.ArgumentParser(description='lig-template pdb generation')
parser.add_argument('--trj_dir', type=str, help='directory of simulation trajectories')
parser.add_argument('--trj_ext', type=str, help='simulation trajectory extension', default='netcdf')
parser.add_argument('--w_ext', type=str, help='written sim trajectory extension', default='dcd')
parser.add_argument('--top', type=str, help='reference pdb topology')
parser.add_argument('--stride1', type=int, help='individual trj subsample rate')
parser.add_argument('--stride2', type=int, help='combined trj subsample rate')
parser.add_argument('--cut', type=int, help='frames before this index will be cut')
parser.add_argument('--vs', action='store_true', help='write trj of voltage sensitive residues')
parser.add_argument('--sr', type=str, help='script root', default='/home/kmckiern/scripts/')
args = parser.parse_args()

sr = args.sr
sys.path.insert(0, sr + 'py_general/')
from toolz import natural_sort

ext_i = args.trj_ext
ext_o = args.w_ext

# combine trajectories
trjs = [f for f in os.listdir(args.trj_dir) if ext in f] 
trjs = natural_sort(trjs)
ts = t.load(trjs[0], top=args.top, stride=args.stride1)
for i in trjs[1:]:
    ts += t.load(i, top=args.top, stride=args.stride1)

if args.vs:
    arg = x.top.select('resid 31')
    lysglu = x.top.select('resid 126 to 127')
    vs = np.concatenate([arg, lysglu])
    ts = ts.atom_slice(vs)

# save combined data
ts[args.cut::args.stride2].save_pdb(args.trj_dir + '/out.pdb')
