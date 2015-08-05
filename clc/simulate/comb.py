#!/usr/bin/env python

import mdtraj
from mdtraj import Trajectory as t
import os, sys
import argparse
import IPython

parser = argparse.ArgumentParser(description='lig-template pdb generation')
parser.add_argument('--trj_dir', type=str, help='directory of simulation trajectories')
parser.add_argument('--top', type=str, help='reference pdb topology')
parser.add_argument('--stride', type=int, help='trj subsample rate')
parser.add_argument('--cut', type=int, help='frames before this index will be cut')
parser.add_argument('--sr', type=str, help='script root', default='/home/kmckiern/scripts/')
args = parser.parse_args()

sr = args.sr
sys.path.insert(0, sr + 'py_general/')
from toolz import natural_sort

# combine trajectories
trjs = [f for f in os.listdir(args.trj_dir) if 'netcdf' in f] 
trjs = natural_sort(trjs)
ts = t.load(trjs[0], top=args.top, stride=args.stride)
for i in trjs[1:]:
    ts += t.load(i, top=args.top, stride=args.stride)

# save combined data
ts[cut:].save_pdb(args.trj_dir + '/out.pdb')
