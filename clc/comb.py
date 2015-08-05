#!/usr/bin/env python

import mdtraj
from mdtraj import Trajectory as t
import os, sys
from toolz import *
import argparse

parser = argparse.ArgumentParser(description='lig-template pdb generation')
parser.add_argument('--trj_dir', type=str, help='directory of simulation trajectories')
parser.add_argument('--top', type=int, help='reference pdb topology')
parser.add_argument('--stride', type=int, help='trj subsample rate')
parser.add_argument('--cut', type=int, help='frames before this index will be cut')
args = parser.parse_args()

# combine trajectories
trjs = [f for f in os.listdir(args.trj_dir) if 'equil' in f and 'netcdf' in f] 
trjs = natural_sort(trjs, top=args.top, stride=args.stride)
ts = t.load(trjs[0])
for i in trjs[1:]:
    ts += t.load(i, top=args.top, stride=args.stride)

# save combined data
ts[cut:].save_pdb(args.trj_dir + '/out.pdb')
