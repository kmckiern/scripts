#!/usr/bin/env python

import numpy as np
import mdtraj
from mdtraj import Trajectory as t
import os, sys
import argparse

parser = argparse.ArgumentParser(description='lig-template pdb generation')
parser.add_argument('--trj_dir', type=str, help='directory of simulation trajectories')
parser.add_argument('--trj_ext', type=str, help='input trajectory extension', default='netcdf')
parser.add_argument('--out', type=str, help='output trajectory', default='out.dcd')
parser.add_argument('--top', type=str, help='reference pdb topology')
parser.add_argument('--stride1', type=int, help='individual trj subsample rate', default=1)
parser.add_argument('--stride2', type=int, help='combined trj subsample rate', default=1)
parser.add_argument('--cut', type=int, help='frames before this index will be cut', default=0)
parser.add_argument('--vs', action='store_true', help='write trj of voltage sensitive residues', default=False)
parser.add_argument('--pro', action='store_true', help='write trj of protein only', default=False)
parser.add_argument('--sr', type=str, help='script root', default='/home/kmckiern/scripts/')
args = parser.parse_args()

sr = args.sr
sys.path.insert(0, sr + 'py_general/')
from toolz import natural_sort

ext_i = '.' + args.trj_ext
td = args.trj_dir

# combine trajectories
trjs = [f for f in os.listdir(td) if ext_i in f] 
trjs = natural_sort(trjs)

ts = t.load(trjs[0], top=args.top, stride=args.stride1)
if args.vs:
    # i'm going to pad these residues by 8
    arg = ts.top.select('resid 23 to 39')
    lysglu = ts.top.select('resid 118 to 135')
    lys = ts.top.select('resid 308 to 324')
    vs = np.concatenate([arg, lysglu, lys])
    ts = ts.atom_slice(vs)
    try:
        ts[0].save_pdb('/home/kmckiern/clc/analysis/vs_dihed/pro/vs_ref.pdb')
    except:
        print 'usual protonation error, probably'
    nt = len(trjs)
    for ndx, i in enumerate(trjs[1:]):
        # for newest trj, cut end just in case write is incomplete
        if ndx + 1 == nt:
            i = i[:-2]
        new = t.load(i, top=args.top, atom_indices=vs, stride=args.stride1)
        ts += new
elif args.pro:
    pro = ts.top.select('protein')
    ts = ts.atom_slice(pro)
    ts[0].save_pdb('pro_ref.pdb')
    for i in trjs[1:]:
        ts += t.load(i, top=args.top, atom_indices=pro, stride=args.stride1)
else:
    for i in trjs[1:]:
        ts += t.load(i, top=args.top, stride=args.stride1)

# save combined data
ts[args.cut::args.stride2].save(args.out)
