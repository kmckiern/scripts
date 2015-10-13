#!/bin/env python

"""
example usage:
        python pore_occupancy.py --trj /home/harrigan/data/trek/processed/p9761/17/14/cnv.xtc --top 4xdk.pdb --out sf_h2o_p9761_17_14.dat
"""

import mdtraj
import numpy as np
import argparse
import IPython

parser = argparse.ArgumentParser(description='get water and ion occupancy of trek pore')
parser.add_argument('--trj', type=str, help='trajectory')
parser.add_argument('--top', type=str, help='reference pdb topology')
parser.add_argument('--out', type=str, help='out file', default='sf_water.dat')
args = parser.parse_args()

# parameters obtained from visual inspection in vmd ... could be better
# help to define the boundary of the selectivity filter
z_offset = 0.0
cyl_length = .5
cyl_radius = .58

# get reference indices
x = mdtraj.load(args.top)
# CA of each GLY located centerally to each selectivity filter strand
res_indxs = [102, 211, 363, 472]
atom_indxs = []
for ref in res_indxs:
    atom_indxs += [atom.index for atom in x.top.atoms if ((atom.residue.resSeq == ref) and (atom.name == 'CA'))]
ri = np.array(atom_indxs)

# reference z coordinate
ref_z = np.average(x.xyz[0][ri][:,-1])
ref_xy = [np.average(x.xyz[0][ri][:,0]), np.average(x.xyz[0][ri][:,1])]

traj = mdtraj.load(args.trj, top=args.top)
wi = [atom.index for atom in traj.topology.atoms if (atom.residue.is_water and (atom.element.symbol == 'O'))]
ki = np.array([i.index for i in traj.top.atoms_by_name('K+')])
nf = len(traj)
# time series of frame, n_h2o, n_k
time_series = np.zeros([nf, 3])
for frame in range(nf):
    time_series[frame][0] = frame
    for ndx, atom_type in enumerate([wi, ki]):
        atz = traj.xyz[frame][atom_type][:,-1]
        atz -= ref_z
        z_filter = [i for i, val in enumerate(atz) if ((val > -cyl_length) & (val < cyl_length))]
        # shift center to reference and calc norm
        atzy = traj.xyz[frame][atom_type][z_filter,0:2]
        atzy -= ref_xy
        norm_xy = np.sum(np.square(atzy), axis=1)
        xy_filter = [i for i, val in enumerate(norm_xy) if (val < cyl_radius)]
        time_series[frame][ndx+1] = len(xy_filter)
# save
np.savetxt('test_' + args.out, time_series)
