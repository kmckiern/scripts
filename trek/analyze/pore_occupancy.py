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
cyl_radius = .6
z_slice = 0.1
cyl_length = 1.0
num_bins = int(cyl_length/z_slice)
nb = 2 * num_bins

# get reference indices
x = mdtraj.load(args.top)
# CA of each GLY located centerally to each selectivity filter strand
# +2 is top, -2 is bottom of strand
res_indxs = [102, 211, 363, 472]
for i in [-2, 2]:
    for j in range(4):
        res_indxs.append(res_indxs[j]+i)
atom_indxs = []
for ref in res_indxs:
    atom_indxs += [atom.index for atom in x.top.atoms if ((atom.residue.resSeq == ref) and (atom.name == 'CA'))]
ri = np.array(atom_indxs)

traj = mdtraj.load(args.trj, top=args.top)
wi = [atom.index for atom in traj.topology.atoms if (atom.residue.is_water and (atom.element.symbol == 'O'))]
ki = np.array([i.index for i in traj.top.atoms_by_name('K+')])
nf = len(traj)
# time series of bin by frame, n_h2o, n_k
time_series = np.zeros([nb, nf, 3])
for frame in range(nf):
    # reference coordinates
    ref_z = np.average(traj.xyz[frame][ri][:,-1])
    ref_xy = [np.average(traj.xyz[frame][ri][:,0]), np.average(traj.xyz[frame][ri][:,1])]
    for ndx, atom_var in enumerate([wi,ki]):
        at_all = traj.xyz[frame][atom_var]
        atz = at_all[:,-1]
        z_min = ref_z - cyl_length
        z_max = ref_z + cyl_length
        z_filter = [i for i, val in enumerate(atz) if ((val > z_min) & (val < z_max))]
        # shift center to reference and calc norm
        atxy = at_all[z_filter,0:2]
        atxy -= ref_xy
        norm_xy = np.sum(np.square(atxy), axis=1)
        xy_filter = [i for i, val in enumerate(norm_xy) if (val < cyl_radius)]
        bins = [(num_bins + int(i)) for i in ((at_all[z_filter][xy_filter][:,-1]-ref_z)/z_slice)]
        for i in range(nb):
            time_series[i][frame][ndx+1] = bins.count(i)

# save
for ndx, bin_vals in enumerate(time_series):
    np.savetxt('out/b' + str(ndx) + '/' + args.out, bin_vals)
