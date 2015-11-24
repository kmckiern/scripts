#!/bin/env python

"""
example usage:
        python pore_occupancy.py --trj whereever/cnv.xtc --top 4xdk.pdb --out sf_h2o_p9761_17_14.dat
"""

import mdtraj
import numpy as np
import argparse
import glob
from natsort import natsorted
from msmbuilder.dataset import dataset
import pandas as pd
import os.path
import IPython

parser = argparse.ArgumentParser(description='get water and ion occupancy of trek sf')
parser.add_argument('--tf', type=str, help='trajectory file')
parser.add_argument('--start', type=int, help='start bin for sum window')
parser.add_argument('--end', type=int, help='end bin for sum window')
args = parser.parse_args()

# parameters obtained from visual inspection in vmd ... could be better
# help to define the boundary of the selectivity filter
cyl_radius = .55
z_slice = 0.2
cyl_length = 1.0
num_bins = round(cyl_length/z_slice)
full_l = 2.0*cyl_length
nb = 2 * num_bins

res_indxs = [101, 210, 362, 471]

def init_trj(trajectory, topology):
    # load trjs
    traj = mdtraj.load(trajectory, top=topology)

    pi = traj.top.select('protein')
    traj = traj.superpose(traj, atom_indices=pi)
    
    # get reference indices
    x = traj[0]
    # CA of each GLY located centerally to each selectivity filter strand
    # +2 is top, -2 is bottom of strand
    for i in [-2, -1, 1, 2]:
        for j in range(4):
            res_indxs.append(res_indxs[j]+i)
    atom_indxs = []
    for ref in res_indxs:
        atom_indxs += [atom.index for atom in x.top.atoms if ((atom.residue.resSeq == ref) and (atom.name == 'CA'))]
    ri = np.array(atom_indxs)
    wi = [atom.index for atom in traj.topology.atoms if (atom.residue.is_water and (atom.element.symbol == 'O'))]
    ki = np.array([i.index for i in traj.top.atoms_by_name('K+')])
    
    # filter by z
    ri_qz = traj.xyz[:, ri, -1]
    z_cent = np.average(ri_qz, axis=1)
    z_min = z_cent - cyl_length
    z_max = z_cent + cyl_length
    # filter by xy
    ri_qxy = traj.xyz[:,ri,:2]
    xy_cent = np.average(ri_qxy, axis=1)

    return traj, wi, ki, z_min, z_max, xy_cent

def discrete_bins(traj, wi, ki, z_min, z_max, xy_cent):
    nf = len(traj)
    # time series of bin by frame, n_h2o, n_k
    time_series = np.zeros([nb, nf, 2])

    z_lbls = []
    for ndx, atom_var in enumerate([wi, ki]):
        zs = traj.xyz[:, atom_var, -1]
        # label above (1), within (0), or below (-1) sf wrt z
        z_scale = (zs.T - z_min).T
        z_scale[(z_scale >= 0.0) & (z_scale <= full_l)] = 0
        z_scale[z_scale < 0.0] = -1
        z_scale[z_scale > full_l] = 1
        z_lbls.append(z_scale)
    
        xys = traj.xyz[:,atom_var,:2]
        # label if within radius of sf (0), or outside (1)
        # broadcasting is a little weird
        sx_scale = (xys[:,:,0].T - xy_cent[:,0]).T
        sy_scale = (xys[:,:,1].T - xy_cent[:,1]).T
        moduli = sx_scale**2 + sy_scale**2
        moduli[moduli <= cyl_radius] = 0
        moduli[moduli > cyl_radius] = 1
    
        for frame in range(nf):
            z_good = np.where(z_scale[frame]==0)[0]
            xy_good = np.where(moduli[frame]==0)[0]
            in_sf = np.array(list(set(z_good).intersection(xy_good)))
    
            final_z = [traj.xyz[frame, atom_var[i], -1] for i in in_sf]
    
            # histogram the results
            tval, tbin = np.histogram(final_z, bins=nb, range=(z_min[frame], z_max[frame]))
            time_series[:, frame, ndx] = tval

    return time_series, z_lbls

# if state = end = None, full matrix will be summed
def sum_bins(start, end, timeseries):
     totaled = np.sum(timeseries[start:end], axis=0)
     return totaled

# tsave
def write_bins(of, time_series, out_pref):
    if len(time_series.shape) == 3:
        for ndx, bin_vals in enumerate(time_series):
            np.save(out_pref + str(ndx) + '/' + of, bin_vals)
    else:
        np.save(out_pref + '_' + of, time_series)
    
def main():
    trj_data = open(args.tf).readlines() 
    for t_ndx, trj in enumerate(trj_data):
        # label stuff
        label = trj.strip()

        # read in and histogram ions and water molecules
        t = 'sf_trjs/' + label + '.xtc'
        top = 'tops/p9761_79_12.pdb'

        # get trj info
        trj_geo = init_trj(t, top)
        # histogram data
        bin_timeseries, zs = discrete_bins(*trj_geo)
        # write bin data
        for ndx, i in enumerate(bin_timeseries):
            if ndx == 0:
                write_bins(label, zs[0], 'atom_resolved_bins/water_')
            if ndx == 1:
                write_bins(label, zs[1], 'atom_resolved_bins/k_')
        # write bin data
        write_bins(label, bin_timeseries, 'bins/b')
        # sum over bins for feature vector of SF occupancy
        first = args.start
        last = args.end
        bin_window = sum_bins(first, last, bin_timeseries)
        write_bins(label, bin_window, 'bins/summed/w_' + str(first) + '_' + str(last))

if __name__ == '__main__':
    main()
