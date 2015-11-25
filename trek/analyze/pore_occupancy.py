#!/bin/env python

"""
example usage:
        python pore_occupancy.py --trj whereever/cnv.xtc --top 4xdk.pdb --out sf_h2o_p9761_17_14.dat

note: use python version >= 3.4
"""

import argparse
import mdtraj
import numpy as np
import pandas as pd
from scipy.stats import multivariate_normal
import IPython

parser = argparse.ArgumentParser(description='get water and ion occupancy of trek sf')
parser.add_argument('--tf', type=str, help='trajectory file')
parser.add_argument('--gaussian', action='store_true', help='gaussian counting', default=False)
parser.add_argument('--histogram', action='store_true', help='histogram counting', default=False)
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

# via: 
# http://stackoverflow.com/questions/25720600/generating-3d-gaussian-distribution-in-python
def numerical_gaussian(xmin, xmax, ymin, ymax, mu, sigma):
    # generate numerical grid
    x, y = np.mgrid[xmin:xmax, ymin:ymax]
    xy = np.column_stack([x.flat, y.flat])
    # get gaussian values over grid
    covariance = np.diag(sigma**2)
    gz = multivariate_normal.pdf(xy, mean=mu, cov=covariance)
    gz = g.reshape(x.shape)
    return gz
    
def relevant_indices(trj):
    # CA of each GLY located centerally to each selectivity filter strand
    # +2 is top, -2 is bottom of strand
    for i in [-3, -2, -1, 1, 2, 3]:
        for j in range(4):
            res_indxs.append(res_indxs[j]+i)
    atom_indxs = []
    for ref in res_indxs:
        atom_indxs += [atom.index for atom in trj.top.atoms if ((atom.residue.resSeq == ref) and (atom.name == 'CA'))]
    ri = np.array(atom_indxs)
    wi = [atom.index for atom in trj.topology.atoms if (atom.residue.is_water and (atom.element.symbol == 'O'))]
    ki = np.array([i.index for i in trj.top.atoms_by_name('K+')])
    return ri, wi, ki

def sf_bc(traj, ri):
    # filter by z
    ri_qz = traj.xyz[:, ri, -1]
    z_cent = np.average(ri_qz, axis=1)
    z_min = z_cent - cyl_length
    z_max = z_cent + cyl_length
    # filter by xy
    ri_qxy = traj.xyz[:, ri, :2]
    xy_cent = np.average(ri_qxy, axis=1)
    return z_min, z_max, xy_cent

def discrete_bins(traj, wi, ki, z_min, z_max, xy_cent):
    nf = len(traj)
    # time series of bin by frame, n_h2o, n_k
    time_series = np.zeros([nb, nf, 2])

    z_lbls = []
    for ndx, atom_var in enumerate([wi, ki]):
        # label above (1), within (0), or below (-1) sf wrt z
        zs = traj.xyz[:, atom_var, -1]
        z_scale = (zs.T - z_min).T
        z_scale[(z_scale >= 0.0) & (z_scale <= full_l)] = 0
        z_scale[z_scale < 0.0] = -1
        z_scale[z_scale > full_l] = 1
        z_lbls.append(z_scale)
    
        # label if within radius of sf (0), or outside (1)
        xys = traj.xyz[:,atom_var,:2]
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
    for t_ndx, trjd in enumerate(trj_data):
        # read, load, align
        t, top = trjd.split()
        label = t.strip().split('/')[-1].split('.')[0]
        traj = mdtraj.load(t, top=top)
        pi = traj.top.select('protein')
        traj = traj.superpose(traj, atom_indices=pi)
        x = traj[0]

        ri, wi, ki = relevant_indices(x)
        # get ions and water in the selectivity filter
        if args.gaussian:
            # fit 10 gaussians to SF 
            gaussians = fit_gaussians()
            # bin gaussian counts
            bin_timeseries, zs = gaussian_bins(*gaussians)

        if args.histogram:
            # get trj info
            z_min, z_max, xy_cent = sf_bc(traj, ri)
            # histogram data
            bin_timeseries, zs = discrete_bins(traj, wi, ki, z_min, z_max, xy_cent)

        # record bin data
        for ndx, i in enumerate(bin_timeseries):
            if ndx == 0:
                write_bins(label, zs[0], 'atom_resolved_bins/water_')
            if ndx == 1:
                write_bins(label, zs[1], 'atom_resolved_bins/k_')
        write_bins(label, bin_timeseries, 'bins/b')
        # sum over bins for feature vector of SF occupancy
        first = args.start
        last = args.end
        bin_window = sum_bins(first, last, bin_timeseries)
        write_bins(label, bin_window, 'bins/summed/w_' + str(first) + '_' + str(last))

if __name__ == '__main__':
    main()
