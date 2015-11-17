#!/bin/env python

"""
example usage:
        python pore_occupancy.py --trj whereever/cnv.xtc --top 4xdk.pdb --out sf_h2o_p9761_17_14.dat
"""

import mdtraj
import numpy as np
import argparse
import IPython

parser = argparse.ArgumentParser(description='get water and ion occupancy of trek pore')
parser.add_argument('--trj', type=str, help='trajectory')
parser.add_argument('--top', type=str, help='reference pdb topology')
parser.add_argument('--out', type=str, help='out file', default='sf_water.dat')
parser.add_argument('--record', action='store_true', help='whether to save first frame', default=False)
args = parser.parse_args()

# parameters obtained from visual inspection in vmd ... could be better
# help to define the boundary of the selectivity filter
cyl_radius = .6
z_slice = 0.1
cyl_length = 1.0
num_bins = round(cyl_length/z_slice)
nb = 2 * num_bins

traj = mdtraj.load(args.trj, top=args.top)
# get reference indices
x = traj[0]
# CA of each GLY located centerally to each selectivity filter strand
# +2 is top, -2 is bottom of strand
res_indxs = [101, 210, 362, 471]
for i in [-2, -1, 1, 2]:
    for j in range(4):
        res_indxs.append(res_indxs[j]+i)
atom_indxs = []
for ref in res_indxs:
    atom_indxs += [atom.index for atom in x.top.atoms if ((atom.residue.resSeq == ref) and (atom.name == 'CA'))]
ri = np.array(atom_indxs)

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
        z_min = ref_z - cyl_length
        z_max = ref_z + cyl_length
        z_filter = [i for i in atom_var if ((traj.xyz[frame][i][-1]) > z_min) & (traj.xyz[frame][i][-1] < z_max)]
        # shift center to reference and calc norm
        filtered = traj.xyz[frame][z_filter][:,:2]
        filtered -= ref_xy
        norm_xy = np.sum(np.square(filtered), axis=1)
        final = [z_filter[ndx] for ndx, i in enumerate(norm_xy) if i < cyl_radius]
        final_z = traj.xyz[frame][final][:,-1]

        # histogram the results
        tval, tbin = np.histogram(final_z, bins=nb, range=(z_min, z_max))
        time_series[:, frame, ndx+1] = tval

        if args.record:
                # this would be useful for tracking waters and ions:
                # sf_atom_names = [x.top.atom(atom_var[i]) for i in xy_filter]
                sf_atom_indxs = final + atom_indxs
                if frame == 0:
                    sf = traj[frame].atom_slice(sf_atom_indxs)
                    sf.save_pdb('record/sf_' + args.out.split('.')[0] + '.pdb')	            	
                    x.save_pdb('record/full_' + args.out.split('.')[0] + '.pdb')	            	

# save
for ndx, bin_vals in enumerate(time_series):
    np.savetxt('out/b' + str(ndx) + '/' + args.out, bin_vals)
