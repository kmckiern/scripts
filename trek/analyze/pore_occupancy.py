#!/bin/env python

"""
example usage:
        python pore_occupancy.py --trj /home/harrigan/data/trek/processed/p9761/17/14/cnv.xtc --top 4xdk.pdb --out sf_h2o_p9761_17_14.dat
"""

import mdtraj
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='get water and ion occupancy of trek pore')
parser.add_argument('--trj', type=str, help='trajectory')
parser.add_argument('--top', type=str, help='reference pdb topology')
parser.add_argument('--out', type=str, help='out file', default='sf_water.dat')
args = parser.parse_args()

# via http://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays
def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.
    """
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype
    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)
    m = int(n / arrays[0].size)
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in range(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

def filter_distances(distances, cutoff, atom_pairs):
    frame_occupancy = []
    for f_num, frame in enumerate(distances):
        h2o = []
        for ndx, pair in enumerate(frame):
            if pair <= cutoff:
                h2o.append(atom_pairs[ndx][-1])
        frame_occupancy.append((f_num, len(set(h2o))))
    return frame_occupancy

x = mdtraj.load(args.top)

# CA of each GLY located centerally to each selectivity filter strand
res_indxs = [102, 211, 363, 472]
atom_indxs = []
for ref in res_indxs:
    atom_indxs += [atom.index for atom in x.top.atoms if ((atom.residue.resSeq == ref) and (atom.name == 'CA'))]
ri = np.array(atom_indxs)

trj = mdtraj.load(args.trj, top=args.top)

wi = trj.top.select('water')
w_atom_pairs = cartesian((ri, wi))
w_distances = mdtraj.compute_distances(trj, w_atom_pairs)
wo = np.array(filter_distances(w_distances, .8000001, w_atom_pairs))
np.savetxt('out/water_' + args.out, wo)

ki = np.array([i.index for i in trj.top.atoms_by_name('K+')])
k_atom_pairs = cartesian((ri, ki))
k_distances = mdtraj.compute_distances(trj, k_atom_pairs)
ko = np.array(filter_distances(k_distances, .8000001, k_atom_pairs))
np.savetxt('out/k_' + args.out, ko)
