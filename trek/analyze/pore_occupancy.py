#!/bin/env python

import mdtraj
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='get water and ion occupancy of trek pore')
parser.add_argument('--trj', type=str, help='trajectory')
parser.add_argument('--out', type=str, help='output trajectory', default='out.dcd')
parser.add_argument('--top', type=str, help='reference pdb topology')
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
    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

# GLYs located in center of each selectivity filter strand
reference_indxs = [101, 210, 362, 471]
ri = np.array(reference_indxs)

trj = mdtraj.load(args.trj, top=args.top)
wi = trj.top.select('water')

atom_pairs = cartesian((ri, wi))

distances = mdtraj.compute_distances(trj, atom_pairs)

frame_occupancy = []
for frame in distances:
    h2o = []
    for ndx, pair in enumerate(frame):
        if pair < .6:
            h2o.append(tp[ndx][-1])
    frame_occupancy.append(len(set(h2o)))
