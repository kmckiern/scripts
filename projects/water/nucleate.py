#!/usr/bin/env python

"""
replicate a crystalline pdb.  currently only works for rectangular boxes
example usage:
    >> python /path/to/script/dir/nucleate.py --i input.pdb --o output.pdb --dims x,y,z
"""

import argparse
import mdtraj

parser = argparse.ArgumentParser(description='replicate crystal pdb to arbitrary size')
parser.add_argument('--i', type=str, help='input pdb file')
parser.add_argument('--o', type=str, help='name of output file')
parser.add_argument('--dims', type=str, help='desired pdb dimensions (comma separated)')
args = parser.parse_args()

def unit_slice(unit_cell, scalings):
    return [sx, sy, sz]

def trans(unit_cell, slices):
    return None

def main():
    # read in input
    p = mdtraj.load(args.i)

    # get unit cell vectors
    d0 = p.unitcell_vectors
    # determine needed scaling
    scale = np.array([int(i) for i in args.dims.split(',')])
    rep = np.diagonal(d0[0]/sc)

    # get atom slice according to specified scaling
    slices = unit_slice(p, rep)
    # shift slice along lattice vector
    trans(p, slices)

    # write output
    p.save_pdb(args.o)

if __name__ == '__main__':
    main()
