#!/usr/bin/env python

"""
reimage a molecule to the center of a pbc and cut beyond a certain radius
example usage:
    >> python /path/to/script/dir/image.py --p input.pdb --t input.dcd --m mol_num --cd cutdist -o out
"""

import argparse
import mdtraj
from copy import copy

parser = argparse.ArgumentParser(description='replicate crystal pdb to arbitrary size')
parser.add_argument('--p', type=str, help='input pdb file')
parser.add_argument('--t', type=str, help='input trj')
parser.add_argument('--m', type=int, help='molecule to center')
parser.add_argument('--cd', type=str, help='distance beyond which molecules are cut')
parser.add_argument('--o', type=str, help='output pdb and trj name')
args = parser.parse_args()

#def center_mol(trj, molnum):
    
def move(trj, ndx):
    tm = copy(trj)
    tm.xyz -= t.unitcell_vectors[0][:,ndx]
    tp = copy(trj)
    tp.xyz += t.unitcell_vectors[0][:,ndx]
    return trj.stack(tm).stack(tp)

def image_trj(trj):
    # stack x for first row
    rank1 = trj.stack(move(trj, 0))
    # stack y for matrix    
    rank2 = rank1.stack(move(rank1, 1))
    # now, 3d
    rank3 = rank2.stack(move(rank2, 2))
    return rank3

def main():
    # get input
    p = mdtraj.load(args.t, top=args.p)

    # image
    p_im = image_trj(p)
    p_im[0].save_pbc('maybe.pdb')

if __name__ == '__main__':
    main()
