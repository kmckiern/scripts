#!/usr/bin/env python

"""
for superposing two atoms over a subset of atom indices
"""

import argparse
import mdtraj
import numpy as np

parser = argparse.ArgumentParser(description='superpose molecules over subset of atoms')
parser.add_argument('--ref', type=str, help='pdb to which we are aligning')
parser.add_argument('--translate', type=str, help='pdb we are moving during alignment')
parser.add_argument('--rselect', type=str, help='atom selection string for ref')
parser.add_argument('--tselect', type=str, help='atom selection string for translate lig')
parser.add_argument('--pout', type=str, help='name of output pdb')
args = parser.parse_args()

def main():
    # load in reference and molecule being superposed
    ref = args.ref
    r = mdtraj.load(ref)
    translate = args.translate
    t = mdtraj.load(translate)

    # find subset of atoms to align
    sub = r.top.select(args.tselect)

    # align
    sposed = t.superpose(r, atom_indices=sub)
    sposed.save_pdb(args.pout)

if __name__ == '__main__':
    main()
