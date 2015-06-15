#!/usr/bin/env python

"""
for superposing two atoms over a subset of atom indices
"""

import argparse
import mdtraj
import numpy as np
import sys
sys.path.insert(0, '/home/kmckiern/scripts/py_general/')
from toolz import call_cl
import IPython

parser = argparse.ArgumentParser(description='superpose molecules over subset of atoms')
parser.add_argument('--ref', type=str, help='pdb to which we are aligning')
parser.add_argument('--translate', type=str, help='pdb we are moving during alignment')
parser.add_argument('--select', type=str, help='atom selection string')
parser.add_argument('--pout', type=str, help='name of output pdb')
args = parser.parse_args()

def main():
    # load in reference and molecule being supoerposed
    ref = args.ref
    r = mdtraj.load(ref)
    translate = args.trans
    t = mdtraj.load(translate)
    sel = args.select

    # find subset of atoms to align
    sub = r.top.select(sel)
    check = t.top.select(sel)

    # align
    if sub == check:
        sposed = t.superpose(r, atom_indices=sub)
        sposed.save_pdb(args.pout)

if __name__ == '__main__':
    main()
