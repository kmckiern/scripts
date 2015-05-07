#!/usr/bin/env python

"""
avoid segfaults bc clc is big af
example usage:
    >> python ~/script/sub_sphere.py --rec receptor.pdb --lig lig.pdb 
           --dist 10.0
"""

import argparse
import mdtraj
sys.path.insert(0, '~/scripts/py_general/')
from toolz import call_cl

parser = argparse.ArgumentParser(description='trim a pdb file based on cutoff 
    from some ligand')
parser.add_argument('--rec', type=str, help='pdb to trim, usually a receptor')
parser.add_argument('--lig', type=str, help='ligand pdb')
parser.add_argument('--dist', type=float, help='dist from lig (/angstrom)')
args = parser.parse_args()

def main():
    rec = args.rec
    lig = args.lig
    d = args.dist
    pref = rec.split('.')[0]

    # can't figure out a non-hacky way to combine pdbs.
    temp = 'temp.pdb'
    cp_temp = ['cp', rec, temp]
    call_cl(cp_temp)
    rm_end = ['sed', '\'$d\'', temp]
    call_cl(rm_end)
    append_lig = ['cat', lig, '>>', temp]
    call_cl(append_lig)
    new_end = ['echo', 'END', '>>', temp]
    

if __name__ == '__main__':
    main()
