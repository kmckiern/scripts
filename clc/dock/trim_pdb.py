#!/usr/bin/env python

"""
avoid segfaults bc clc is big af
example usage:
    >> python ~/script/sub_sphere.py --rec receptor.pdb --lig lig.pdb 
           --dist 10.0
"""

import argparse
import mdtraj
import numpy as np
import sys
sys.path.insert(0, '/home/kmckiern/scripts/py_general/')
from toolz import call_cl

parser = argparse.ArgumentParser(description='trim a pdb file based on cutoff \
    from some ligand')
parser.add_argument('--rec', type=str, help='pdb to trim, usually a receptor')
parser.add_argument('--lig', type=str, help='ligand pdb')
parser.add_argument('--temp', type=str, help='temp combined pdb', default='temp.pdb')
parser.add_argument('--dist', type=float, help='dist from lig (/nm)')
args = parser.parse_args()

def main():
    rec = args.rec
    r = mdtraj.load(rec)
    lig = args.lig
    l = mdtraj.load(lig)
    d = args.dist
    temp = args.temp
    pref = rec.split('.')[0]

    # can't figure out a non-hacky way to combine pdbs.
    cp_receptor = ['head', '-n', '-1', rec]
    receptor, r_err = call_cl(cp_receptor)
    cp_ligand = ['tail', '-n', '+2', lig]
    ligand, l_err = call_cl(cp_ligand)
    tf = open(temp, 'w')
    tf.write(receptor)
    tf.write(ligand)
    tf.close()

    # get indices of receptor within distance d of ligand
    comb = mdtraj.load(temp)
    # get ligand indices
    ca = comb.n_atoms
    la = l.n_atoms
    li = np.arange(ca-la, ca)
    # find neighbors
    neighbors = mdtraj.compute_neighbors(comb, d, li)[0]
    # remove ligand from neighbor list (symmetric diff)
    n = np.setxor1d(li, neighbors)
    # easier to reset sometimes
    comb = mdtraj.load(temp)
    comb.restrict_atoms(n)
    comb.save_pdb(pref + '_trim.pdb')

if __name__ == '__main__':
    main()
