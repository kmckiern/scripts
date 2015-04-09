#!/usr/bin/env python

"""
This script removes a set of atoms specified in an external file from a pdb file.

Wrote this mostly to convert AA lipids to UA lipids.
"""

import argparse
import os
import mdtraj

parser = argparse.ArgumentParser(description='remove atoms in rf from pdb and write to pdb_o')
parser.add_argument('--rf', type=str, help='file of atomnames to remove')
parser.add_argument('--pdb', type=str, help='pdb file')
parser.add_argument('--pdb_o', type=str, help='name of output pdb file')
opts = parser.parse_args()

def main():
    # get list of atoms to remove.
    rm = []
    for line in open(opts.rf, 'r').readlines():
        rm.append(line.strip())

    # read in pdb
    f = mdtraj.load(opts.pdb) 

    # remove atoms in rm
    atoms_to_keep = [a.index for a in f.topology.atoms if a.name not in rm]
    f.restrict_atoms(atoms_to_keep)
    f.save(opts.pdb_o)

if __name__ == '__main__':
    main()
