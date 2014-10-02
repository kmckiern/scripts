#!/usr/bin/env python

"""
this script pulls all relevant parameter specifications from a forcefield directory in order to generate a ForceBalance .itp file.

parser makes lots of assumptions ... will make more general with time.
"""

import argparse
import os

parser = argparse.ArgumentParser(description='remove atoms in rf from pdb and write to pdb_o')
parser.add_argument('--af', type=str, help='file of atoms in system')
parser.add_argument('--ff_dir', type=str, help='path to ff parameter directory')
parser.add_argument('--out', type=str, help='name of output parameter file')
opts = parser.parse_args()

ff_categories = {'defaults': [], 'atomtypes': [], 'bondtypes': [], 'angletypes': [], 'dihedraltypes': [], 'pairtypes': []}

def main():
    # get list of atoms for which we would like to find parameters
    atms = []
    for line in open(opts.af, 'r').readlines():
        atms.append(line.strip())
    num_atms = len(atms)

    # read in pdb
    f = mdtraj.load(opts.pdb) 

    # remove atoms in rm
    atoms_to_keep = [a.index for a in f.topology.atoms if a.name not in rm]
    f.restrict_atoms(atoms_to_keep)
    f.save(opts.pdb_o)

if __name__ == '__main__':
    main()
