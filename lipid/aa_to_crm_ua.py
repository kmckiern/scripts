#!/usr/bin/env python
"""
this script converts an all-atom lipid monomer to a charmm-style united atom representation.
"""
import argparse
import os
from forcebalance.molecule import Molecule
from parse_tools import *

parser = argparse.ArgumentParser(description='parse ff params from atom list file to out file')
parser.add_argument('--xyz_file', type=str, help='path to coordinate file')
parser.add_argument('--atm_list', type=str, help='path to list of atoms to keep')
opts = parser.parse_args()

def get_ls(list_file):
    a_map = {}
    # get list of atoms to remove.
    for line in open(list_file, 'r').readlines():
        l = line.split()
        a_map[l[1].strip()] = l[1].strip()
    return a_map

def main():
    # create molecule
    m = Molecule(opts.xyz_file)
    # get elements
    ele = np.array(m.elem)
    # get atom map
    am = get_ls(opts.atm_list)
    # get current atoms
    a_0 = am.keys()

    # get indices of atoms to keep
    keep = []
    for itr, e in enumerate(ele):
        i = itr + 1
        if e + str(i) in a_0:
            keep.append(i)

    # throw out some hydrogens
    xx = m.atom_select(np.where(keep))
    x = Molecule(opts.xyz_file)
    x.xyzs = xx.xyzs

    x.write('idk.gro')

if __name__ == '__main__':
    main()
