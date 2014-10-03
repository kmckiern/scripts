#!/usr/bin/env python
"""
this script converts an all-atom lipid monomer to a charmm-style united atom representation.
"""
import argparse
from forcebalance.molecule import Molecule
import numpy as np
from parse_tools import *

parser = argparse.ArgumentParser(description='parse ff params from atom list file to out file')
parser.add_argument('--xyz_file', type=str, help='path to coordinate file')
parser.add_argument('--atm_list', type=str, help='path to list of atoms to keep')
opts = parser.parse_args()

def main():
    # create molecule
    m = Molecule(opts.xyz_file)
    # get elements
    ele = np.array(m.elem)
    # get atom map
    am = get_ls(opts.atm_list)
    indx = [int(i)-1 for i in am]

    print indx

    # throw out some hydrogens
    xx = m.atom_select(np.where(indx))
    x = Molecule(opts.xyz_file)
    x.xyzs = xx.xyzs

    x.write('idk.gro')

if __name__ == '__main__':
    main()
