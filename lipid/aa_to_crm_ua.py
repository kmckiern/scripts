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

def main():
    M = Molecule(opts.xyz_file)

    atm_ls = get_ls(opts.rf)

    

if __name__ == '__main__':
    main()
