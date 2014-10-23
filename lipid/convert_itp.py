#!/usr/bin/env python

"""
g53a6b := g
c36ua := c

write conversion scripts
have:
- g dppc itp file (50 atms)
- c dpuc itp file (72 atms)
- c dpuc pdb

1. add hydrogens and connectivity information from c dpuc to an
2. convert charmm gui dpuc pdb to have g atomnames (matching the itp written above).
    aa: /Users/kerimckiernan/Downloads/step5_assembly.pdb
    parsed:
3. maybe think about how to add parameters for hydrogens. start with some initial guess and scan?

atomname map file: /Users/kerimckiernan/Downloads/aa_gb/at.dat
"""

import argparse
from parse_tools import get_ls

parser = argparse.ArgumentParser(description='parse ff params from atom list file to out file')
parser.add_argument('--map_file', type=str, help='name of file with atom map specifications')
parser.add_argument('--model_i', type=str, help='initial model')
parser.add_argument('--model_f', type=str, help='final model')
opts = parser.parse_args()

def main():
    # get atom type and name map
    atom_map = get_ls(opts.map_file)
    # build description dictionary for each model
    mi = model_dict(atom_map, opts.model_i)
    mf = model_dict(atom_map, opts.model_f)

    # build map between models
    

if __name__ == '__main__':
    main()
