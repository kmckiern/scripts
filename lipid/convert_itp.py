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
3. maybe think about how to add parameters for hydrogens. start with some initial guess and scan?
"""

import argparse
from parse_tools import get_ls
import sys

parser = argparse.ArgumentParser(description='parse ff params from atom list file to out file')
parser.add_argument('--map_file', type=str, help='name of file with atom map specifications')
parser.add_argument('--model_i', type=str, help='initial model')
parser.add_argument('--model_f', type=str, help='final model')
opts = parser.parse_args()

def model_dict(am, mdl_i, mdl_f):
    li, lf, ca, cat = {}, {}, {}, {}
    for i, line in enumerate(am):
        l = line.split()
        if i == 0:
            if l[0] == mdl_i and l[3] == mdl_f:
                i_ndx = (1,2)
                f_ndx = (-2,-1)
            elif l[3] == mdl_i and l[0] == mdl_f:
                f_ndx = (1,2)
                i_ndx = (-2,-1)
            else:
                print 'fix your map file'
                sys.exit()
        # kind of gross but w/e
        li[l[i_ndx[0]]] = l[i_ndx[1]]
        lf[l[f_ndx[0]]] = l[f_ndx[1]]
        ca[l[i_ndx[0]]] = l[f_ndx[0]]
        cat[l[i_ndx[1]]] = l[f_ndx[1]]
    return li, lf, ca, cat

def main():
    # get atom type and name map
    atom_map = get_ls(opts.map_file)
    # build description dictionary for each model
    mi, mf, cross_an, cross_at = model_dict(atom_map, opts.model_i, opts.model_f)

    # build map between models

if __name__ == '__main__':
    main()
