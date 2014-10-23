#!/usr/bin/env python

"""
script for converting files between models

mostly written for the following two:
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
parser.add_argument('--model_i', type=str, help='initial model name')
parser.add_argument('--model_f', type=str, help='final model name')
parser.add_argument('--structure_file', action='store_true', help='convert pdb or gro file')
parser.add_argument('--itp', action='store_true', help='convert itp file')
parser.add_argument('--f_i', type=str, help='initial model file')
parser.add_argument('--f_f', type=str, help='initial model file', default=None)
parser.add_argument('--write_out', type=str, help='name for converted file')
opts = parser.parse_args()

def model_dict(am, mdl_i, mdl_f):
    li, lf, ca, cat = {}, {}, {}, {}
    for i, line in enumerate(am):
        l = line.split()
        if i == 0:
            # if atom name or type dne for mapped model, n/t is kept the same (eg hg hydrogens)
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

def map_sf(in_file, atomname_map, out_file):
    old_keys = atomname_map.keys()
    with open(in_file, "r") as template:
        lines = template.readlines()
    with open(out_file, "w") as of:
        for line in lines:
            l = line.split()
            if len(l) > 1:
                if l[1] in old_keys:
                    l[1] = atomname_map[l[1]]
                    of.write("\t" + "\t".join(l) + "\n")
                else:
                    of.write(line)
            else:
                of.write(line)

def main():
    # get atom type and name map
    atom_map = get_ls(opts.map_file)
    # build description dictionary for each model, as well as cross model reference
    mi, mf, cross_an, cross_at = model_dict(atom_map, opts.model_i, opts.model_f)

    # pdb or gro conversion just maps atom names (not types)
    if opts.structure_file:
        map_sf(opts.f_i, cross_an, opts.write_out)

if __name__ == '__main__':
    main()
