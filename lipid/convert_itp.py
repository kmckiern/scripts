#!/usr/bin/env python

"""
script for converting files between lipid models

example usage:
    >> python /Users/kerimckiernan/Dropbox/scripts/lipid/convert_itp.py --map_file at.dat --model_i DPUC --model_f DPPC --structure_file --f_i cg_cent.gro --write_out test.gro
"""

import argparse
from parse_tools import get_ls
import sys
import re

parser = argparse.ArgumentParser(description='parse ff params from atom list file to out file')
parser.add_argument('--map_file', type=str, help='name of file with atom map specifications')
parser.add_argument('--model_i', type=str, help='initial model name')
parser.add_argument('--model_f', type=str, help='final model name')
parser.add_argument('--structure_file', action='store_true', help='convert pdb or gro file', default=False)
parser.add_argument('--itp', action='store_true', help='convert itp file', default=False)
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

# get nr for for each atom
def get_nr(itp_file, an):
    print 'nothing yet'

def map_file(in_file, atomname_map, out_file, sf=False, itp=False):
    old_keys = atomname_map.keys()
    an_field = 4
    with open(in_file, "r") as template:
        lines = template.readlines()
    with open(out_file, "w") as of:
        for line in lines:
            l = re.split(r'(\s+)', line)
            if sf:
                if len(l) > an_field + 1:
                    if l[an_field] in old_keys:
                        # whitespace conservation, right align
                        conserve_length = len(l[an_field-1] + l[an_field])
                        sub = atomname_map[l[an_field]]
                        l_diff = len(sub) - len(l[an_field])
                        # replace atom name
                        l[an_field] = sub
                        # subtract whitespace diff from prev field
                        if l_diff > 0:
                            l[an_field-1] = l[an_field-1][:-(l_diff)]
                        else:
                            for space in range(abs(l_diff)):
                                l[an_field-1] += ' '
                        of.write("".join(l))
                    else:
                        of.write(line)
                else:
                    of.write(line)
            if itp:
                print 'making sure it works so far'

def main():
    # get atom type and name map
    atom_map = get_ls(opts.map_file)
    # build description dictionary for each model, as well as cross model reference
    mi, mf, cross_an, cross_at = model_dict(atom_map, opts.model_i, opts.model_f)

    # pdb or gro conversion just maps atom names (not types)
    map_file(opts.f_i, cross_an, opts.write_out, opts.structure_file, opts.itp)

if __name__ == '__main__':
    main()
