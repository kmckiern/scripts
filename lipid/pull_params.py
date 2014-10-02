#!/usr/bin/env python
"""
this script pulls all relevant parameter specifications from a forcefield directory in order to generate a ForceBalance .itp file

parser makes lots of assumptions ... will make more general with time
"""
import argparse
import os
import subprocess as subp

parser = argparse.ArgumentParser(description='parse ff params from atom list file to out file')
parser.add_argument('--af', type=str, help='file of atoms in system')
parser.add_argument('--ff_dir', type=str, help='path to ff parameter directory')
parser.add_argument('--out', type=str, help='name of output parameter file')
opts = parser.parse_args()

def main():
    # get list of atoms for which we would like to find parameters
    atms = []
    for line in open(opts.af, 'r').readlines():
        atms.append(line.strip())
    # remove duplicates, and determine now many unique elements
    uniq_atms = set(atms)
    num_atms = len(atms)
    
    # buffer by 3 for parameter type specification and comments
    parse_length = num_atms + 3

    ff_categories = {'defaults': 3, 'atomtypes': parse_length, 'bondtypes': parse_length, 'angletypes': parse_length, 'dihedraltypes': parse_length, 'pairtypes': parse_length}

    path = opts.ff_dir

    # iterate over parameter categories
    p_cats = ff_categories.keys()
    for p_type in p_cats:
        data = subp.Popen(['grep', '-A', str(gg_categories[p_type]), ptype, path + '*'], stdout=subp.PIPE).communicate()[0])
        print 'mer', data

if __name__ == '__main__':
    main()
