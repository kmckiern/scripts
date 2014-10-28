#!/usr/bin/env python

"""
post process titan opioid simulations
"""

import argparse
import os
import subprocess as subp
import numpy as np
import copy

parser = argparse.ArgumentParser(description='ptraj and subscript to postprocess data')
parser.add_argument('--map_file', type=str, help='map between raw trjs and out')
parser.add_argument('--pt_app', type=str, help='ptraj.in append file')
parser.add_argument('--sub', type=str, help='sub script')
parser.add_argument('--root', type=str, help='root dir')
args = parser.parse_args()

def get_map(list_file):
    # get list of atoms to remove.
    data = []
    for line in open(list_file, 'r').readlines():
        ls.append(line.split())
    return ls

def fill_template(template_file, map_data):
    with open(template_file, "r") as template:
        lines = template.readlines()
    out = []
    for line in lines:
        l = line.split('$$')
        # if variables are present, replace them using the replacement_map
        if len(l) > 1:
            for p, entry in enumerate(l):
                if entry in map_data.keys():
                    l[p] = map_data[entry]
        out.append(''.join(l))
    return out

def main():
    # init vars
    rd = args.root
    p_app = args.pti
    ss = args.sub

    # read in map file
    map = get_map(args.map_file)
    
    for pair in map:
        lig, trj_out = pair

        replace = {'LIG': lig, 'TRJ': trj_out, 'ROOT_DIR': rd}
        
        # gen sub script
        # replace
        subs = fill_template(ss, replace)
        print subs
        # write to correct dir

        # append ligand dependent ptraj data
        # replace
        ptraj = fill_template(p_app, replace)
        # append to ptraj.in in correct dir
        print ptraj

if __name__ == '__main__':
    main()
