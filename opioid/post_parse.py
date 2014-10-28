#!/usr/bin/env python

"""
post process titan opioid simulations
"""

import argparse
import os
import subprocess as subp
import numpy as np
import copy
from fmt_path import format_path

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
        data.append(line.split())
    return data

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
    p_app = args.pt_app
    ss = args.sub

    # read in map file
    map = get_map(args.map_file)
    
    for pair in map:
        lig, trj_out = pair

        replace = {'LIG': lig, 'TRJ': trj_out, 'ROOT_DIR': rd}

        lig_dir = format_path(rd, lig)
        sub_dir = format_path(rd, 'sub')
        
        # replace sub script
        subs = fill_template(ss, replace)
        # write to sub dir
        with open(sub_dir + lig + '.sh', 'w') as out_sub:
            for ln in subs:
                out_sub.write(ln)

        # replace ligand dependent ptraj data
        ptraj = fill_template(p_app, replace)
        # append

if __name__ == '__main__':
    main()
