#!/usr/bin/env python

"""
given a swissdock cluster pdb, pull highest rank lig from each cluster and write pdb
"""

import argparse
import os
import sys
import mdtraj

parser = argparse.ArgumentParser(description='pull best ligs from sw output')
parser.add_argument('--clust', type=str, help='sw cluster pdb file')
parser.add_argument('--out', type=str, help='directory for output pdbs (have trailing /)')
parser.add_argument('--sr', type=str, help='script root', default='/home/kmckiern/scripts/')
args = parser.parse_args()

sr = args.sr
sys.path.insert(0, sr + 'py_general/')
from toolz import call_cl

def get_num(cmd):
    return int(call_cl(cmd)[0].strip())

def line_range(clust_line, len_lig, clust_position):
    b4 = clust_line - clust_position + 1
    aftr = len_lig - clust_position + clust_line
    return (b4, aftr)

def main():
    ligs = args.clust
    od = args.out
    pref = ligs.split('/')[-1]

    # get length of one pdb w comments
    ll = 'grep -n \"TER\" ' + ligs + ' | head -n 1 | awk -F \':\' \'{print $1}\''
    lenl = get_num(ll)-1

    # get number of clusters
    cn = 'grep \"Cluster:\" ' + ligs + ' | tail -n 1 | awk -F \':\' \'{print $2}\''
    nc = get_num(cn)

    # loop over clusters
    for i in range(nc):
        # remove pdb if it exists
        new_pdb = od + 'c' + str(i) + '_' + pref
        try:
            os.remove(new_pdb)
        except OSError:
            pass
        # create new pdb
        place = 'grep -n \"Cluster: ' + str(i) + '$\" ' + ligs + ' | head -n 1 | awk -F \':\' \'{print $1}\''
        clust = get_num(place)
        if i == 0:
            clust0 = clust
        b4, aftr = line_range(clust, lenl, clust0)
        pull = 'sed -n \'' + str(b4) + ',' + str(aftr) + 'p;' + str(aftr) + 'q\' ' + ligs + '>> ' + new_pdb
        call_cl(pull)
        end = 'echo \"END\" >> ' + new_pdb
        call_cl(end)

if __name__ == '__main__':
    main()
