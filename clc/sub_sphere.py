#!/usr/bin/env python

"""
sub-select docking spheres (bc the dock default is hella slow)
example usage:
    >> python ~/script/sub_sphere.py --i input.pdb --o output.pdb --dims x,y,z
"""

import argparse
import os
import mdtraj
import sys
sys.path.insert(0, '/'.join(os.getcwd().split('/')[:-1]) + '/py_general/')
from toolz import call_cl

parser = argparse.ArgumentParser(description='replicate crystal pdb to arbitrary size')
parser.add_argument('--sph', type=str, help='full sphere file')
parser.add_argument('--lig', type=str, help='ligand pdb')
parser.add_argument('--dist', type=float, help='dist from lig')
args = parser.parse_args()

def get_close(sphere, ligand, dist):
    

def main():
    sph = args.sph
    lig = args.lig
    d = args.dist
    pref = sph.split('.')[0]

    # run showsphere to get sphere pdb
    # assumes i'm on not0rious (only place i have dock installed)
    show = ['~/src/dock6/bin/showsphere']
    show_in = [sph, '-1', 'N', pref + '_show', 'N']
    call_cl(show, show_in)

    # read in sphere pdb
    sphurz = mdtraj.load(pref + '_show.pdb')
    # read in lig pdb
    l = mdtraj.load(lig)
    # get list of spheres within distance d of lig   
    close = get_close(sphurz, l, d)

if __name__ == '__main__':
    main()
