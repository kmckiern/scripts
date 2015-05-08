#!/usr/bin/env python

"""
steps:
0. throw out pdb atoms w d_rl > 30 A (prevent seg faulting), save mol2
1. remove hydrogens from trimmed pdb, save mol2
2. use chimera to generate receptor surface file
3. generate docking spheres
4. edit number of spheres (cut those > 10 A from ligands)
5. generate box/grid
6. run dock program
"""

import argparse
import os
import mdtraj
import sys

parser = argparse.ArgumentParser(description='automate clc docking')
parser.add_argument('--rec', type=str, help='receptor file')
parser.add_argument('--lig', type=str, help='ligand pdb')
parser.add_argument('--dist', type=float, help='dist from lig for spheres')
parser.add_argument('--sr', type=str, help='script root', default='/home/kmckiern/scripts/')
parser.add_argument('--trim_rec', action='store_true', help='trim pdb?')
args = parser.parse_args()

sr = args.sr
sys.path.insert(0, sr + 'py_general/')
from toolz import call_cl

d6bin = '/home/kmckiern/src/dock6/bin/'

def call_chimera(args):
    command = ['chimera', '--nogui', '--script', args]
    call_cl(command)

def main():
    rec = args.rec
    lig = args.lig
    d = args.dist
    pref = rec.split('.')[0]

    # if pdb is gt 15000 atoms, prob need to trim.
    if args.trim_rec:
        rt = ['python', sr + 'clc/dock/trim_pdb.py', '--rec', rec, '--lig', lig, '--d', str(d*3.0)]
        call_cl(rt)
        rec = pref + '_trim.pdb'

    # generate ligand and receptor (h and no h) mol2 files
    gm2 = sr + 'clc/dock/gen_mol2.py '
    call_chimera(gm2 + lig)
    call_chimera(gm2 + rec)
    call_chimera(gm2 + rec + ' --noH')

    # write receptor surface file
    rnoH = rec.split('.')[0] + '_noH.mol2'
    gsurf = sr + 'clc/dock/gen_surf.py '
    call_chimera(gsurf + rnoH)
    
    """# run showsphere to get sphere pdb
    # assumes i'm on not0rious (only place i have dock installed)
    show = [d6bin + 'showsphere']
    show_in = [sph, '-1', 'N', pref + '_show', 'N']
    call_cl(show, show_in)

    # read in sphere pdb
    sphurz = mdtraj.load(pref + '_show.pdb')
    # read in lig pdb
    l = mdtraj.load(lig)
    # get list of spheres within distance d of lig   
    close = get_close(sphurz, l, d)"""

if __name__ == '__main__':
    main()
