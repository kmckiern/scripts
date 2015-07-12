#!/usr/bin/env python

import argparse
import os
import sys
import mdtraj

parser = argparse.ArgumentParser(description='use chimera to combine protein and ligand pdbs')
parser.add_argument('--pro', type=str, help='protein pdb')
parser.add_argument('--lig', type=str, help='directory for output pdbs (have trailing /)')
parser.add_argument('--sr', type=str, help='script root', default='/home/kmckiern/scripts/')
parser.add_argument('--clean_lig', action='store_true', help='clean up ligand pdb', default=False)
parser.add_argument('--ref', type=str, help='cutoff reference residue id', default=None)
parser.add_argument('--cut', type=float, help='cutoff dist between pore and ligand (/nm)', default=1.5)
args = parser.parse_args()

sr = args.sr
sys.path.insert(0, sr + 'py_general/')
from toolz import call_cl

def call_chimera(args):
    command = 'chimera --nogui --script \'' + args + '\''
    call_cl(command)

def main():
    p = args.pro
    ppref = p.split('/')[-1].split('.')[0]
    l = args.lig
    lpref = l.split('.')[0]

    if args.clean_lig:
        mdtraj.load(l).save_pdb('tidy_' + l)

    comb = sr + 'clc/dock/analysis/comb_pl.py ' + p + ' ' + l
    call_chimera(comb)

    # keep remarks w binding info
    fn = ppref + '_' + lpref + '.pdb'
    get_remarks = 'grep REMARK ' + l + '>> ' + fn
    call_cl(get_remarks)

    ref = args.ref
    if ref != None:
        import numpy as np
        # read in combined pdb
        lp = mdtraj.load(fn)
        # find distance between reference residue and ligand
        pi = lp.topology.select('resid ' + ref + ' and name C')
        li = lp.topology.select('not protein and name S')
        lq = lp.atom_slice(li).xyz[0][0]
        pq = lp.atom_slice(pi).xyz[0][0]
        dist = np.linalg.norm(pq-lq)
        # log results
        with open('dists.dat', 'a') as f:
            f.write(lpref + ': ' + str(dist) + '\n')
        f.close()
        # remove comb pdb if beyond cutoff distance
        if dist > args.cut:
            os.remove(fn)

if __name__ == '__main__':
    main()
