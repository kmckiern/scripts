#!/usr/bin/env python

"""
extract patchdock results
find new lig closest to reference
"""

import os
import sys
import mdtraj
import argparse
sys.path.insert(0, '/home/kmckiern/scripts/py_general/')
from toolz import *

parser = argparse.ArgumentParser(description='get proper docking lig')
parser.add_argument('--zf', type=str, help='zip file')
parser.add_argument('--ref', type=str, help='reference ligand')
args = parser.parse_args()

def main():
    rstr = args.ref
    ref = mdtraj.load(rstr)
    fn = args.zf
    pth = '/'.join(fn.split('/')[:-1])

    # extract zip file
    xtract(fn)

    # easier to move into working dir
    os.chdir(pth)

    # set some final params
    od = './relevant/'
    try:
        os.mkdir(od)
    except OSError:
        pass
    rf = 'receptor.pdb'
    ligf = 'lig.pdb'

    og_pdbs = [f for f in os.listdir('.') if '.pdb' in f]
    # pull clc config from lig 0
    call_write('grep -v LIG ' + og_pdbs[0], rf)
    # rewrite to dock dir in consistent format
    recptr = mdtraj.load(rf)
    recptr.save_pdb(od + rf)

    # parse ligs out of dock pdbs
    for p, pdb in enumerate(og_pdbs):
        call_write('grep LIG ' + pdb, 'lig_' + str(p) + '.pdb')
    # get updated list of lig pdbs
    lig_pdbs = [f for f in os.listdir('.') if 'lig_' in f]
    lig_pdbs.sort()

    # join lig pdbs into one trj
    lig_trj = mdtraj.load(lig_pdbs[0])
    for lig in lig_pdbs[1:]:
        l = mdtraj.load(lig)
        lig_trj += l

    # compare distances of ligand pdbs to reference
    dists = mdtraj.lprmsd(lig_trj, ref)
    winnar = dists.argmin()

    # move receptor and closest lig to dock dir
    lig_trj[winnar].save_pdb(od + ligf)

    # write some reference info
    results = open('results.dat', 'w+')
    results.write('path: ' + pth + '\n')
    results.write('ref: ' + rstr + '\n')
    results.write('w1nn4r: ' + lig_pdbs[winnar] + '\n')
    
if __name__ == '__main__':
    main()
