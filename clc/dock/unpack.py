#!/usr/bin/env python

"""
extract patchdock results
find new ligs closest to a reference
"""

import os
import sys
import argparse
sys.path.insert(0, '/home/kmckiern/scripts/py_general/')
from toolz import call_cl
from toolz import xtract

parser = argparse.ArgumentParser(description='get proper docking ligs')
parser.add_argument('--zf', type=str, help='zip file')
parser.add_argument('--ref', type=str, help='reference ligand')
args = parser.parse_args()

def main():
    fn = args.zf
    pth = '/'.join(fn.split('/')[:-1])

    # extract zip file
    # xtract(fn)

    og_pdbs = [f for f in os.listdir(pth) if '.pdb' in f]
    # pull clc config from lig 0
    call_cl(['grep', '-v', 'LIG', og_pdbs[0], '>>', 'recep.pdb'])
    # parse ligs out of docking pdbs
    for p, pdb in enumerate(og_pdbs):
        write = open('lig_' + str(p) + '.pdb', 'w+')
        call_cl(['grep', 'LIG', pdb], write)

if __name__ == '__main__':
    main()
