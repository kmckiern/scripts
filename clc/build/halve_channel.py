#!/usr/bin/env python

import argparse
import numpy as np
import mdtraj

parser = argparse.ArgumentParser(description='use chimera to combine protein and ligand pdbs')
parser.add_argument('--pro', type=str, help='protein pdb')
args = parser.parse_args()

def main():
    pro = args.pro
    p = mdtraj.load(pro)

    na = p.n_atoms
    first = np.arange(na/2)
    second = np.arange(na/2, na)

    p.atom_slice(first).save_pdb('p0_' + pro)
    p.atom_slice(second).save_pdb('p1_' + pro)

if __name__ == '__main__':
    main()