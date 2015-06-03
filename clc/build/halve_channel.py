#!/usr/bin/env python

import argparse
import numpy as np
import mdtraj

parser = argparse.ArgumentParser(description='use chimera to combine protein and ligand pdbs')
parser.add_argument('--pdb', type=str, help='system pdb')
args = parser.parse_args()

def main():
    system = mdtraj.load(args.pdb)
    pi = system.topology.select('protein')
    p = system.atom_slice(pi)

    na = p.n_atoms
    first = np.arange(na/2)
    second = np.arange(na/2, na)

    p.atom_slice(first).save_pdb('p0_' + pro)
    p.atom_slice(second).save_pdb('p1_' + pro)

if __name__ == '__main__':
    main()
