#!/usr/bin/env python

"""
because the monomer is hard to align in the membrane given it's symmetry
"""

import argparse
import mdtraj
import numpy as np

parser = argparse.ArgumentParser(description='align 3nmo')
parser.add_argument('--pdb', type=str, help='system pdb')
parser.add_argument('--ref', type=str, help='reference halved dimer pdb')
parser.add_argument('--out', type=str, help='out file')
parser.add_argument('--trim_ref', action='store_true', help='trim hydrogens from ref pdb?')
args = parser.parse_args()

def main():
    monomer = args.pdb
    ref = args.ref
    m = mdtraj.load(monomer)
    r = mdtraj.load(ref)

    if args.trim_ref:
        r_not_H = r.topology.select('name != hydrogen')
        r = r.atom_slice(r_not_h)

    m_sup = m.superpose(r, atom_indices=np.arange(r.n_atoms))
    m_sup.save_pdb(args.out)

if __name__ == '__main__':
    main()
