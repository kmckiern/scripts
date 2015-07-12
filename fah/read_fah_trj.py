#!/bin/env python

import argparse
import mdtraj
import numpy as np

parser = argparse.ArgumentParser(description='get rmsd if fah WU to reference')
parser.add_argument('--top', type=str, help='topology file')
parser.add_argument('--ref', type=str, help='reference file')
parser.add_argument('--out', type=str, help='name of output file', default='rmsd.dat')
args = parser.parse_args()

def main():
    # read in fah trj
    t = args.top
    trj = mdtraj.load('positions.xtc', top=t)

    # strip water and ions from simulation
    pi = trj.top.select('protein')
    trj = trj.atom_slice(pi)

    # calculate rmsd of trj to xtal structure
    ref = mdtraj.load(args.ref)
    trj.superpose(ref)
    dists = mdtraj.rmsd(trj, ref)
    
    # print final rmsd value
    print dists[-1]
    # np.save(args.out, dists)

if __name__ == '__main__':
    main()
