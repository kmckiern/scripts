#!/bin/env python

import argparse
import mdtraj
import numpy as np

parser = argparse.ArgumentParser(description='align pdbs')
parser.add_argument('--ref', type=str, help='reference file')
parser.add_argument('--align', type=str, help='file to align')
parser.add_argument('--out', type=str, help='out file')
args = parser.parse_args()

x = np.arange(6)
y = np.arange(9,14)
z = np.concatenate((x,y))

ref = args.ref
template = mdtraj.load(ref)

align = args.align
ta = mdtraj.load(align)

tad = ta.superpose(template, atom_indices=z)
tad.save_pdb(args.out)
