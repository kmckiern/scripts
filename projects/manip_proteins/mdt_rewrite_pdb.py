#!/bin/env python

import mdtraj
import argparse

parser = argparse.ArgumentParser(description='for rewriting pdbs by mdtraj')
parser.add_argument('--p0', type=str, help='pdb in', default='protein.pdb')
parser.add_argument('--pf', type=str, help='pdb out', default='out.pdb')
args = parser.parse_args()

p0 = args.p0
pf = args.pf

mdtraj.load_pdb(p0).save_pdb(pf)
