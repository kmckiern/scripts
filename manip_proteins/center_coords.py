#!/bin/env python

import argparse

parser = argparse.ArgumentParser(description='center the coordinates of a \
     small molecule')
parser.add_argument('--pdb', type=str, help='pdb file')
parser.add_argument('--of', type=str, help='pdb out file', default='centered.pdb')
args = parser.parse_args()

def main():
    p = args.pdb
    pdb = mdtraj.load(p)
    pdb.center_coordinates()
    pdb.save_pdb(args.of)

if __name__ == '__main__':
    main()
