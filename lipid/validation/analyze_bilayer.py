#!/usr/bin/env python

"""
feed in bilayer trj and generate output file with thermodynamic and structural characterization
example usage:
    >> python /path/to/script/dir/analyze_bilayer.py --trj /path/to/trj.trr --o ./meh.dat
"""

import argparse
import mdtraj

parser = argparse.ArgumentParser(description='full analysis of bilayer trajectories')
parser.add_argument('--trj', type=str, help='path to trj')
parser.add_argument('--o', type=str, help='name of output file')
args = parser.parse_args()

# area per lipid
def get_al():

# deuterium order parameter
def get_al():

# volume per lipid
def get_vl():

# isothermal compressibility modulus
def get_al():

# x-ray structure factor
def get_fq(): 

# diffusion constant
def get_dl():

def main():
    trj = args.trj
    out = args.o

if __name__ == '__main__':
    main()
