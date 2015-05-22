#!/usr/bin/env python

import argparse
import os
import sys

parser = argparse.ArgumentParser(description='use chimera to combine protein and ligand pdbs')
parser.add_argument('--pro', type=str, help='protein pdb')
parser.add_argument('--lig', type=str, help='directory for output pdbs (have trailing /)')
parser.add_argument('--sr', type=str, help='script root', default='/home/kmckiern/scripts/')
args = parser.parse_args()

sr = args.sr
sys.path.insert(0, sr + 'py_general/')
from toolz import call_cl

def call_chimera(args):
    command = 'chimera --nogui --script \'' + args + '\''
    call_cl(command)

def main():
    p = args.pro
    l = args.lig

    comb = st + 'clc/dock/analysis/comb_pl.py ' + p + ' ' + l
    call_cl(comb)

if __name__ == '__main__':
    main()
