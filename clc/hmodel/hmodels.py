#!/usr/bin/env python

"""
wrapper script for calling chimera
"""

import argparse
import os
import sys

parser = argparse.ArgumentParser(description='use chimera to combine protein and ligand pdbs')
parser.add_argument('--seq', type=str, help='target fasta sequence')
parser.add_argument('--templates', type=str, help='directory of template pdbs (have trailing /)')
parser.add_argument('--sr', type=str, help='script root', default='/home/kmckiern/scripts/')
parser.add_argument('--od', type=str, help='out dir')
args = parser.parse_args()

sr = args.sr
sys.path.insert(0, sr + 'py_general/')
from toolz import call_cl

def call_chimera(args):
    # MAV requires gui
    command = 'chimera --script \'' + args + '\''
    call_cl(command)

def main():
    seq = args.seq
    pref = seq.split('/')[-1]
    templ8s = args.templates
    out = args.od + pref 

    templates = [f for f in os.listdir(templ8s) if '.pdb' in f]

    for t in templates:
        model = sr + 'clc/hmodel/align_model.py ' + seq + ' ' + templ8s + t + ' ' + out
        call_chimera(model)

if __name__ == '__main__':
    main()
