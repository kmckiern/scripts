#!/usr/bin/env python

"""
Use chimera to convert a pdb to a mol2
"""

import os
from chimera import runCommand as rc
import argparse

parser.add_argument('--pdb', type=float, help='pdb to convert')
parser.add_argument('--rmh', action='store_true', help='remove hydrogens')
args = parser.parse_args()

pdb = args.pdb

for fn in pdb:
    pref = fn.split('.')[0]
    rc('open ' + fn)
    if args.rmh:
        rc('delete H')
    rc('write format mol2 0 ' + pref + '.mol2')
    rc('close all')
rc('stop now')
