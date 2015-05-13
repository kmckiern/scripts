#!/usr/bin/env python

"""
extract patchdock results
find new ligs closest to a reference
"""

import os
import sys
import argparse
sys.path.insert(0, '/home/kerimckiernan/Dropbox/scripts/py_general/')
from toolz import xtract

parser = argparse.ArgumentParser(description='get proper docking ligs')
parser.add_argument('--ref', type=str, help='reference ligand')
args = parser.parse_args()

fn = sys.argv[1]
xtract(fn)
