#!/usr/bin/env python

"""
extract patchdock results
find new ligs closest to a reference
"""

import os
import sys
import argparse
sys.path.insert(0, '/home/kmckiern/scripts/py_general/')
from toolz import xtract

parser = argparse.ArgumentParser(description='get proper docking ligs')
parser.add_argument('--zf', type=str, help='zip file')
parser.add_argument('--ref', type=str, help='reference ligand')
args = parser.parse_args()

fn = args.zf
xtract(fn)
