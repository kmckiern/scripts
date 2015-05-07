#!/usr/bin/env python

"""
Use chimera to convert a pdb to a mol2
"""

import os
from chimera import runCommand as rc
import sys

fn = sys.argv[1]

pref = fn.split('.')[0]
rc('open ' + fn)
if len(sys.argv) > 2:
    rc('delete H')
    rc('write format mol2 0 ' + pref + '_noH.mol2')
else:
    rc('write format mol2 0 ' + pref + '.mol2')
rc('close all')
rc('stop now')
