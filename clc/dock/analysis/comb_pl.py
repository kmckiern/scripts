#!/usr/bin/env python

"""
use chimera to combine lig and receptor into one pdb
"""

import os
from chimera import runCommand as rc
import sys

p = sys.argv[1]
ppref = p.split('/')[-1].split('.')[0]
l = sys.argv[2]
lpref = l.split('.')[0]

rc('open ' + p)
rc('open ' + l)
rc('combine #0,#1')
rc('write #2 ' + ppref + '_' + lpref + '.pdb')
rc('close all')
rc('stop now')
