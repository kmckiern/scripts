#!/usr/bin/env python

"""
use chimera to generate a receptor surface file
example usage:
    >> chimera --nogui --script "/home/kmckiern/scripts/clc/dock/gen_surf.py \
           receptor.pdb"
"""

import os
from chimera import runCommand as rc
from chimera import openModels, MSMSModel
from WriteDMS import writeDMS
import sys

fn = sys.argv[1]
pref = fn.split('.')[0]

rc('open ' + fn)
rc("surf")
surf = openModels.list(modelTypes=[MSMSModel])[0]
writeDMS(surf, pref + '.dms')
rc('close all')
rc('stop now')
