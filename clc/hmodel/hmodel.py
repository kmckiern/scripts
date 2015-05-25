#!/usr/bin/env python

"""
use chimera for homology modeling!
example usage:
    >> chimera --nogui --script "~/scripts/clc/hmodel/hmodel.py \
           target_fasta_seq template.pdb path/out.pdb"
"""

import os
import sys
from chimera import runCommand as rc

target = sys.argv[0]
templ8 = sys.argv[1]
out = sys.argv[2]

# open target sequence
rc('open ' + target)
# open template
rc('open ' + templ8)


rc('close all')
rc('stop now')
