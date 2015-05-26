#!/usr/bin/env python
"""
use chimera for homology modeling!
example usage:
    >> chimera --script "~/scripts/clc/hmodel/hmodel.py \
           target_fasta_seq template.pdb path/out.pdb"
"""
import sys
from chimera import runCommand as rc
import time

# read in args
_, target, template, od = sys.argv

# open target sequence
rc('open ' + target)
# open template structure, and generate sequence
rc('open ' + template)
rc('sequence #0')

# run alignment and modelling script
rc('open /Users/kerimckiernan/Dropbox/scripts/clc/hmodel/align_model.py; wait')

# save models
rc('preset apply publication 4')
#rc('write 1 ' + od + '$number_' + template)
#rc('close all')
#rc('stop now')
