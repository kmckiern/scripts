#!/usr/bin/env python
"""
use chimera for homology modeling!
example usage:
    >> chimera --script "~/scripts/clc/hmodel/hmodel.py \
           target_fasta_seq template.pdb path/out.pdb"
"""
import sys
from chimera import runCommand as rc

# read in args
_, target, template, od = sys.argv

# open target sequence
rc('open ' + target)
# open template
rc('open ' + templ8)
rc('sequence #0')

# run alignment and modelling script
rc('open /Users/kerimckiernan/Dropbox/scripts/clc/hmodel/align.py')

# save models
rc('sel 1-5')
rc('write selected ' + od + '$number_' + template)
rc('close all')
rc('stop now')
