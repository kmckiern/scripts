#!/bin/env python

import sys
import mdtraj

_, p0, pf = sys.argv
mdtraj.load(p0).save_pdb(pf)
