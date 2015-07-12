import sys
from chimera import runCommand as rc

_, pdb_in, pdb_out = sys.argv

rc('open ' + pdb_in)
rc('write #0 ' + pdb_out)
rc('close all')
rc('stop now')
