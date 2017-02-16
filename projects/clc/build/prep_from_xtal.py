from chimera import runCommand as rc
import sys

_, pdb_id, of = sys.argv

rc('open pdb:' + pdb_id)
rc('split #0')
rc('select #0.1-2')
rc('write selected #0 ' + of)
rc('close all')
rc('stop now')
