#!/bin/env python

'''
ex: python renumber.py random.pdb test.pdb
'''

import sys
import re
import mdtraj
import IPython

f = open(sys.argv[1], 'r').readlines()
fo = open(sys.argv[2], 'w')
for line in f:
    l = re.split(r'(\s+)', line)
    if len(l) == 0:
        fo.write(line)
    elif l[0] == 'ATOM':
        current = int(l[10])
        new = current + 83
        spacep = l[9]
        spaces = l[11]
        current = str(current)
        new = str(new)
        strip = ''.join([spacep, current, spaces])
        space_diff = len(current)-len(new)
        if space_diff == 0:
            upd8 = ''.join([spacep, new, spaces])
        else:
            ds = space_diff - 1
            new_space = spacep[:ds]
            upd8 = ''.join([new_space, new, spaces])
        line = line.replace(strip, upd8)
        fo.write(line)
        continue
    else:
        fo.write(line)
fo.close()

x = mdtraj.load(sys.argv[2])
x.save_pdb(sys.argv[2])

