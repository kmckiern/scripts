#!/bin/env python

from subprocess import Popen, PIPE, STDOUT
import json
import sys
import os

def call_cl(command_lst, std=PIPE, pipe_args=[]):
    p = Popen(command_lst, stdout=PIPE, stdin=PIPE, stderr=STDOUT, shell=True)
    out, err = p.communicate(input='\n'.join(pipe_args))
    return out, err

trj_file = open(sys.argv[1], 'r')
tfs = [line.strip() for line in trj_file.readlines()]

for t in tfs:
    trj, top = t.split()

    tsplit = trj.split('/')
    label = '_'.join(tsplit[6:-2]) + '_' + tsplit[-1].split('.')[0] + '.dat'

    cmd = 'python pore_occupancy.py --trj ' + trj + ' --top ' + top + ' --out ' + label

    # sanity
    # print ('cmd: ', cmd)
    out, err = call_cl(cmd)
