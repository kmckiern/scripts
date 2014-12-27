#!/usr/bin/env python

# parse ncbi alignment file for modeller

import subprocess as subp
import os
import argparse
import numpy as np
from fmt_path import format_path

parser = argparse.ArgumentParser(description='parse ncbi ms alignment file for modeller')
parser.add_argument('--a', type=str, help='ncbi alignment file')
parser.add_argument('--s', type=str, help='sequence file')
args = parser.parse_args()

# this prob isn't pythonic but w/e.
def structures(af):
    entry_data = subp.Popen(['grep', '-n', 'pdb', 'jt.pir'], stdout=subp.PIPE).communicate()[0].split('\n')
    lines = []
    labels = []
    for l in entry_data:
        lines.append(l.split(':')[0])
        labels.append(l.split('|')[3])
    return lines, labels

def gen_ali(struc, nl):
    algmt = subp.Popen(['grep', '-A', nl, struc, 'jt.pir'], stdout=subp.PIPE).communicate()[0].split('\n')
    algmt.insert(1, 'structure:%s.pdb::::::::' % struc)
    f = open('%s.ali' % struc, 'w')
    for info in algmt[:-1]:
        f.write('%s\n' % info)
    return f

def main():
    seq = args.s.readlines()
    ln, pdbs = structures(args.a)
    for i, p in enumerate(pdbs):
        after = ln[i+1] - ln[i] - 1
        ali = gen_ali(p, after)
        with open(ali, 'a') as f:
            f.write(seq)

if __name__ == '__main__':
    main()
