#!/usr/bin/env python

# parse ncbi alignment file for modeller

import subprocess as subp
import argparse
import sys
from toolz import call_cml

parser = argparse.ArgumentParser(description='parse ncbi ms alignment file for modeller')
parser.add_argument('--a', type=str, help='ncbi alignment file')
parser.add_argument('--s', type=str, help='sequence file')
args = parser.parse_args()

def get_stuff(call):
     return call_cml(call).split('\n')[:-1]

# this prob isn't pythonic but w/e.
def structures(af):
    entry_data = get_stuff('grep -n pdb %s' % af) 
    lines = []
    labels = []
    for l in entry_data:
        lines.append(int(l.split(':')[0]))
        labels.append(l.split('|')[3])
    last = get_stuff('wc %s' % af)
    lines.append(int(last[0].split()[0])+1)
    return lines, labels

def gen_ali(af, struc, nl):
    algmt = get_stuff('grep -A %s %s %s' % (nl, struc, af)) 
    algmt.insert(1, 'structure:%s.pdb::::::::' % struc)
    f = open('%s.ali' % struc, 'w')
    for info in algmt:
        f.write('%s\n' % info)
    return f

def main():
    seq = open(args.s, 'r').read()
    ln, pdbs = structures(args.a)
    for i, p in enumerate(pdbs):
        after = str(ln[i+1] - ln[i] - 1)
        ali = gen_ali(args.a, p, after)
        ali.write(seq)

if __name__ == '__main__':
    main()
