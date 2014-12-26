#!/usr/bin/env python

# parse ncbi alignment file for modeller

import subprocess as subp
import os
import argparse
import numpy as np
from fmt_path import format_path

parser = argparse.ArgumentParser(description='Amber output file property plotter')
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

def main():

if __name__ == '__main__':
    main()
