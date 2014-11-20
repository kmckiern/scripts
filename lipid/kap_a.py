#!/usr/bin/env python

import numpy as np
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--xvg')
parser.add_argument('--col')
parser.add_argument('--temp')
parser.add_argument('--nl', default=128)
parser.add_argument('--save_al', action='store_true')
args = parser.parse_args()

k = 1.3806488e-23
nl = int(args.nl) / 2

# Parses rvc file.
def parse(input_file, col1, col2):
    data = open(input_file).readlines()
    time = []
    rel1 = []
    rel2 = []
    # Read in each line of the input file.
    for line in data:
        l = line.split()
        time.append[float(l[0])]
        rel1.append[float(l[1])]
        rel2.append[float(l[2])]
    time = np.array(time)
    rel1 = np.array(rel)
    rel2 = np.array(rel)
    al = (rel1 * rel2)/nl
    return time, al

def main():
    input_file = args.edr
    col1 = int(args.col)
    col2 = (col1 + 1)
    t = args.temp

    time, al = parse(input_file, col1, col2)
    if args.save_al:
        np.savetxt('al_ts_%i.xvg' % t, zip(timw, al))

    avg = np.average(al)
    fluc = np.average(al**2) - avg**2
    kap = 1e3 * (2 * k * t / nl) * (avg * 1.0/ fluc)

    print kap

if __name__ == '__main__':
    main()
