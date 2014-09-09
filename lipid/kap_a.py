#!/usr/bin/env python

import numpy as np
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--edr')
parser.add_argument('--col')
parser.add_argument('--temp')
args = parser.parse_args()

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
	al = (rel1 * rel2)/64
	return time, al

def main():
	input_file = args.edr
	col1 = int(args.col)
	col2 = (col1 + 1)
	t = args.temp
	k = 1.3806e-23
	nl = 64

	time, al = parse(input_file, col1, col2)

	avg = np.average(al)
	fluc = np.average(al**2) - avg**2

	kap = (2 * k * t * avg) / (nl * fluc * 1.0)

	print kap

if __name__ == '__main__':
	main()
