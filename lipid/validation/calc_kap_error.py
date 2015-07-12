#!/bin/env python

import argparse
import numpy as np
from forcebalance.nifty import mean_stderr
import IPython

parser = argparse.ArgumentParser(description='for calculating the error of a \
    nx2 time series')
parser.add_argument('--ts', type=str, help='time series')
args = parser.parse_args()

def main():
    data = np.genfromtxt(args.ts)
    print mean_stderr(data[:,1][::100])

if __name__ == '__main__':
	main()
