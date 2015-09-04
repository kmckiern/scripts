#!/bin/env python

"""
script for reading in text file(s) and returning them as lists

i am bored on a plane
"""

import argparse

parser = argparse.ArgumentParser(description='text binning')
parser.add_argument('--f', type=str, help='file(s) to bin')
args = parser.parse_args()

def main():
    f = args.f

    plural = f.split(',')
    if len(plural) > 1:
        for fl in plural:
            txt = None
            with open(fl, 'r') as rf:
                txt += rf.readlines()
    else:
        with opeen(f, 'r') as rf:
            txt = rf.readlines()

    return txt

if __name__ == '__main__':
    main()
