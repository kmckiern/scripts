#!/bin/env python

import argparse
import glob
import numpy as np
from collections import Counter
import IPython

parser = argparse.ArgumentParser(description='get water and ion flux of trek sf')
parser.add_argument('--ad', type=str, help='array directory', 
    default='/home/server/trek/flux/')
parser.add_argument('--opref', type=str, help='path to output directory', 
    default='/home/server/trek/flux/')
args = parser.parse_args()

def main():
    # just worry about potassium for now
    files = glob.glob(args.ad + 'wa*npy')

    for f in files: 
        label = f.split('/')[-1].split('.npy')[0]

        trj = np.load(f)

        # break label symmetry for faster calculation
        trj[trj==-1] = -2

        dtrj = np.diff(trj, axis=0)
        net_flux = np.zeros(1+dtrj.shape[0], dtype=int)
        for ndx, frame in enumerate(dtrj):
            counts = Counter(frame)
            flux_in = counts[2.0] + counts[-1.0]
            flux_out = counts[-2.0] + counts[1.0]
            net_flux[ndx+1] = int(flux_in - flux_out)

        np.savetxt(args.opref + label, net_flux, fmt='%i')

if __name__ == '__main__':
    main()
