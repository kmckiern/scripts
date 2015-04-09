#!/usr/bin/env python

from optparse import OptionParser
import numpy as np

parser = OptionParser()
parser.add_option('--ts', help='Input time series')
parser.add_option('--nb', help='Max number of blocks', type=int)
(opts, args) = parser.parse_args()

def divide_ts(ts, n_blocks, len_block):
    block_avgs = np.zeros(n_blocks)

    for b in range(n_blocks):
        block_avgs[b] = np.mean(ts[(b*len_block):((b+1)*len_block)])

    return np.std(block_avgs)

def main():
    df = opts.ts
    ts = np.genfromtxt(df)
    file_pref = df.split('.')[0]

    n_dp = ts.shape[0]
    n_b = opts.nb

    sigmas = np.zeros((n_b + 1, 2))

    for block in range(1, n_b + 1):
        block_size = n_dp / block
        sigmas[block, 0] = block
        sigmas[block, 1] = divide_ts(ts[:,1], block, block_size)

    np.savetxt(file_pref + '_block_avg_' + str(n_b) + '.dat', sigmas)

if __name__ == '__main__':
    main()
