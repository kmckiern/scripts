#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

data = sys.argv[1:]

def plot_pmf(element, timeseries):
    symbol = element['symbol']
    nb = element['num_bins']
    col = element['data_column']

    plt.clf()

    bins = np.arange(0, nb)
    for i, ts in enumerate(timeseries):
        t = np.genfromtxt(ts)
        hist, bin_edges = np.histogram(t[:,col], bins=bins, normed=1, label=data[i])
        plt.plot(bin_edges[:-1], hist, '-o', label=i)

    plt.legend(loc=0)
    plt.title('selectivity filter ' + symbol + ' PMF')
    plt.xlabel('n_' + symbol)
    plt.ylabel('P')
    plt.savefig('compare/' + symbol + '_hist.png')

def main():
    water_specs = {
        'symbol': 'h2o',
        'num_bins': 40,
        'data_column': 1
    }
    k_specs = {
        'symbol': 'k',
        'num_bins': 6,
        'data_column': 2
    }

    elements = [water_specs, k_specs]
    for j in elements:
        plot_pmf(j, data)

if __name__ == '__main__':
    main()
