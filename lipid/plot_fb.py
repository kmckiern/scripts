#!/usr/bin/env python

"""
This script parses a forcebalance output file and plots the data.

example usage:
    >> ./amber_equil.py --root_dir path/to/root/dir --ligand oxy --rst --e1
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser(description='Parse and plot ForceBalance output data.')
parser.add_argument('--of', type=str, help='ForceBalance out file')
parser.add_argument('--ir', type=str, help='Range of iterations to plot.  enter as follows: first,last')
parser.add_argument('--pdf_out', type=str, help='PDF out file name')
# parser.add_argument('--nt', type=int, help='Number of condensed phase temperature points')
opts = parser.parse_args()

script_root = '/'.join(os.path.realpath(__file__).split('/')[0:-1])

# currently don't feel like getting this information in a general way.
properties = {'LipidDPPC Average Area per Lipid': 1, 'LipidDPPC Deuterium Order Parameter': 28, 'LipidDPPC Average Volume per Lipid': 1}
scd_key = properties.keys()[1]
scd_val = properties[scd_key]
property_files = script_root + '/possible_props.dat'

def parse_out(fb_out, ref_dict, data_dict, tp, iter_ls, selected_props):
    # read and record data
    i0 = iter_ls[0]
    record = False
    scd = False
    pi = None
    iteration = 0
    ln = 0
    boo = 0
    for line in open(fb_out, 'r'):
        # skip blank lines
        if not line.strip():
            continue
        # parse non-blank
        l = line.split()
        if record:
            # reset a bunch of variables
            if line.startswith('-'):
                record = False
                pi = None
                ln = 0
                boo += 1      
                if boo == len(selected_props):
                    boo = 0
                    iteration += 1
                continue
            elif l[0] == '#|':
                continue
            # deuterium requires special handling
            elif scd:
                if l[0] in tp:
                    scd_temp = float(l[0])
                    t = tp.index(l[0])
                    continue
                else:
                    data_dict[pi][iteration - i0][t][ln][1:] = float(l[1]), float(l[3])
                    if iteration == i0:
                        ref_dict[pi][t][ln][1:] = float(l[0])
                    if ln == properties[pi] - 1:
                        if scd_temp == float(tp[-1]):
                            scd = False
                            continue
                        else:
                            ln = 0
                    else:
                        ln += 1
            elif l[0] in tp:
                ln = tp.index(l[0])
                if iteration in iter_ls:
                    # s_cd is formatted really weirdly compared to the other properties
                    if pi not in scd_key:
                        data_dict[pi][iteration - i0][0][ln] = float(l[0]), float(l[4]), float(l[6])
                        if iteration == i0:
                            # don't record the exp data with zero weight
                            if abs(float(l[8])) < 1e-9:
                                ref_dict[pi][0][ln] = float(l[0]), np.NAN
                            else:
                                ref_dict[pi][0][ln] = float(l[0]), float(l[3])
                    else:
                        scd = True
                        scd_temp = float(l[0])
                        t = tp.index(l[0])
                        continue
        else:
            # look for property flag
            if l[0] == '#|':
                for prop in selected_props:
                    # if property of interest is read, begin recording data
                    if prop in line:
                        record = True
                        pi = prop
                        continue
            else:
                continue
    return data_dict, ref_dict

def gen_plots(ID, data_arrs, ref_arrs, props, iter_ls, tp):
    '''
    ID: for naming the PDF
    data_arrays: time series matrices for all properties
    ref_arr: experiemental data
    props: list of properties we wish to plot
    iter_ls: list of iterations to plot
    '''
    # Pull trj time information.
    pdfs = PdfPages(ID + '.pdf')
    n_colors = len(plt.rcParams['axes.color_cycle'])
    for p in props:
        plt.ylabel(p)
        plt.grid(True)
        # plot s_cd by carbon node
        if scd_key in p:
            """
            plt.figure(1)
            # iterate over iterations...
            for itr, arr in enumerate(data_arrs[p]):
                # nodes
                for cn, ds in enumerate(arr):
                    # temperature points
                    for t, temp in enumerate(tp):
                        plt.subplot(22t)
                        plt.title(temp)
                        plt.xlabel('hc node')
                        
                        
                    plt.plot(cn, ds[]
                    c_nodes = 
            """
            continue
        else:
            plt.xlabel('T / K')
            for itr, arr in enumerate(data_arrs[p]):
                data_label = 'iteration ' + str(iter_ls[itr])
                temp = arr[0][:,0]
                vals = arr[0][:,1]
                plt.plot(temp, vals, 'o', alpha = .75, label=data_label)
            plt.legend(prop={'size':6}, loc=4)
            plt.plot(temp, ref_arrs[p][0][:,1], 'ks', label='experiment')
        pdfs.savefig()
        plt.clf()
    pdfs.close()

def main():
    # get properties of interest from CL
    choices = raw_input('Which properties would you like to plot?\n' + open(property_files,'r').read() + '(separate choices using one comma)\n')
    inp = choices.split(',')
    # map user input to corresponding parsing strings
    props = [properties.keys()[int(i)] for i in inp]

    # figure out which iterations to care about
    i_init, i_final = opts.ir.split(',')
    if i_init == i_final:
        iters = [int(i_init)]
    else:
        iters = [i for i in range(int(i_init), int(i_final) + 1)]
    # this is temporary.  don't feel like being general
    temps = ['323.15', '333.15', '338.15', '353.15']
    
    c_nodes = np.arange(scd_val)

    """
    build data dictionary
        keys are referenced by property
        values are data matrices
    rows: temperature
    columns: t, val, stdev
    depth: vals at a given temp.  only useful for s_cd

    datums: dictionary, indexed by property, of calcuated property arrays from fb iterations
    ref_data: property dict of reference data arrays
    """
    nt = len(temps)
    datums = {}
    ref_data = {}
    for p in props:
        for i in iters:
            if not p in datums:
                if scd_key in p:
                    datums[p] = [np.zeros((nt, properties[p], 3))]
                    ref_data[p] = np.zeros((nt, properties[p], 2))
                else:
                    datums[p] = [np.zeros((properties[p], nt, 3))]
                    ref_data[p] = np.zeros((properties[p], nt, 2))
            else:
                if scd_key in p:
                    datums[p].append(np.zeros((nt, properties[p], 3)))
                else:
                    datums[p].append(np.zeros((properties[p], nt, 3)))

    if scd_key in datums:
        for itr in datums[scd_key]:
            for temp in itr:
                temp[:,0] = c_nodes

    # parse forcebalance output file.
    data, ref = parse_out(opts.of, ref_data, datums, temps, iters, props)

    # gen_plots(opts.pdf_out, datums, ref_data, props, temps, iters)

if __name__ == '__main__':
    main()
