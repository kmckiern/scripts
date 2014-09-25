#!/usr/bin/env python

"""
This script parses a forcebalance output file and plots the data.

example usage:
    >> ./amber_equil.py --root_dir path/to/root/dir --ligand oxy --rst --e1
"""

import subprocess as subp
import os
import shutil
import argparse
import math
import copy
import numpy as np

parser = argparse.ArgumentParser(description='Parse and plot ForceBalance output data.')
parser.add_argument('--of', type=str, help='ForceBalance out file')
parser.add_argument('--ir', type=str, help='Range of iterations to plot.  enter as follows: first,last')
# parser.add_argument('--nt', type=int, help='Number of condensed phase temperature points')
opts = parser.parse_args()

script_root = '/'.join(os.path.realpath(__file__).split('/')[0:-1])

# currently don't feel like getting this information in a general way.
properties = {'LipidDPPC Average Area per Lipid': 1, 'LipidDPPC Deuterium Order Parameter': 28, 'LipidDPPC Average Volume per Lipid': 1}
property_files = script_root + '/possible_props.dat'

def parse_out(fb_out, ref_dict, data_dict, tp, iter_ls, selected_props):
    # read and record data
    i0 = iter_ls[0]
    record = False
    scd = False
    pi = None
    iteration = 0
    c_node = 0
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
                c_node = 0
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
                    ln = tp.index(l[0])
                    continue
                else:
                    data_dict[pi][iteration - i0][c_node][ln] = scd_temp, float(l[1]), float(l[3])
                    if iteration == i0:
                        ref_dict[pi][c_node][ln] = scd_temp, float(l[0])
                    if c_node == properties[pi] - 1:
                        if scd_temp == float(tp[-1]):
                            scd = False
                            continue
                        else:
                            c_node = 0
                    else:
                        c_node += 1
            elif l[0] in tp:
                ln = tp.index(l[0])
                if iteration in iter_ls:
                    # s_cd is formatted really weirdly compared to the other properties
                    if pi not in properties.keys()[1]:
                        data_dict[pi][iteration - i0][0][ln] = float(l[0]), float(l[4]), float(l[6])
                        if iteration == i0:
                            ref_dict[pi][0][ln] = float(l[0]), float(l[3])
                    else:
                        scd = True
                        scd_temp = float(l[0])
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

def gen_plots(ID, data_arrays, props, file_labels):
    '''
    ID: for naming the PDF.
    data_arrays: time series matrix for all files.
    props: list of properties we wish to plot.
    file_labels: file names (for legend).
    '''
    # Pull trj time information.
    pdfs = PdfPages(ID + '.pdf')
    for p in props:
        if p in ('TIME(PS)', 'NSTEP', 'EHBOND'):
            continue
        else:
            n_colors = len(plt.rcParams['axes.color_cycle'])
            plt.ylabel(p)
            plt.xlabel('t / ps')
            plt.grid(True)
            for ndx, arr in enumerate(data_arrays):
                if args.min:
                    time_ax = np.array(arr[:, props['NSTEP']])
                else:
                    time_ax = np.array(arr[:, props['TIME(PS)']])
                data_label = file_labels[ndx].split('/')[-1].split('.')[0]
                # if 'equil' in data_label:
                #    plt.plot(time_ax[:300], arr[:, props[p]][:300], alpha=.5, label=data_label)
                # else:
                plt.plot(time_ax, arr[:, props[p]], alpha=.5, label=data_label)
                if args.final:
                    plt.plot(time_ax[-1], arr[:, props[p]][-1], 'o', ms=12, color=plt.rcParams['axes.color_cycle'][ndx % n_colors])
            if args.legend:
                plt.legend(prop={'size':6}, loc=4)
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
                datums[p] = [np.zeros((properties[p], nt, 3))]
                ref_data[p] = np.zeros((properties[p], nt, 2))
            else:
                datums[p].append(np.zeros((properties[p], nt, 3)))

    # parse forcebalance output file.
    data, ref = parse_out(opts.of, ref_data, datums, temps, iters, props)

    # write_pref = out_data + args.pdf_name
    # gen_plots(write_pref, data_mtrxs, prop_arr, out_files)

if __name__ == '__main__':
    main()
