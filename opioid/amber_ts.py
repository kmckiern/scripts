#!/usr/bin/env python

"""
Reads an AMBER mdcrd.out file, and plots specified properties over a specified range.

For now, it just plot all of the property data, over all time, from the mdcrd.out file.
In the future, I'd like to add the ability to:
    - compare different out files on same plots
    - take CL input for a subset of the properties to be plotted
        data_arr does not need to be built fully
    - take in arguments for the frame range.
    - cut out os calls.
        can't think of a better way to get number of frames without reading in file.
"""

import subprocess as subp
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser(description='Amber output file property plotter')
parser.add_argument('--out_file', type=str, help='Amber output file')
args = parser.parse_args()

# Generalize this to represent a specified subset.
property_ndxs = {'ANGLE': 8, 'BOND': 7, 'DIHED': 9, 'Density': 20, 'EEL': 11, 'EELEC': 13, 'EHBOND': 14, 'EKCMT': 16, 'EKtot': 5, 'EPtot': 6, 'Etot': 4, 'NB': 10, 'NSTEP': 0, 'PRESS': 3, 'RESTRAINT': 15, 'SURFTEN': 19, 'TEMP(K)': 2, 'TIME(PS)': 1, 'VDWAALS': 12, 'VIRIAL': 17, 'VOLUME': 18}

def parse_out(of, props):
    # initialize
    num_frames = int(subp.Popen(['grep', '-c', 'TIME', of], stdout=subp.PIPE).communicate()[0])
    data_arr = np.zeros((num_frames, len(props)))
    frame = 0
    record = False
    # populate data_arr
    for line in open(of, 'r'):
        # skip blank lines
        if not line.strip():
            continue
        l = line.split()
        if record:
            if l[0] in props:
                for ndx, ent in enumerate(l):
                    if ent in props:
                        data_arr[frame, property_ndxs[ent]] = float(l[ndx + 2])
            elif l[0].startswith('======='):
                frame += 1
        # skip over everything before results section
        elif l[-1] == 'RESULTS':
            record = True
    return data_arr

def gen_plots(ID, data_arr, props):
    '''
    props: list of properties strings we wish to plot.
        index is mapped?
    time_ax: linear array of time points over specified range.
    ts: list of numpy time serieses 
    '''
    # Pull trj time information.
    time_ax = np.array(data_arr[:, props['TIME(PS)']])
    pdfs = PdfPages('%s_ts.pdf' % ID)
    for p in props:
        plt.ylabel(p)
        plt.xlabel('t / ps')
        plt.grid(True)
        plt.plot(time_ax, data_arr[:, props[p]])
        pdfs.savefig()
        plt.clf()
    pdfs.close()

def main():
    out_data = args.out_file
    property_list = property_ndxs.keys()

    # Parse amber output file, return data matrix.
    data_mtrx = parse_out(out_data, property_ndxs)
    # np.savetxt('dm.dat', data_mtrx, fmt='%1.4f')

    # Plot 
    write_pref = out_data.split('.')[0]
    gen_plots(write_pref, data_mtrx, property_ndxs) 

if __name__ == '__main__':
    main()
