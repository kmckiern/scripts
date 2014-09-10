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
parser.add_argument('--out_dir', type=str, help='Amber output file directory')
parser.add_argument('--pdf_name', type=str, default='time_series', help='Name of final pdf file')
parser.add_argument('--legend', action='store_true', help='Show legend on plots')
parser.add_argument('--final', action='store_true', help='Highlight time series final point')
args = parser.parse_args()

# Generalize this to represent a specified subset.
property_ndxs = {'ANGLE': 8, 'BOND': 7, 'DIHED': 9, 'Density': 20, 'EEL': 11, 'EELEC': 13, 'EHBOND': 14, 'EKCMT': 16, 'EKtot': 5, 'EPtot': 6, 'Etot': 4, 'NB': 10, 'NSTEP': 0, 'PRESS': 3, 'RESTRAINT': 15, 'SURFTEN': 19, 'TEMP(K)': 2, 'TIME(PS)': 1, 'VDWAALS': 12, 'VIRIAL': 17, 'VOLUME': 18}

def parse_out(of, props):
    """
    of: mdcrd.out file for parsing.
    props: system properties of interest.
    """
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
                time_ax = np.array(arr[:, props['TIME(PS)']])
                data_label = file_labels[ndx].split('/')[-1].split('.')[0]
                plt.plot(time_ax, arr[:, props[p]], alpha=.66, label=data_label)
                if args.final:
                    plt.plot(time_ax[-1], arr[:, props[p]][-1], 'o', ms=12, color=plt.rcParams['axes.color_cycle'][ndx % n_colors])
            if args.legend:
                plt.legend(prop={'size':6})
            pdfs.savefig()
            plt.clf()
    pdfs.close()

def main():
    out_data = args.out_dir
    out_files = [f for f in os.listdir(out_data) if '.out' in f]

    # Parse amber output file, return data matrix.
    data_mtrxs = []
    for mdcrd in out_files:
        data_mtrxs.append(parse_out(out_data + '/' + mdcrd, property_ndxs))

    # Plot 
    write_pref = out_data + args.pdf_name
    gen_plots(write_pref, data_mtrxs, property_ndxs, out_files) 

if __name__ == '__main__':
    main()
