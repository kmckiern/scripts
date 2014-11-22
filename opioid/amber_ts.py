#!/usr/bin/env python

"""
Reads an AMBER mdcrd.out file, and plots specified properties over a specified range.
For now, it just plot all of the property data, over all time, from the mdcrd.out file.

In the future, I'd like to add the ability to:
    - take CL input for a subset of the properties to be plotted
        data_arr does not need to be built fully
    - take in arguments for the frame range.
    - cut out os calls.
        can't think of a better way to get number of frames without reading in file.
    - have final points print on top of ts data.
    - figure out what to do with the legend.
    - if job is terminated, some vals aren't recorded.
        exception handling.
"""

import subprocess as subp
import os
import argparse
import numpy as np
from fmt_path import format_path
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser(description='Amber output file property plotter')
parser.add_argument('--out_dir', type=str, help='Amber output file directory')
parser.add_argument('--pdf_name', type=str, default='time_series', help='Name of final pdf file')
parser.add_argument('--legend', action='store_true', help='Show legend on plots')
parser.add_argument('--final', action='store_true', help='Highlight time series final point')
parser.add_argument('--min', action='store_true', help='Specify if the out files are from a minimization trajectory')
args = parser.parse_args()

# Generalize this to represent a specified subset.
property_ndxs = {'ANGLE': 8, 'BOND': 7, 'DIHED': 9, 'Density': 20, 'EEL': 11, 'EELEC': 13, 'EHBOND': 14, 'EKCMT': 16, 'EKtot': 5, 'EPtot': 6, 'Etot': 4, 'NB': 10, 'NSTEP': 0, 'PRESS': 3, 'RESTRAINT': 15, 'SURFTEN': 19, 'TEMP(K)': 2, 'TIME(PS)': 1, 'VDWAALS': 12, 'VIRIAL': 17, 'VOLUME': 18, 'EAMBER (non-restraint)': 21}
min_ndxs = {'BOND': 0, 'ANGLE': 1, 'DIHED': 2, 'VDWAALS': 3, 'EEL': 4, 'HBOND': 5, '1-4 VDW': 6, '1-4 EEL': 7, 'RESTRAINT': 8, 'NSTEP': 9, 'ENERGY': 10, 'RMS': 11, 'GMAX': 12, 'NUMBER': 13}

def parse_out(of, props):
    """
    of: mdcrd.out file for parsing.
    props: system properties of interest.
    """
    # initialize
    if args.min:
        num_frames = int(subp.Popen(['grep', '-c', 'NSTEP', of], stdout=subp.PIPE).communicate()[0])
        num_frames -= 1
        blerg = False
    else:
        num_frames = int(subp.Popen(['grep', '-c', 'TIME', of], stdout=subp.PIPE).communicate()[0])
        summary_frames = int(subp.Popen(['grep', '-c', 'A V E R A G E S', of], stdout=subp.PIPE).communicate()[0])
        summary_frames += int(subp.Popen(['grep', '-c', 'F L U C T U A T I O N S', of], stdout=subp.PIPE).communicate()[0])
        num_frames -= summary_frames
    data_arr = np.empty((num_frames, len(props)))
    data_arr[:] = np.nan
    frame = 0
    record = False
    # populate data_arr
    for line in open(of, 'r'):
        # skip blank lines
        if not line.strip():
            continue
        line = line.replace('=-', '= -')
        l = line.split()
        if record and frame < num_frames: 
            if args.min:
                if l[0] == 'NSTEP':
                    blerg = True
                    min_ndx = []
                    min_ent = []
                    for ndx, ent in enumerate(l):
                        if ent in props:
                            min_ndx.append(ndx)
                            min_ent.append(ent)
                    continue
                if blerg:
                    for ndx, ent in enumerate(min_ent):
                        data_arr[frame, props[ent]] = float(l[min_ndx[ndx]])
                    blerg = False
                    frame += 1
                elif l[0] in props:
                    for ndx, ent in enumerate(l):
                        if ent in props:
                            data_arr[frame, props[ent]] = float(l[ndx + 2])
            else:
                if l[0] in props:
                    for ndx, ent in enumerate(l):
                        if ent in props:
                            data_arr[frame, props[ent]] = float(l[ndx + 2])
                elif l[0].startswith('-------'):
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
                if args.min:
                    time_ax = np.array(arr[:, props['NSTEP']])
                else:
                    time_ax = np.array(arr[:, props['TIME(PS)']])
                data_label = file_labels[ndx].split('/')[-1].split('.')[0]
                # if 'equil' in data_label:
                #    plt.plot(time_ax[:300], arr[:, props[p]][:300], alpha=.5, label=data_label)
                # else:
                plt.plot(time_ax[:-2], arr[:, props[p]][:-2], alpha=.5, label=data_label)
                if args.final:
                    plt.plot(time_ax[-1], arr[:, props[p]][-1], 'o', ms=12, color=plt.rcParams['axes.color_cycle'][ndx % n_colors])
            if args.legend:
                plt.legend(prop={'size':6}, loc=4)
            pdfs.savefig()
            plt.clf()
    pdfs.close()

def main():
    out_data = format_path(args.out_dir)
    out_files = [f for f in os.listdir(out_data) if '.out' in f]

    if args.min:
        prop_arr = min_ndxs
    else:
        prop_arr = property_ndxs

    # Parse amber output file, return data matrix.
    data_mtrxs = []
    for mdcrd in out_files:
        data_mtrxs.append(parse_out(out_data + '/' + mdcrd, prop_arr))

    # Plot 
    write_pref = out_data + args.pdf_name
    gen_plots(write_pref, data_mtrxs, prop_arr, out_files) 

if __name__ == '__main__':
    main()
