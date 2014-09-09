#!/bin/env python

"""
Reads an AMBER mdcrd.out file, and plots specified properties over a specified range.

For now, it just plot all of the property data, over all time, from the mdcrd.out file.
In the future, I'd like to add the ability to:
    - take CL input for a subset of the properties to be plotted
        data_arr does not need to be built fully
    - take in arguments for the frame range.
    - make the parser better.  building general parsers seems hard in general as the files are fairly unique.
    - cut out os calls.
        don't need temp file.  just have python read in file?
"""

from subprocess import call
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser(description='Amber output file property plotter')
parser.add_argument('--out_file', type=str, help='Amber output file')
args = parser.parse_args()

# Global variables?
temp_of = 'temp_props.dat'
property_ndxs = {'ANGLE': 8, 'BOND': 7, 'DIHED': 9, 'Density': 20, 'EEL': 11, 'EELEC': 13, 'EHBOND': 14, 'EKCMT': 16, 'EKtot': 5, 'EPtot': 6, 'Etot': 4, 'NB': 10, 'NSTEP': 0, 'PRESS': 3, 'RESTRAINT': 15, 'SURFTEN': 19, 'TEMP(K)': 2, 'TIME(PS)': 1, 'VDWAALS': 12, 'VIRIAL': 17, 'VOLUME': 18}

def parse_out(of, props):
    # parse amber output file, write data file of step/property information
    subprocess.call('grep -A 9 '========' %s.out >> %s.dat' % (of, temp_of), shell=True)
    # initialize
    num_frames = int(subprocess.call('grep -c TIME %s.dat' % temp_of, shell=True))
    data_arr = np.zeros((num_frames, len(props)))
    frame = 0
    # populate data_arr
    for line in open(temp_of):
        if line.startswith('====='):
            if not line.startswith('--'):
                if len(line) > 1:
                    l = line.split()
                    for ent in l:
                        if l in props:
                            data_arr[frame, props[ent]]= float(l[ent + 2])
            else:
                frame += 1
    return data_arr

def plot(ID, data_arr, props):
    '''
    props: list of properties strings we wish to plot.
        index is mapped?
    time_ax: linear array of time points over specified range.
    ts: list of numpy time serieses 
    '''
    # Pull trj time information.
    time_ax = np.array(data_arr[:, props['TIME(PS)']])
    pdfs = PdfPages('%s_ts.pdf' % ID)
    for ndx, p in enumerate(props):
        plt.title('%s --- %s(t)' % (ID, p))
        plt.ylabel(p)
        plt.xlabel('t / ps')
        plt.plot(time_ax, ts[ndx])
        pp.savefig()
    pp.close()

main():
    out_data = args.out_file
    property_list = property_ndxs.keys()

    # Parse amber output file, return data matrix.
    data_mtrx = parse_out(out_data, property_list) 
    # Format 
    write_pref = out_data.split('.')[0]
    plot(write_pref, data_mtrx, property_list) 

if __name__ == '__main__':
    main()
