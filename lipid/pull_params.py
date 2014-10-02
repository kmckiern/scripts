#!/usr/bin/env python
"""
this script pulls all relevant parameter specifications from a forcefield directory in order to generate a ForceBalance .itp file

parser makes lots of assumptions ... will make more general with time
"""
import argparse
import os
import subprocess as subp
import glob

parser = argparse.ArgumentParser(description='parse ff params from atom list file to out file')
parser.add_argument('--af', type=str, help='file of atoms in system')
parser.add_argument('--ff_dir', type=str, help='path to ff parameter directory')
parser.add_argument('--out', type=str, help='name of output parameter file')
opts = parser.parse_args()

def main():
    # get list of atoms for which we would like to find parameters
    atms = []
    for line in open(opts.af, 'r').readlines():
        atms.append(line.strip())
    # remove duplicates, and determine now many unique elements
    uniq_atms = set(atms)
    num_atms = len(atms)
    
    flz = glob.glob(opts.ff_dir + '*')
    of = opts.out

    # figure out with files have needed data
    p_cats = ['defaults', 'atomtypes', 'bondtypes', 'angletypes', 'dihedraltypes', 'pairtypes']
    cat_to_file = {}
    for p_type in p_cats:
        data = None
        for f in flz:
            data = subp.Popen(['grep', p_type, f], stdout=subp.PIPE).communicate()[0]
            # write relevant info to file
            if data is not None:
                cat_to_file[f] = p_type
                continue

    print uniq_atms

    # parse relevant files and record data to new file
    for fl in cat_to_file.keys():
        record = False
        with open(fl, 'r') as ref:
            lines = ref.readlines()
        with open(of, 'a') as out_file:
            for line in lines:
                # skip blank lines
                if not line.strip():
                    continue
                # look for parameter flags
                if line.startswith('['):
                    if any(key in line for key in p_cats):
                        out_file.write(line)
                        record = True
                    continue
                if record:
                    # stop recording when section ends (blank line)
                    if not line.strip():
                        out_file.write(line)
                        record = False
                        continue
                    elif any(atm in line.split() for atm in uniq_atms): 
                        out_file.write(line)
                        continue

if __name__ == '__main__':
    main()
