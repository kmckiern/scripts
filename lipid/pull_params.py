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

    # sanity check
    print uniq_atms

    # figure out with files have needed data
    cats = {'defaults': 0, 'atomtypes': 1, 'bondtypes': 2, 'angletypes': 3, 'dihedraltypes': 4, 'pairtypes': 2}
    p_cats = cats.keys()
    cat_to_file = {}
    for p_type in p_cats:
        data = None
        for f in flz:
            data = subp.Popen(['grep', p_type, f], stdout=subp.PIPE).communicate()[0]
            # write relevant info to file
            if data is not None:
                cat_to_file[f] = p_type
                continue

    # parse relevant files and record data to new file
    for fl in cat_to_file.keys():
        record = False
        prop = None
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
                        prop = line.split()[1].strip()
                        record = True
                    continue
                if record:
                    # skip switches
                    if line.startswith('#'):
                        continue
                    # stop recording when section ends (blank line)
                    if not line.strip():
                        out_file.write(line)
                        record = False
                        continue
                    else:
                        l = line.split()
                        atms = set(l[0:cats[prop]])
                        # check if atms list is a subset of uniq_atm set
                        if atms & uniq_atms == atms:
                            out_file.write(line)
                            continue

if __name__ == '__main__':
    main()
