#!/usr/bin/env python

"""
feed in bilayer trj and generate output file with thermodynamic and structural characterization
example usage:
    >> python /path/to/script/dir/analyze_bilayer.py --trj /path/to/trj.trr --o meh.dat --xyz tp/
"""

import argparse
from subprocess import Popen, PIPE, STDOUT

parser = argparse.ArgumentParser(description='full analysis of bilayer trajectories')
parser.add_argument('--prod_dir', type=str, help='path to directory with trj and tpr')
parser.add_argument('--prod_pref', type=str, help='file preface for gmx files')
parser.add_argument('--valid_dir', type=str, help='name of directory for gmx output xyz files')
parser.add_argument('--o', type=str, help='name of output file')
args = parser.parse_args()

def cl_gmx(command_lst, pipe_args):
    p = Popen(command_lst, stdout=PIPE, stdin=PIPE, stderr=STDOUT)  
    p.communicate(input='\n'.join(pipe_args))

"""def parse_gmx(prop):

# area per lipid
def get_al():

# deuterium order parameter
def get_al():

# volume per lipid
def get_vl():

# isothermal compressibility modulus
def get_al():

# x-ray structure factor
def get_fq(): 

# diffusion constant
def get_dl():"""

def get_struc(pdir, ppref, vdir):
    props = ['Volume', 'Box-X', 'Box-Y']
    ge = ['g_energy', '-f', pdir + ppref + '.edr', '-xvg', 'no', '-o', vdir + 'struc.xvg']
    cl_gmx(ge, props)

def main():
    pdir = args.prod_dir
    ppref = args.prod_pref
    vdir = args.valid_dir
    out = args.o

    get_struc(pdir, ppref, vdir)

if __name__ == '__main__':
    main()
