#!/bin/env python

import os
import pandas as pd
import IPython
import sys
import numpy as np
import argparse

sr = '/home/kmckiern/scripts/'
sys.path.insert(0, os.path.join(sr, 'py_general'))
from toolz import call_cl

parser = argparse.ArgumentParser(description='pull simulation')
parser.add_argument('--state', type=str, help='micro state id')
parser.add_argument('--gnum', type=int, help='gro file number')
args = parser.parse_args()

state = args.state
gnum = args.gnum

gmx_pref = '/home/kmckiern/opt/gromacs-5.0.5/bin/gmx '
make_ndx = gmx_pref + 'make_ndx'
grompp = gmx_pref + 'grompp'
mdrun = gmx_pref + 'mdrun'
trjconv = gmx_pref + 'trjconv'

pref_pref = '/home/kmckiern/sims/trek/'
ic_pref = os.path.join(pref_pref, 'ICs')
pmf_pref = os.path.join(pref_pref, 'pull')
top_pref = os.path.join(pref_pref, 'tops')

pull_mdp = os.path.join(pmf_pref, 'pull.mdp')
push_mdp = os.path.join(pmf_pref, 'push.mdp')

# get data from df
df = os.path.join(pmf_pref, 'gmxing', state, 'meta.p')
data = pd.read_pickle(df)
rele = data.loc[gnum]
gro = rele['new']
top = rele['top']
top = os.path.join(top_pref, top + '.top')
indices = np.array(rele[['b0','b1','b2','b3','bp1']], dtype=int)

# add 1 for indexing difference
def make_ion_ndx(indices, gro, hold=False):
    ndx_file = gro.replace('.gro', '.ndx')
    # ion_top = ' '.join([str(indices[i]+1) for i in ion_indices])
    # ipa = ['a ' + ion_top, 'name 13 ion_top']
    ipa = ['name 13 ion_pull']
    ipa = ipa + ['q']
    run_make_ndx = 'echo \"' + ' \n'.join(ipa) + '\" | ' + make_ndx + \
            ' -f ' + gro + ' -o ' + ndx_file
    return run_make_ndx, ndx_file

def grompp_str(mdp_file, in_gro, out_tpr, top, ndx_file=None):
    run_grompp = grompp + ' -f ' + mdp_file + ' -c ' + in_gro + ' -p ' + top + \
        ' -o ' + out_tpr + '.tpr -maxwarn 1'
    if ndx_file != None:
        run_grompp = run_grompp + ' -n ' + ndx_file
    return run_grompp

def trjconv_str(top, trj, fn, of, cent=False):
    run_trjconv = trjconv + ' -f ' + trj + ' -s ' + top.replace('.top','.pdb') \
        + ' -o ' + of + ' -b ' + fn + ' -e ' + fn
    if cent:
        run_trjconv += ' -pbc mol -center'
    return run_trjconv

def pull_sim():
    # make ndx
    run_make_ndx, ndx_file = make_ion_ndx(indices, gro)
    print(run_make_ndx)
    call_cl(run_make_ndx)

    # grompp
    out_name = gro.replace('.gro', '_pull_1')
    run_grompp = grompp_str(pull_mdp, gro, out_name, top, ndx_file)
    print(run_grompp)
    call_cl(run_grompp)

    # mdrun
    run_mdrun = mdrun + ' -deffnm ' + out_name
    print(run_mdrun)
    call_cl(run_mdrun)

    return None
    
pull_sim()
