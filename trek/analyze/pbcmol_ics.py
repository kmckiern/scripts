#!/bin/env python

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description='get water and ion occupancy of trek sf')
parser.add_argument('--df',type=str,help='trajectory file')
args = parser.parse_args()

import sys
sr = '/home/kmckiern/scripts/'
sys.path.insert(0, os.path.join(sr, 'py_general'))
from toolz import call_cl

top_pref = '/home/kmckiern/sims/trek/pmf/gmx/top/tops-p9761'
gmx_pref = '/home/kmckiern/opt/gromacs/bin'

def main():
    df = pd.read_pickle(args.df)
    dirs = df.index
    for state in dirs:
        frames = df.loc[state]
        for ndx, IC in enumerate(frames):
            if IC == None:
                continue
            else:
                print (state, IC)
                # prep everything
                odir = os.path.join('done', 'ICs', state)
                if not os.path.isdir(odir):
                    os.makedirs(odir)
                fr, top = IC
                trace = fr.split('.xtc')[0]
                topology = os.path.join(top_pref, top + '.top')
                trajectory = os.path.join('ICs', state, str(ndx) + '.gro')
                t_out = os.path.join(odir, str(ndx) + '.gro')
                # run grompp
                g_out = os.path.join('tpr', trace + '_' + str(ndx) + '.tpr')
                run_grompp = os.path.join(gmx_pref, 'grompp') + ' -f files/npt.mdp -c ' \
                    + trajectory + ' -p ' + topology + ' -o ' +  g_out + ' -maxwarn 1'
                call_cl(run_grompp)
                # run trjconv
                run_trjconv = os.path.join(gmx_pref, 'trjconv') + ' -f ' + trajectory \
                    + ' -s ' + g_out + ' -o ' +  t_out + ' -pbc mol -center'
                call_cl(run_trjconv, pipe_args=[b'Protein', b'System'])

if __name__ == '__main__':
    main()
