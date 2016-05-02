#!/bin/env python

import pandas as pd
import argparse
import os
import numpy as np
import IPython
from subprocess import Popen, PIPE, STDOUT
import multiprocessing

def call_cl(command_lst, std=PIPE, pipe_args=[]):
    p = Popen(command_lst, stdout=PIPE, stdin=PIPE, stderr=STDOUT, shell=True)
    out, err = p.communicate(input=b'\n'.join(pipe_args))
    return out, err

parser = argparse.ArgumentParser(description='generate tica trj')
parser.add_argument('--df', type=str, help='micro state id')
parser.add_argument('--od', type=str, help='output dir')
args = parser.parse_args()

gmx_pref = '/home/server/opt/gromacs-5.0.5/bin/gmx '
trjconv = gmx_pref + 'trjconv'

top_root = '/home/server/clc/tops/verbose'
gmxt_root = '/home/server/clc/tops/gmx_tops'

def raw_trjs(trj_in, trj_out, top, frame):
    run_trjcnv = trjconv + ' -f ' + trj_in + ' -o ' + trj_out + ' -s ' + top \
        + ' -b ' + str(frame) + ' -e ' + str(frame) + ' -tu ps -pbc nojump'
    print (trj_out)
    call_cl('echo -ne \"0\\n\" | ' + run_trjcnv)

def wrap_raw_trjs(df_vals):
    p, tr, trj_in, frame, trj_out = df_vals

    if tr % 2 == 0:
        run = '0'
    else:
        run = '1'
    top = os.path.join(top_root, 'p' + str(p) + '_r' + run + '.pdb')

    raw_trjs(trj_in, trj_out, top, frame)

    return None

def main():
    df = pd.read_pickle(args.df)
    
    projs = np.array(df['PROJ'])
    runs = np.array(df['RUN'])
    raw_trjs = np.array(df['raw_trj'])
    frames = np.array([i*200.0 for i in df['fn']])

    nf = df.shape[0]
    out = [os.path.join(args.od, str(i) + '.pdb') for i in range(nf)]
    
    input_data = []
    for f in range(nf):
        input_data.append([i[f] for i in [projs, runs, raw_trjs, frames, out]])

    p = multiprocessing.Pool(14)
    p.map(wrap_raw_trjs, input_data)
    p.close()

if __name__ == '__main__':
    main()
