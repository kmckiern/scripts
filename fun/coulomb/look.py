#!/bin/env python

import glob
import natsort
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import IPython

sns.set_style(style='white')

central_gly = [101,210,362,471]
frame_rate = .16

f_pref = 'PROJ9712_RUN3_CLONE18_'
pklz = glob.glob('forces/' + f_pref + '*.pkl')
pklz = natsort.natsorted(pklz)
lpklz = [pd.read_pickle(i) for i in pklz]
num_pts = len(lpklz)

def main():
    forces = ['PeriodicTorsionForce','NonbondedForce']
    for force in forces:
        for a_type in ['K+','GLY']: # ['ASN','ASP','GLY','ILE','PHE','THR','TYR','VAL']:
            plt.clf()
            # k ions only have nontrivial values for the nonbonded force
            if ((a_type == 'K+') & (force != forces[-1])):
                continue
            for ndx, df in enumerate(lpklz):
                # remove dataframe duplicates (hack for now)
                df = df.reset_index().drop_duplicates(subset='index',
                        keep='last').set_index('index')
                # filter data frame
                if a_type == 'GLY':
                    if force == forces[0]:
                        atm = ['C']
                    if force == forces[1]:
                        atm = ['O']
                    meta = df.loc[(df['resName'] == a_type) &
                            (df['resSeq'].isin(central_gly)) &
                            (df['name'].isin(atm))]
                if a_type == 'K+':
                    meta = df.loc[(df['resName'] == a_type)]
                atoms = meta[force]
                # get magnitude of force for filtered df
                f_mag = np.array([np.linalg.norm(x) for x in atoms])
                if ndx == 0:
                    n_atom = len(f_mag)
                    data = np.zeros((num_pts, n_atom))
                data[ndx] = f_mag
            resids = list(meta['resSeq'])
            names = list(meta['resName'])
            labels = [str(resids[rndx]) + '_' + names[rndx] for rndx in range(len(resids))]
            df = pd.DataFrame(data, columns=labels)
            df['frame'] = df.index
            dfm = pd.melt(df, id_vars=['frame'], value_vars=labels)
            dfm['unit'] = ['frame' for i in range(len(dfm.index))]
            ax = sns.tsplot(data=dfm, time='frame', value='value', 
                    condition='variable', unit='unit')
            plt.savefig(f_pref + force + '_' + a_type + '.png')

if __name__ == '__main__':
    main()
