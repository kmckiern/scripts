#!/bin/env python

"""
featurize clc trjs using the following order parameters:

note regarding indexing:
    suppose we'd like the atom indices of residue i (wrt pdb numbering)
    these are equivalent:
        - [atom.index for atom in trj.top.atoms if atom.residue.resSeq == i]
        - trj.top.select('resid i-1')
"""

import argparse
import mdtraj
import numpy as np
import sys
sys.path.insert(0, '/Users/kerimckiernan/Dropbox/scripts/py_general')
from toolz import cartesian
import IPython

# voltage sensitive residues of interest
vs_ri = {
    'lys127': 127,
    'glu128': 128,
    'asp32': 32,
    'lys317': 317,
    'glu318': 318,
}

# get CA atom indices for these residues
def CA(trj):
    ca_dict = {}
    for resid in vs:
        ca_dict[resid] = [atom.index for atom in trj.top.atoms if ((atom.residue.resSeq == resid) and (atom.name == 'CA'))][0]
    return ca_dict

# bc numbering is weird
def res_ndxs(trj, resid):
    resid -= 1
    return trj.top.select('resid ' + str(resid))

def ele_ndxs(trj, resid, ele_ls):
    return [atom.index for atom in trj.top.atoms if ((atom.residue.resSeq == resid) and (atom.name in ele_ls))]

def main():
    parser = argparse.ArgumentParser(description='custom featurization of clc fah trjs')
    parser.add_argument('--ref', type=str, help='homology model pdb file')
    parser.add_argument('--trj', type=str, help='trajectory file')
    parser.add_argument('--mol2', type=str, help='homology model mol2 file (charges needed for dipole calc)')
    args = parser.parse_args()
     
    # load system data
    trj = mdtraj.load(args.trj, top=args.ref)
    hmodel = mdtraj.load(args.ref)
    
    ### feature 0: protein RMSD from hmodel ###
    pi_noh = [atom.index for atom in trj.top.atoms if ((atom.residue.is_protein) and (atom.element.symbol != 'H'))]
    p_rmsd = mdtraj.rmsd(trj, hmodel, atom_indices=pi_noh)

    ### feature 1: GLU128 RMSD from hmodel ###
    e128 = res_ndxs(hmodel, vs_ri['glu128'])
    e128_rmsd = mdtraj.rmsd(trj, hmodel, atom_indices=e128)
    
    ### feature 2: LYS317 and GLU318 RMSD from hmodel ###
    tl = np.concatenate((res_ndxs(hmodel, vs_ri['lys317']), res_ndxs(hmodel, vs_ri['glu318'])))
    tl_rmsd = mdtraj.rmsd(trj, hmodel, atom_indices=tl)

    ### feature 2: distance between ASP32 and LYS127 ###
    a32 = ele_ndxs(hmodel, vs_ri['asp32'], ['OD1', 'OD2'])
    l127 = ele_ndxs(hmodel, vs_ri['lys127'], ['NZ'])
    al_pairs = cartesian([a32, l127])
    # i think the asp oxygens are degenerate, so i'll look at the min here
    al_dist = np.amin(al_pairs, axis=1)

    ### feature 4: dipole moment over voltage sensitive residues ***
    # atoms, bonds = mdtraj.formats.mol2.mol2_to_dataframes(args.mol2)

if __name__ == '__main__':
    main()
