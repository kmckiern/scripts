#!/usr/bin/env python

"""
reimage a molecule to the center of a pbc and cut beyond a certain radius
example usage:
    >> python /path/to/script/dir/image.py --p input.pdb --t input.dcd --m mol_num --cd cutdist -o out
"""

import argparse
import mdtraj
import copy

parser = argparse.ArgumentParser(description='replicate crystal pdb to arbitrary size')
parser.add_argument('--p', type=str, help='input pdb file')
parser.add_argument('--b', type=str, help='big input pdb file')
parser.add_argument('--t', type=str, help='input trj')
parser.add_argument('--m', type=int, help='molecule to center')
parser.add_argument('--cd', type=str, help='distance beyond which molecules are cut')
parser.add_argument('--o', type=str, help='output pdb and trj name')
args = parser.parse_args()

#def center_mol(trj, molnum):
    
def sandwich(trj, ndx):
    tm = copy.deepcopy(trj)
    tm.xyz -= trj.unitcell_vectors[0][:,ndx]
    tp = copy.deepcopy(trj)
    tp.xyz += trj.unitcell_vectors[0][:,ndx]
    return tm.stack(tp)

def image_trj(trj):
    # stack x for first row
    rank1 = trj.stack(sandwich(trj, 0))
    # stack y for matrix    
    rank2 = rank1.stack(sandwich(rank1, 1))
    # now, 3d
    rank3 = rank2.stack(sandwich(rank2, 2))
    return rank3

def indx_dists(trj, center_mol, cut):
    nmols = trj.xyz.shape[1]
    pairs = np.zeros((nmols, 2))
    pairs[:,0] = np.arange(nmols)
    pairs[:,1] = center_mol*np.ones(nmols)
    ds = mdtraj.compute_distances(trj, pairs)
    return [i for i in ds if i < cut]

def main():
    # get input
    p = mdtraj.load(args.t, top=args.p)

    # image
    p_im = image_trj(p)
    p_im[0].save_pdb('zero.pdb')
    p_im.save_dcd('big.dcd')

    # couldn't figure out how to merge chains w/o chimera
    p_big = mdtraj.load('big.dcd', top='ftlog.pdb')
    # fix first oxygen
    msel = p_big.xyz[:,0,:]
    # shift coords
    for t in range(p_big.xyz.shape[0]):
        p_big.xyz[t] -= msel[i]
    # write new dcd
    p_big.write_dcd('shifted.dcd')

if __name__ == '__main__':
    main()
