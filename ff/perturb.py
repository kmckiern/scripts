#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np

parser = argparse.ArgumentParser(description='perturb pdb coords')
parser.add_argument('--pdb', type=str, help='pdb file')
parser.add_argument('--pert', type=float, help='perturbation')
args = parser.parse_args()

def main():
    zero = md.load(args.pdb)
    pert = args.pert
    ext = 6

    na = zero.top.n_atoms
    # coords
    new_coords = np.zeros((na*ext, na, 3))
    # box vecs
    z_box = zero.unitcell_vectors[0]
    vecs = np.zeros((na*ext, 3, 3))
    # for each atom, move 1/1000 of an angstrom in each direction
    for atom in range(na):
        for i in range(ext):
            t = (atom*6) + i
            # initialize
            new_coords[t] = zero.xyz[0]
            # perturb
            # x
            if i == 0:
                new_coords[t][atom][0] -= pert 
            if i == 1:
                new_coords[t][atom][0] += pert 
            # y
            if i == 2:
                new_coords[t][atom][1] -= pert 
            if i == 3:
                new_coords[t][atom][1] += pert 
            # z
            if i == 4:
                new_coords[t][atom][2] -= pert 
            if i == 5:
                new_coords[t][atom][2] += pert 
            vecs[t] = z_box
    # coords
    zero.xyz = new_coords
    # box vecs
    zero.unitcell_vectors = vecs
    # time arr
    zero.time = np.arange(na*ext)
    zero[0].save_pdb('perturb0.pdb')
    zero.save_xtc('perturb.xtc')

if __name__ == '__main__':
    main()
