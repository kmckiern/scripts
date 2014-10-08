#!/usr/bin/env python
"""
this script converts an all-atom lipid monomer to a charmm-style united atom representation.
"""
import argparse
from forcebalance.molecule import Molecule
import mdtraj
import numpy as np

parser = argparse.ArgumentParser(description='parse ff params from atom list file to out file')
parser.add_argument('--xyz_file', type=str, help='path to coordinate file')
parser.add_argument('--out_xyz', type=str, help='path to list of atoms to keep')
opts = parser.parse_args()

tf = 'temp.pdb'
c_strip = ['C18', 'C19', 'C20', 'C21', 'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30', 'C31', 'C37', 'C38', 'C39', 'C40', 'C41', 'C42', 'C43', 'C44', 'C45', 'C46', 'C47', 'C48', 'C49', 'C50']
op = 'out.pdb'

def main():
    # create molecule
    m = Molecule(opts.xyz_file)
    # write to pdb
    m.write(tf)
    # read into mdtraj
    p = mdtraj.load(tf)
    # get list of hydrogens to remove
    c_ndx = [i.index for i in p.topology.atoms if i.name in c_strip]
    tail_atms = [y for (x, y) in m.Data['bonds'] if x in c_ndx]
    tail_h = [i for i in tail_atms if i not in c_ndx]
    # throw out some hydrogens
    atoms_to_keep = [b.index for b in p.topology.atoms if b.index not in tail_h]
    p.restrict_atoms(atoms_to_keep)
    p.save(op)
    # hella circular 
    x = Molecule(op)
    x.write(opts.out_xyz)

if __name__ == '__main__':
    main()
