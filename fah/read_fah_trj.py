#!/bin/env python

import argparse
import mdtraj
import numpy as np

parser = argparse.ArgumentParser(description='get rmsd if fah WU to reference')
parser.add_argument('--pn', type=str, help='project number', default=None)
parser.add_argument('--top', type=str, help='topology file')
parser.add_argument('--ref', type=str, help='reference file')
parser.add_argument('--out', type=str, help='name of output file', default='rmsd.dat')
args = parser.parse_args()

pmap = {}
pmap['9608'] = '2JOF_a99sc-v2-r1_tip3pfb'
pmap['9609'] = '2JOF_a99sc-v2-r1_tip4pfb'
pmap['9610'] = '1FME_a99sc-v2-r1_tip3pfb'
pmap['9611'] = '1FME_a99sc-v2-r1_tip4pfb'
pmap['9612'] = '2F4K_a99sc-v2-r1_tip3pfb'
pmap['9613'] = '2F4K_a99sc-v2-r1_tip4pfb'
pmap['9614'] = '2F21_a99sc-v2-r1_tip3pfb'
pmap['9615'] = '2F21_a99sc-v2-r1_tip4pfb'
pmap['9616'] = '2P6J_a99sc-v2-r1_tip3pfb'
pmap['9617'] = '2P6J_a99sc-v2-r1_tip4pfb'
pmap['9621'] = 'C025_a99sc-v2-r1_tip3pfb'
pmap['9622'] = 'C025_a99sc-v2-r1_tip4pfb'

def main():
    prj = args.pn
    if prj != None:
       t = '../tops/' + pmap[prj] + '.pdb' 
    else:    
        t = args.top
    # read in fah trj
    trj = mdtraj.load('positions.xtc', top=t)

    # strip water and ions from simulation
    pi = trj.top.select('protein')
    trj = trj.atom_slice(pi)

    # calculate rmsd of trj to xtal structure
    ref = mdtraj.load(args.ref)
    trj.superpose(ref)
    dists = mdtraj.rmsd(trj, ref)
    
    # print final rmsd value
    print dists[-1]
    # np.save(args.out, dists)

if __name__ == '__main__':
    main()
