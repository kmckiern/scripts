#!/bin/env python

import argparse

parser = argparse.ArgumentParser(description='renumber membrane residues')
parser.add_argument('--pdb', type=str, help='pdb file')
args = parser.parse_args()

"""
POPC: pa, pc, ol
POPE: pa, pe, ol
"""
mem_res = ['PA', 'PC', 'OL', 'PE']

def main():
    pdb = args.pdb
    pdb_o = 'edit_' + pdb
    of = open(pdb_o, 'w+')
    newid = 0

    # read in og pdb
    with open(pdb, 'r') as infile:
        og = infile.readlines()
    # parse party
    for line in og:
        l = line.split()
        if len(l) == 11:
            residue = l[3]
            resid = int(l[4])
            if residue in mem_res:
                # if new residue, note index
                if residue == mem_res[0]:
                    if newid == 0:
                        newid = resid
                    else:
                        newid += 1
                line.replace(str(resid), str(newid))
        of.write(line)

if __name__ == '__main__':
    main()
