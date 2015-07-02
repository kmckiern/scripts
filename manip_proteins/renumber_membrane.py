#!/bin/env python

import argparse
import IPython

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
    iterd = False

    # read in og pdb
    with open(pdb, 'r') as infile:
        og = infile.readlines()
    # parse party
    for line in og:
        l = line.split()
        # if line is a coordinate specification
        if len(l) == 11:
            # get residue name and number
            residue = l[3]
            resid = int(l[4])
            # if residue is a membrane residue
            if residue in mem_res:
                # renumbering mapped newid
                if residue == mem_res[0] and iterd == False:
                    if newid == 0:
                        newid = resid
                    else:
                        newid += 1
                        iterd = True
                if residue == mem_res[2]:
                    iterd = False
                print ('before: ', line, residue, resid, newid) 
                line = line.replace(str(resid), str(newid))
                print ('after: ', line) 
        of.write(line)

if __name__ == '__main__':
    main()
