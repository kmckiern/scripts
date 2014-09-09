#!/bin/bash/env python

import Avogadro as a
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--file', help='UA file to add hydrogens to')
(optz, args) = parser.parse_args()

def ua_to_aa(options):
    options = options.__dict__
    df = options['file']
    pref = df.split('.')[0]
    mol = a.MoleculeFile.readMolecule(df)
    mol.addHydrogens()
    return mol, pref

def main():
    mol, out = ua_to_aa(optz)
    of = out + '_h.pdb'
    a.MoleculeFile.writeMolecule(mol, of)

if __name__ == '__main__':
    main()
