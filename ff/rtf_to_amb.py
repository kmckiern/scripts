#!/usr/bin/env python

import os, sys
import numpy as np
from collections import OrderedDict
import xml.etree.ElementTree as ET
import parmed
import argparse
import IPython

parser = argparse.ArgumentParser(description='rewrite charmm rtf to have amber atom names')
parser.add_argument('--rtf', type=str, help='starting charmm rtf')
parser.add_argument('--op', type=str, help='reference amber param dir')
parser.add_argument('--ox', type=str, help='name of output rtf file')
args = parser.parse_args()

def parse_rtf(rtf):
    AtAm = OrderedDict()
    AtAc = OrderedDict()
    # read in original rft
    lines = open(rtf).readlines()
    for line in lines:
        l = line.split('!')[0].strip()
        s = l.split()
        if len(s) == 0: continue
        specifier = s[0]
        # build map between atom types and masses
        if specifier == 'MASS':
            spec, Anum, At, mass = s
            AtAm[At] = mass
        # profile atom name to atom type on a per residue basis
        elif specifier.startswith("RESI"):
            spec, AA, AACharge = s
            AtAc[AA] = dict()
        # populate residue indexed map from atom names to atom types
        elif specifier == 'ATOM':
            spec, An, At, ACharge = s
            AtAc[AA][An] = At
    return AtAc, lines

# i think these can only differ by one
def format_at(old, new):
    lo = len(old)
    ln = len(new)
    if lo == ln:
        return old, new
    diff = lo-ln
    if diff == 1:
        return old, new + ' '
    elif diff == -1:
        return old + ' ', new
    else:
        print 'formatting assumption incorrect: ', old, new
        return None

def rewrite_rtf(lines, rmap):
    res = None
    for i, line in enumerate(lines):
        l = line.split('!')[0].strip()
        s = l.split()
        if len(s) == 0: continue
        specifier = s[0]
        # determine which residue we're parsing
        if specifier.startswith("RESI"):
            spec, res, AACharge = s
            if res in rmap.keys():
                res_rm = rmap[res]
            else:
                res_rm = None
        # replace old atom types
        elif specifier == 'ATOM':
            spec, At, Ac, ACharge = s
            if res_rm != None and At in res_rm:
                old, new = format_at(At, res_rm[At])
                line = line.replace(old, new)
            lines[i] = line
    return lines

def parse_amb_ref(off_path):
    # Load AMBER OFF libraries for amino acids
    amb_off = ['all_aminofb15.lib', 'all_aminoctfb15.lib', 'all_aminontfb15.lib']
    amb_resatoms = OrderedDict()
    for fin in amb_off:
        ambroot = off_path
        AOff = parmed.modeller.offlib.AmberOFFLibrary.parse(os.path.join(ambroot, 'dat', 'leap', 'lib', fin))
        for k, v in AOff.items():
            for a in v.atoms:
                if k not in amb_resatoms:
                    amb_resatoms[k] = {}
                amb_resatoms[k][a.name] = a.type
    return amb_resatoms

def main():
    c_AtAc, lines = parse_rtf(args.rtf)
    a_AtAc = parse_amb_ref(args.op)

    same_res = set(c_AtAc).intersection(a_AtAc)
    rmap = {}
    for res in same_res:
        cr = c_AtAc[res]
        ar = a_AtAc[res]
        c_uniq = list(set(cr)-set(ar))
        a_uniq = list(set(ar)-set(cr))
        if len(c_uniq) == len(a_uniq):
            # for each unique atom type, try to match by atom class
            match = {}
            for cu in c_uniq:
                for au in a_uniq:
                    if cr[cu] == ar[au]:
                        match[cu] = au
            rmap[res] = match
        else:
            print 'something is wrong: ', res
 
    out = rewrite_rtf(lines, rmap)
    with open('test.rtf', 'wb') as f:
        for line in out:
            f.write(line)

if __name__ == '__main__':
    main()
