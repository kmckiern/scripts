#!/usr/bin/env python

import os, sys
import numpy as np
from collections import OrderedDict
import xml.etree.ElementTree as ET
import parmed
import argparse
from operator import itemgetter
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
        elif specifier.startswith('RESI') or specifier.startswith('PRES'):
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
    for space in range(abs(diff)):
        if diff < 0:
                old += ' '
        else:
                new += ' '
    return old, new

def rewrite_rtf(lines, rmap):
    res = None
    for i, line in enumerate(lines):
        l = line.split('!')[0].strip()
        s = l.split()
        if len(s) == 0: continue
        specifier = s[0]
        # determine which residue we're parsing
        if specifier.startswith('RESI') or specifier.startswith('PRES'):
            spec, res, AACharge = s
            if res in rmap.keys():
                res_rm = rmap[res]
            else:
                res_rm = None
        # replace old atom types
        if specifier in ['ATOM', 'BOND', 'IMPROPER', 'DONOR', 'ACCEPTOR', 'IC', 'DELETE']:
            if res_rm != None:
                # if patching, need atom types from non terminal version of the residue
                if specifier == 'DELETE':
                    res_rm.update(rmap[res[1:]])
                for At in s:
                    if At in res_rm:
                        old, new = format_at(At, res_rm[At])
                        # bond v justification is tricky, so here's a hack
                        if len(new) > 3 and specifier == 'BOND':
                            new += ' '
                            if s.index(At) == len(s) - 1:
                                old = old[:-1]
                        # weirdly picky about these
                        if 'OT1' in old or 'OT2' in old:
                            continue
                        else:
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

    same_res = set(c_AtAc).intersection(set(a_AtAc))
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
                overlap = []
                for a_ndx, au in enumerate(a_uniq):
                    if cr[cu] == ar[au]:
                        overlap.append(len(set(cu).intersection(set(au))))
                    else:
                        overlap.append(0)
                if len(overlap) > 0:
                    # replace most similar corresponding AT
                    a_ndx = max(enumerate(overlap), key=itemgetter(1))[0]
                    au = a_uniq[a_ndx]
                    match[cu] = au
                    a_uniq.pop(a_ndx)
            rmap[res] = match
        else:
            print 'something is wrong: ', res
 
    out = rewrite_rtf(lines, rmap)
    with open('test.rtf', 'wb') as f:
        for line in out:
            f.write(line)

    print 'residues not parsed: ', set(c_AtAc) - set(a_AtAc)

if __name__ == '__main__':
    main()
