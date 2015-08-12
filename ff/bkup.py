#!/usr/bin/env python

import os, sys
import numpy as np
from collections import OrderedDict
import xml.etree.ElementTree as ET
import argparse
import IPython

parser = argparse.ArgumentParser(description='convert between amber99sb atom names and gmx/charmm equivalent')
parser.add_argument('--rtf', type=str, help='charmm rtf to parse')
parser.add_argument('--rxml', type=str, help='reference amber xml file')
parser.add_argument('--ox', type=str, help='name of output xml file')
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
    return AtAm, AtAc

def parse_amb_ref(off_path):
    AtAc = OrderedDict()
    et_xml = ET.parse(xml)
    root = et_xml.getroot()
    for force in root:
        if force.tag == 'Residues':
            for elem in force:
                res = elem.attrib['name']
                for subelem in elem:
                    if subelem.tag == 'Atom':
                        IPython.embed()
                        if res not in AtAc:
                            AtAc[res] = OrderedDict()
                        AtAc[res][subelem.attrib['name']] = subelem.attrib['type']
    return AtAc

def RTFtoXML(rtf_lines):
    print 'replacing'

def main():
    r_AtAm, r_AnAt = parse_rtf(args.rtf)
    IPython.embed()
    x_AnAt = parse_ref_xml(args.rxml)
    r_xml = RTFtoXML(rtf_l)

    with open(args.ox, 'wb') as f:
        r_xml.write(f)

if __name__ == '__main__':
    main()
