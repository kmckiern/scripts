#!/usr/bin/env python

import os, sys
from collections import OrderedDict
import numpy as np
import parmed as pmd
from parmed.amber import parameters
import xml.etree.ElementTree as ET
import argparse
import IPython

parser = argparse.ArgumentParser(description='convert xml file to ipython xml object')
parser.add_argument('--xml', type=str, help='path to xml file')
args = parser.parse_args()

def xml_parameters(xml_parsed):
    root = xml_parsed.getroot()

    # OpenMM Atom Types to Atom Class
    OAtAc = OrderedDict()

    # OpenMM Atom Classes to Masses
    OAcMass = OrderedDict()

    # OpenMM Residue-Atom Names to Atom Class
    AA_OAc = OrderedDict()

    # OpenMM Atom Class to Parameter Mapping
    # (vdW sigma and epsilon in AKMA)
    ONbPrm = OrderedDict()
    OBondPrm = OrderedDict()
    OAnglePrm = OrderedDict()
    ODihPrm = OrderedDict()
    OImpPrm = OrderedDict()

    Params = parameters.ParameterSet()

    # Stage 1 processing: Read in force field parameters
    for force in root:
        # Top-level tags in force field XML file are:
        # Forces
        # Atom types
        # Residues
        if force.tag == 'AtomTypes':
            for elem in force:
                OAtAc[elem.attrib['name']] = elem.attrib['class']
                mass = float(elem.attrib['mass'])
                if elem.attrib['class'] in OAcMass and mass != OAcMass[elem.attrib['class']]:
                    print "Atom class mass not consistent"
                    raise RuntimeError
                OAcMass[elem.attrib['class']] = mass
            # printcool_dictionary(OAtAc)
        # Harmonic bond parameters
        if force.tag == 'HarmonicBondForce':
            for elem in force:
                att = elem.attrib
                BC = (att['class1'], att['class2'])
                BCr = (att['class2'], att['class1'])
                acij = tuple(sorted(BC))
                if acij in OBondPrm:
                    print acij, "already defined in OBndPrm"
                    raise RuntimeError
                b = float(att['length'])*10
                k = float(att['k'])/10/10/2/4.184
                # sub C* and N* for CS and NS
                # acij = sub_key(acij, CType99)
                OBondPrm[acij] = (k, b)
                # Pass information to ParmEd
                Params._add_bond(acij[0], acij[1], rk=k, req=b)
                # New Params object can't write frcmod files.
                # Params.bond_types[acij] = pmd.BondType(k, b)
        # Harmonic angle parameters.  Same as for harmonic bonds.
        if force.tag == 'HarmonicAngleForce':
            for elem in force:
                att = elem.attrib
                AC = (att['class1'], att['class2'], att['class3'])
                ACr = (att['class3'], att['class2'], att['class1'])
                if AC[2] >= AC[0]:
                    acijk = tuple(AC)
                else:
                    acijk = tuple(ACr)
                if acijk in OAnglePrm:
                    print acijk, "already defined in OAnglePrm"
                    raise RuntimeError
                t = float(att['angle'])*180/np.pi
                k = float(att['k'])/2/4.184
                # acijk = sub_key(acijk, CType99)
                OAnglePrm[acijk] = (k, t)
                # Pass information to ParmEd
                Params._add_angle(acijk[0], acijk[1], acijk[2], thetk=k, theteq=t)
                # New Params object can't write frcmod files.
                # Params.angle_types[acijk] = pmd.AngleType(k, t)
        # Periodic torsion parameters.
        if force.tag == 'PeriodicTorsionForce':
            for elem in force:
                att = elem.attrib
                def fillx(strin):
                    if strin == "" : return "X"
                    else: return strin
                c1 = fillx(att['class1'])
                c2 = fillx(att['class2'])
                c3 = fillx(att['class3'])
                c4 = fillx(att['class4'])
                DC = (c1, c2, c3, c4)
                DCr = (c4, c3, c2, c1)
                if c1 > c4:
                    # Reverse ordering if class4 is alphabetically before class1
                    acijkl = DCr
                elif c1 < c4:
                    # Forward ordering if class1 is alphabetically before class4
                    acijkl = DC
                else:
                    # If class1 and class4 are the same, order by class2/class3
                    if c2 > c3:
                        acijkl = DCr
                    elif c3 > c2:
                        acijkl = DC
                    else:
                        acijkl = DC
                keylist = sorted([i for i in att.keys() if 'class' not in i])
                dprms = OrderedDict()
                for p in range(1, 7):
                    pkey = "periodicity%i" % p
                    fkey = "phase%i" % p
                    kkey = "k%i" % p
                    if pkey in keylist:
                        dprms[int(att[pkey])] = (float(att[fkey])*180.0/np.pi, float(att[kkey])/4.184)
                dprms = OrderedDict([(p, dprms[p]) for p in sorted(dprms.keys())])
                # ParmEd dihedral list
                dihedral_list = []
                for p in dprms.keys():
                    f, k = dprms[p]
                    dihedral_list.append(pmd.DihedralType(k, p, f, 1.2, 2.0, list=dihedral_list))
                    # Pass information to ParmEd
                    dtyp = 'normal' if (elem.tag == 'Proper') else 'improper'
                    Params._add_dihedral(acijkl[0], acijkl[1], acijkl[2], acijkl[3],
                                         pk=k, phase=f, periodicity=p, dihtype=dtyp)
                if elem.tag == 'Proper':
                    if acijkl in ODihPrm:
                        print acijkl, "already defined in ODihPrm"
                        raise RuntimeError
                    # sub C* and N* for CS and NS
                    # acijkl_0 = acijkl
                    # acijkl = sub_key(acijkl, CType99)
                    ODihPrm[acijkl] = dprms
                    # acijkl = acijkl_0
                elif elem.tag == 'Improper':
                    if acijkl in OImpPrm:
                        print acijkl, "already defined in OImpPrm"
                        raise RuntimeError
                    # acijkl_0 = acijkl
                    # acijkl = sub_key(acijkl, CType99)
                    OImpPrm[acijkl] = dprms
                    # acijkl = acijkl_0
                    if len(dihedral_list) > 1:
                        print acijkl, "more than one interaction"
                        raise RuntimeError
                else:
                    raise RuntimeError
        # Nonbonded parameters
        if force.tag == 'NonbondedForce':
            for elem in force:
                sigma = float(elem.attrib['sigma'])*10
                epsilon = float(elem.attrib['epsilon'])/4.184
                atype = elem.attrib['type']
                aclass = OAtAc[atype]
                amass = OAcMass[aclass]
                msigeps = (amass, sigma, epsilon)
                if aclass in ONbPrm:
                    if ONbPrm[aclass] != msigeps:
                        print 'mass, sigma and epsilon for atom class %s not uniquely defined:' % aclass
                        print msigeps, ONbPrm[aclass]
                        raise RuntimeError
                else:
                    ONbPrm[aclass] = (amass, sigma, epsilon)

        # Residue definitions
        if force.tag == 'Residues':
            resnode = force

    return OBondPrm, OAnglePrm, ODihPrm, OImpPrm, ONbPrm

def main():
    xml = ET.parse(args.xml)
    xml_BondPrm, xml_AnglePrm, xml_DihPrm, xml_ImpPrm, xml_NbPrm = xml_parameters(xml)
    
    IPython.embed()

if __name__ == '__main__':
    main()
