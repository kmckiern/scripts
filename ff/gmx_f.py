from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import os
import time
import numpy as np
import json
import argparse

"""
read in g_energy files for an original configuration, and a perturbation trajectory, 
    and calculates force components for each atom
ex usage
    >> python ~/scripts/ff/gmx_f.py --up zero_ge.xvg --p tip3p.xml --h .001
"""

parser = argparse.ArgumentParser(description='get force components')
parser.add_argument('--up', type=str, help='unperturbed g_energy')
parser.add_argument('--p', type=str, help='perturbed g_energy')
parser.add_argument('--h', type=float, help='perturbation amount (angstroms)')
args = parser.parse_args()

def parse_zero(ef):
    return None

def parse_pertbd(ef):
    return None

fg = forcegroupify(system)
erngz, perngz = getEnergyDecomposition(simulation, fg)
f = open(od + '/energies.dat', 'w')
f.write('total energy: ' + str(pe) + '\n')
f.write('component analysis: \n')
f.write(json.dumps(perngz, indent=2) + '\n')
f.close()

forces, pforce = getForceDecomposition(simulation, fg)
for fg in pforce.keys():
    name = fg.split(';')[0].split('openmm.')[-1]
    f = open(od + '/' + name + '.dat', 'w')
    for i in pforce[fg].split('),'):
        f.write(i.replace('[(', '').replace(')]', '').replace(' (', '') + '\n')
    f.close()
