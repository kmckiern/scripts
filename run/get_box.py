from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from simtk.unit import angstrom
from sys import stdout
import os
import time
import numpy as np
import argparse

"""
example usage:
while read LINE; do cd $LINE; python ~/scripts/run/get_box.py --suff 32ns; cd ..; done < dimer_dirs.dat
"""

parser = argparse.ArgumentParser(description='get box vectors from state xml')
parser.add_argument('--suff', type=str, help='state file suffix')
args = parser.parse_args()

suff = args.suff
state = mm.XmlSerializer.deserialize(open('state_' + suff + '.xml').read())

box_vecs = np.diag(state.getPeriodicBoxVectors())
qx, qy, qz = box_vecs
qvs = [i.value_in_unit(angstrom) for i in box_vecs]

of = 'box_' + suff + '.dat'
with open(of + '.box', 'w') as file:
    for item in qvs:
        file.write('{}\n'.format(item))
