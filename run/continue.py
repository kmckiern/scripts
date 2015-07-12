from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import os
import time
import numpy as np
import argparse
from equil import setup_sim, dynamix

parser = argparse.ArgumentParser(description='equilibrate structures')
parser.add_argument('--sys', type=str, help='system pdb preface')
parser.add_argument('--pdb', type=str, help='IC pdb')
parser.add_argument('--nmin', type=int, help='number of minimization steps', default=50)
parser.add_argument('--nstep', type=int, help='number of steps')
args = parser.parse_args()

systm = args.sys
ns = args.nstep

# load initial parameters and geometry
prmtop = app.AmberPrmtopFile(systm + '.prmtop')
pdb = app.PDBFile(args.pdb)

# eq temp
temp = 300.0
# timestep
ts = 2.0

qs = pdb.positions
top = pdb.topology
unit_cell = top.getUnitCellDimensions()
box = unit_cell*np.eye(3)

# run it!
sim = setup_sim(prmtop, temp, ts, qs, 'gpu', top, box)
dynamix(systm, sim, ns, prmtop, temp, ts, 'gpu', min=args.nmin)
