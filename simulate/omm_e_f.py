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
ex usage
    >> python ~/scripts/ff/omm_energies_forces.py --pff whatever.xml --wff tip3p.xml --pdb FD/JIC/last.pdb --od omm
"""

parser = argparse.ArgumentParser(description='get energy and force components of a system for arbitrary FF')
parser.add_argument('--pff', type=str, help='protein FF')
parser.add_argument('--wff', type=str, help='water FF')
parser.add_argument('--pdb', type=str, help='pdb')
parser.add_argument('--od', type=str, help='output directory')
args = parser.parse_args()

# these are slightly modified version of some of Robert's code i saw on a github issue
def forcegroupify(system):
    forcegroups = {}
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        force.setForceGroup(i)
        forcegroups[force] = i
    return forcegroups

def getForceDecomposition(sim, forcegroups):
    forces = {}
    prettier = {}
    for f, i in forcegroups.items():
        forces[f] = sim.context.getState(getForces=True, groups=2**i).getForces()
        prettier[str(f)] = str(forces[f])
    return forces, prettier

def getEnergyDecomposition(sim, forcegroups):
    energies = {}
    prettier = {}
    for f, i in forcegroups.items():
        energies[f] = sim.context.getState(getEnergy=True, groups=2**i).getPotentialEnergy()
        name = str(f).split(';')[0].split('openmm.')[-1]
        prettier[name] = str(energies[f])
    return energies, prettier

ff = app.ForceField(args.pff, args.wff)
platform = mm.Platform.getPlatformByName('CPU')
pdb = app.PDBFile(args.pdb)
integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds,
    2.0*unit.femtoseconds)
integrator.setConstraintTolerance(0.00001)
system = ff.createSystem(pdb.topology, nonbondedMethod=app.PME,
    nonbondedCutoff=2.0*unit.nanometers, constraints=app.HBonds, rigidWater=True,
    ewaldErrorTolerance=0.0005, vdwCutoff=2.0*unit.nanometer)
system.addForce(mm.MonteCarloBarostat(1*unit.atmospheres, 300*unit.kelvin))
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

state = simulation.context.getState(getEnergy=True, getForces=True)
pe = state.getPotentialEnergy()
forces = state.getForces()

od = args.od

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
        f.write(i.replace(' kJ/(nm mol)', '').replace('[(', '').replace(')]', '').replace(' (', '') + '\n')
    f.close()
