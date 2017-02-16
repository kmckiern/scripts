from __future__ import print_function
import time
from simtk.openmm import app
from simtk.openmm.app import AmberPrmtopFile
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import os
import numpy as np
import mdtraj
import pandas as pd
import argparse
import IPython

parser = argparse.ArgumentParser(description='get time series of nonbonded force for an input trajectory')
parser.add_argument('--prmtop', type=str, help='prmtop file')
parser.add_argument('--trj', type=str, help='trajectory file', default='raw_01.xtc')
parser.add_argument('--ri', type=str, help='residues of interest')
parser.add_argument('--od', type=str, help='output directory', default='./forces/')
args = parser.parse_args()

od = args.od
useful = ['HarmonicAngleForce', 'HarmonicBondForce', 'PeriodicTorsionForce', 'NonbondedForce']

res_indxs = args.ri
for j in range(4):
    for i in np.arange(-3,4):
        res_indxs.append(res_indxs[j]+i)

def force_groups(system):
    forcegroups = {}
    for i in range(system.getNumForces()):
        force = system.getForce(i)

        force.setForceGroup(i)
        forcegroups[force] = i
    return forcegroups

def force_components(sim, forcegroups):
    forces = {}
    prettier = {}
    for f, i in forcegroups.items():
        forces[f] = sim.context.getState(getForces=True, groups=2**i).getForces()
        prettier[str(f)] = str(forces[f])
    return forces, prettier

def energy_components(sim, forcegroups):
    energies = {}
    prettier = {}
    for f, i in forcegroups.items():
        energies[f] = sim.context.getState(getEnergy=True, groups=2**i).getPotentialEnergy()
        name = str(f).split(';')[0].split('openmm.')[-1]
        prettier[name] = str(energies[f])
    return energies, prettier

def relevant_indices(trj):
    atom_indxs = []
    for ref in res_indxs:
        atom_indxs += [atom.index for atom in trj.top.atoms if atom.residue.resSeq == ref]
    ri = np.array(atom_indxs)
    wi = np.array([atom.index for atom in trj.topology.atoms if (atom.residue.is_water and (atom.element.symbol == 'O'))])
    ki = np.array([i.index for i in trj.top.atoms_by_name('K+')])
    all_i = np.concatenate((ri,wi,ki))

    table, bonds = trj.topology.to_dataframe()
    data = table.T[all_i].T
    return data

start_time = time.time()

prmtop = AmberPrmtopFile(args.prmtop)
platform = mm.Platform.getPlatformByName('CPU')
system = prmtop.createSystem(nonbondedMethod=app.PME,
    nonbondedCutoff=2.0*unit.nanometers, constraints=app.HBonds, rigidWater=True,
    ewaldErrorTolerance=0.0005)
system.addForce(mm.MonteCarloBarostat(1*unit.atmospheres, 300*unit.kelvin))

integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds,
    2.0*unit.femtoseconds)
integrator.setConstraintTolerance(0.00001)
simulation = app.Simulation(prmtop.topology, system, integrator, platform)

trj = mdtraj.load(args.trj, top=args.prmtop)
datumz = []
for tndx, frame in enumerate(trj):
    if tndx == 0:
        indices = relevant_indices(frame)
        num_indx = indices.index

    qs = frame.xyz[0]
    qst = qs[num_indx]
    indices['q'] = qst.tolist()

    simulation.context.setPositions(qs)
    state = simulation.context.getState(getEnergy=True, getForces=True)
    pe = state.getPotentialEnergy()

    fg = force_groups(system)
    # erngz, perngz = energy_components(simulation, fg)
    forces, pforce = force_components(simulation, fg)
    for fg in pforce.keys():
        name = fg.split(';')[0].split('openmm.')[-1]
        if name in useful:
            force_data = pforce[fg].replace(' kJ/(nm mol)','').replace('[(','').replace(')]','').replace(' (','').split('),')
            rows = [l.rstrip(',').split(',') for l in force_data]
            indices[name] = np.array(rows,float)[num_indx].tolist()
    end_time = time.time()
    print ('% .2f' % (end_time - start_time))
    start_time = time.time()
    # write data to pickle
    indices.to_pickle(od + args.trj.split('/')[-1].split('.')[0] + '_' + str(tndx) + '.pkl')

IPython.embed()
