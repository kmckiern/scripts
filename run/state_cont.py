from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import os
import time
import numpy as np
import argparse
import IPython

# serializing
def serializeObject(obj, filename):
    objfile = open(filename,'w')
    objfile.write(mm.XmlSerializer.serialize(obj))
    objfile.close()

parser = argparse.ArgumentParser(description='continue system equil from state xmls')
parser.add_argument('--pdb', type=str, help='IC pdb')
parser.add_argument('--nmin', type=int, help='number of minimization steps', default=100)
parser.add_argument('--nstep', type=int, help='number of steps')
parser.add_argument('--sint', type=int, help='save interval', default=500000)
parser.add_argument('--prev', type=str, help='file specification from which we are resuming')
parser.add_argument('--save', type=str, help='file specification to which we are saving')
args = parser.parse_args()

ns = args.nstep
prev = args.prev
save_spec = args.save

system = mm.XmlSerializer.deserializeSystem(open('system_' + prev + '.xml').read())
state = mm.XmlSerializer.deserialize(open('state_' + prev + '.xml').read())
integrator = mm.XmlSerializer.deserialize(open('integrator_' + prev + '.xml').read())

platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}

pdb = app.PDBFile(args.pdb)
simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)
simulation.context.setPositions(state.getPositions())
simulation.context.setVelocities(state.getVelocities())
simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
for key, value in state.getParameters().iteritems():
    simulation.context.setParameter(key, value)

# save every 1 ns
save_interval = args.sint
simulation.reporters.append(app.DCDReporter('trajectory_' + save_spec + '.dcd', save_interval))
simulation.reporters.append(app.StateDataReporter(stdout, save_interval, step=True,
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
    speed=True, totalSteps=ns, separator='\t'))

# run
simulation.minimizeEnergy(maxIterations=args.nmin)
simulation.step(ns)

# save
serializeObject(system, 'system_' + save_spec + '.xml')
serializeObject(integrator,'integrator_' + save_spec + '.xml')
state = simulation.context.getState(getPositions=True, getVelocities=True,
    getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
serializeObject(state, 'state_' + save_spec + '.xml')
