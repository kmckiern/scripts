from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from simtk.unit import Quantity, nanometer, angstrom
from sys import stdout
import os
import sys
import numpy as np
import argparse

# serializing
def serializeObject(obj, filename):
    objfile = open(filename,'w')
    objfile.write(mm.XmlSerializer.serialize(obj))
    objfile.close()

parser = argparse.ArgumentParser(description='continue system equil from state xmls')
parser.add_argument('--prmtop', type=str, help='IC prmtop')
parser.add_argument('--nmin', type=int, help='number of minimization steps', default=100)
parser.add_argument('--nstepc', type=int, help='number of cpu steps', default=5000)
parser.add_argument('--nstep', type=int, help='number of steps')
parser.add_argument('--sint', type=int, help='save interval', default=500000)
parser.add_argument('--prev', type=str, help='file specification from which we are resuming')
parser.add_argument('--save', type=str, help='file specification to which we are saving')
args = parser.parse_args()
# globals read in from CL
prmtop = app.AmberPrmtopFile(args.prmtop)
ns = args.nstep
# 1 ns on cpu
tcpu = args.nstepc
# remainder on gpu
tgpu = ns - tcpu
prev = args.prev
save_spec = args.save

## CPU ##
# read in from xml files
system = mm.XmlSerializer.deserializeSystem(open('system_' + prev + '.xml').read())
state = mm.XmlSerializer.deserialize(open('state_' + prev + '.xml').read())
integrator = mm.XmlSerializer.deserialize(open('integrator_' + prev + '.xml').read())
# setup simulation
platform = mm.Platform.getPlatformByName('CPU')
simulation = app.Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(state.getPositions())
simulation.context.setVelocities(state.getVelocities())
simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
for key, value in state.getParameters().iteritems():
    simulation.context.setParameter(key, value)
# monitor every 1 ns
save_interval = args.sint
simulation.reporters.append(app.DCDReporter('trajectory_' + save_spec + '.dcd', save_interval))
simulation.reporters.append(app.StateDataReporter(stdout, save_interval, step=True,
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
    speed=True, totalSteps=ns, separator='\t'))
# look at initial periodic box dimensions
state = simulation.context.getState(getPositions=True, getVelocities=True,
    getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
box = simulation.context.getState().getPeriodicBoxVectors()
box_vecs = np.diag(box)
qx, qy, qz = box_vecs
qvs = [i.value_in_unit(angstrom) for i in box_vecs]
print('initial box')
print('box: ' + ''.join(str(qvs)))
# run after short minimization
simulation.minimizeEnergy(maxIterations=50)
simulation.step(tcpu)
# check box after simulation
state = simulation.context.getState(getPositions=True, getVelocities=True,
    getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
box = simulation.context.getState().getPeriodicBoxVectors()
box_vecs = np.diag(box)
qx, qy, qz = box_vecs
qvs = [i.value_in_unit(angstrom) for i in box_vecs]
print('cpu sim completed')
print('box: ' + ''.join(str(qvs)))

## GPU ##
# use cpu box vector state, reread system and integrator
system.setDefaultPeriodicBoxVectors(*box)
integrator = mm.XmlSerializer.deserialize(open('integrator_' + prev + '.xml').read())
# set up sim
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}
top = prmtop.topology
bv = np.diag(box)
qvs = [i.value_in_unit(nanometer) for i in bv]
top.setUnitCellDimensions(Quantity(value=tuple(qvs), unit=nanometer))
simulation = app.Simulation(top, system, integrator, platform, properties)
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
simulation.step(tgpu)
# look at box after steps
box = simulation.context.getState().getPeriodicBoxVectors()
box_vecs = np.diag(box)
qx, qy, qz = box_vecs
qvs = [i.value_in_unit(angstrom) for i in box_vecs]
print('gpu sim completed')
print('box: ' + ''.join(str(qvs)))
# save xmls
system.setDefaultPeriodicBoxVectors(*box)
serializeObject(system, 'system_' + save_spec + '.xml')
serializeObject(integrator,'integrator_' + save_spec + '.xml')
state = simulation.context.getState(getPositions=True, getVelocities=True,
    getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
serializeObject(state, 'state_' + save_spec + '.xml')
