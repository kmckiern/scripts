from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import os
import time
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='equilibrate structures')
parser.add_argument('--pdb', type=str, help='system pdb preface')
args = parser.parse_args()
systm = args.pdb

# benchmarking
current_milli_time = lambda: int(round(time.time() * 1000))
def milli_to_hr(millis):
    return (1.0 * millis) / (3600000.0)
def t_diff(start, finish, ts, steps):
    elapsed = milli_to_hr(finish-start)
    ns = (ts * steps * 1.0) / (1000000.0)
    return elapsed, elapsed / ns

# setting up simulation
def setup_sim(temperature, timestep, coordinates, which_pu, box, top):
    system = prmtop.createSystem(nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True,
        ewaldErrorTolerance=0.0005)
    system.setDefaultPeriodicBoxVectors(*box)
    integrator = mm.LangevinIntegrator(temperature*unit.kelvin, 1.0/unit.picoseconds,
        timestep*unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)
    system.addForce(mm.MonteCarloBarostat(1*unit.atmospheres, temperature*unit.kelvin))
    if which_pu == 'cpu':
        platform = mm.Platform.getPlatformByName('CPU')
        simulation = app.Simulation(top, system, integrator, platform)
    elif which_pu == 'gpu':
        platform = mm.Platform.getPlatformByName('CUDA')
        properties = {'CudaPrecision': 'mixed'}
        simulation = app.Simulation(top, system, integrator, platform, properties)
    simulation.context.setPositions(coordinates)
    return simulation

# heating temperatures
t0 = 100.0
t1 = 200.0
t2 = 300.0
# timesteps
ts0 = 1.0
ts1 = 2.0
# number of steps for each process
minstep = 1000
nstep0 = 1000
nstep1 = 2000
nstep3 = 2500000
# load starting parameters and geometry
prmtop = app.AmberPrmtopFile(systm + '.prmtop')
inpcrd = app.AmberInpcrdFile(systm + '.inpcrd')
box = inpcrd.getBoxVectors()

## simulating ##
# cpu / 1 fs timestep
simulation = setup_sim(t0, ts0, inpcrd.positions, 'cpu', box, prmtop.topology)
# minimize
simulation.minimizeEnergy(maxIterations=minstep)
positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open(systm + '_min.pdb', 'w'))
# heating, 0 to 100 K over 1 ps
simulation.context.setVelocitiesToTemperature(t0*unit.kelvin)
simulation.step(nstep0)

# gpu / 2 fs timestep
# heating, 100 to 200 K over 4 ps
box = simulation.system.getDefaultPeriodicBoxVectors()
top = simulation.topology
simulation = setup_sim(t1, ts1, positions, 'gpu', box, top)
simulation.context.setVelocitiesToTemperature(t1*unit.kelvin)
simulation.step(nstep1)
# heating, 200 to 300 K for 5 ns
box = simulation.system.getDefaultPeriodicBoxVectors()
top = simulation.topology
simulation = setup_sim(t2, ts1, positions, 'gpu', box, top)
simulation.context.setVelocitiesToTemperature(t2*unit.kelvin)
pre_eq = current_milli_time()
simulation.step(nstep3)
post_eq = current_milli_time()
total_eq, eq_rate = t_diff(pre_eq, post_eq, ts1, nstep3)
print('timing: ', total_eq, eq_rate)
positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open(systm + '_eq_5ns.pdb', 'w'))
