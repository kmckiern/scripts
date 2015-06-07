from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import os
import time
import numpy as np
import argparse
import mdtraj

parser = argparse.ArgumentParser(description='equilibrate structures')
parser.add_argument('--pdb', type=str, help='system pdb preface')
args = parser.parse_args()
systm = args.pdb
print(systm)

# benchmarking
current_milli_time = lambda: int(round(time.time() * 1000))
def milli_to_hr(millis):
    return (1.0 * millis) / (3600000.0)
def t_diff(start, finish, ts, steps):
    elapsed = milli_to_hr(finish-start)
    ns = (ts * steps * 1.0) / (1000000.0)
    return elapsed, elapsed / ns

# simulation setup and running
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
def dynamix(simulation, ns, temperature, timestep, which_pu, min=False):
    if min:
        positions = inpcrd.positions
        box = inpcrd.getBoxVectors()
        top = prmtop.topology
        simulation = setup_sim(temperature, timestep, positions, which_pu, box, top)
        simulation.minimizeEnergy(maxIterations=minstep)
    positions = simulation.context.getState(getPositions=True).getPositions()
    box = simulation.system.getDefaultPeriodicBoxVectors()
    top = simulation.topology
    simulation = setup_sim(temperature, timestep, positions, which_pu, box, top)
    simulation.context.setVelocitiesToTemperature(temperature*unit.kelvin)
    simulation.step(ns)
    return simulation

# heating temperatures, 50 K increments
t0_0, t1_0, t1_1, t1_2, t1_3, t2_0 = np.arange(50.0, 301.0, 50.0)
# timesteps
ts0 = 1.0
ts1 = 2.0
# number of steps for each process
minstep = 5000
nstep0 = 10000
nstep1 = 50000
nstep3 = 5000000

# load initial parameters and geometry
prmtop = app.AmberPrmtopFile(systm + '.prmtop')
inpcrd = app.AmberInpcrdFile(systm + '.inpcrd')

#### simulations ####
## cpu / 1 fs timestep ##
# minimize and 0 to 50 K for 50 ps
s0 = dynamix(None, nstep0, t1_0, ts0, 'cpu', min=True)
## gpu / 2 fs timestep ##
# heating, 50 to 250 K for 100 ps each
# 50-100
s1 = dynamix(s0, nstep1, t1_0, ts1, 'gpu')
# 100-150
s2 = dynamix(s1, nstep1, t1_1, ts1, 'gpu')
# 150-200
s3 = dynamix(s2, nstep1, t1_2, ts1, 'gpu')
# 200-250
s4 = dynamix(s3, nstep1, t1_3, ts1, 'gpu')
# heating, 250 to 300 K for 10 ns
# 250-300
pre_eq = current_milli_time()
s5 = dynamix(s4, nstep2, t2_0, ts1, 'gpu')
post_eq = current_milli_time()
total_eq, eq_rate = t_diff(pre_eq, post_eq, ts1, nstep3)
print('bench: ' + eq_rate + ' hr/ns')

# record results
positions = s6.context.getState(getPositions=True).getPositions()
of = systm + '_eq_' + str(nstep3 * ts1 / 1000.0) + 'ns.pdb'
app.PDBFile.writeFile(s6.topology, positions, open(of, 'w'))
# for bookkeeping
x = mdtraj.load(of)
print('na: ' + x.n_atoms)
