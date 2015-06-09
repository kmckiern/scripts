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
import IPython

# benchmarking
current_milli_time = lambda: int(round(time.time() * 1000))
def milli_to_hr(millis):
    return (1.0 * millis) / (3600000.0)
def t_diff(start, finish, ts, steps):
    elapsed = milli_to_hr(finish-start)
    ns = (ts * steps * 1.0) / (1000000.0)
    return elapsed, elapsed / ns
# simulation setup and running
def setup_sim(prmtop, temperature, timestep, coordinates, which_pu, top, box=None):
    system = prmtop.createSystem(nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True,
        ewaldErrorTolerance=0.0005)
    # if a pdb is loaded, box vectors do not need to be loaded (included in topology)
    # this is not the case for an inpcrd, i think
    if box != None:
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
def dynamix(simulation, ns, prmtop, temperature, timestep, which_pu, min=False, print_box=False):
    if min:
        positions = inpcrd.positions
        box = inpcrd.getBoxVectors()
        top = prmtop.topology
        simulation = setup_sim(temperature, timestep, positions, which_pu, box, top)
        simulation.minimizeEnergy(maxIterations=minstep)
        positions = simulation.context.getState(getPositions=True).getPositions()
        of = systm + '_min.pdb'
        app.PDBFile.writeFile(top, positions, open(of, 'w'))
    positions = simulation.context.getState(getPositions=True).getPositions()
    box = simulation.context.getState().getPeriodicBoxVectors()
    of = systm + '_eq_' + str(temperature).split('.')[0] + '_' + str(ns*timestep/1000.0) + '.pdb'
    if print_box:
        with open(of + '.box', 'w') as file:
            for item in box:
                file.write('{}\n'.format(item))
    top = simulation.topology
    simulation = setup_sim(temperature, timestep, positions, which_pu, box, top)
    simulation.context.setVelocitiesToTemperature(temperature*unit.kelvin)
    simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, 
        potentialEnergy=True, temperature=True, volume=True, progress=True, 
        remainingTime=True, speed=True, totalSteps=1000, separator='\t'))
    simulation.step(ns)
    # record results
    positions = simulation.context.getState(getPositions=True).getPositions()
    box = simulation.context.getState().getPeriodicBoxVectors()
    top = simulation.topology
    app.PDBFile.writeFile(top, positions, open(of, 'w'))
    return simulation

def main():
    parser = argparse.ArgumentParser(description='equilibrate structures')
    parser.add_argument('--pdb', type=str, help='system pdb preface')
    args = parser.parse_args()
    systm = args.pdb
    print(systm)
    # load initial parameters and geometry
    prmtop = app.AmberPrmtopFile(systm + '.prmtop')
    inpcrd = app.AmberInpcrdFile(systm + '.inpcrd')
    
    # heating temperatures, 50 K increments
    t0_0, t1_0, t1_1, t1_2, t1_3, t2_0 = np.arange(50.0, 301.0, 50.0)
    # timesteps
    ts0 = 1.0
    ts1 = 2.0
    # number of steps for each process
    minstep = 1000
    nstep0 = 10000
    nstep1 = 50000
    nstep2 = 5000000
    
    #### simulations ####
    ## cpu / 1 fs timestep ##
    # minimize and 0 to 50 K for 10 ps
    s0 = dynamix(None, nstep0, prmtop, t0_0, ts0, 'cpu', min=True)
    ## gpu / 1 fs timestep ##
    # heating, 50 to 100 K for 20 ps
    s1 = dynamix(s0, 2*nstep1, prmtop, t1_0, ts0, 'gpu')
    ## gpu / 2 fs timestep ##
    # heating, 100 to 250 K for 100 ps each
    # 100-150
    s2 = dynamix(s1, nstep1, prmtop, t1_1, ts1, 'gpu')
    # 150-200
    s3 = dynamix(s2, nstep1, prmtop, t1_2, ts1, 'gpu')
    # 200-250
    s4 = dynamix(s3, nstep1, prmtop, t1_3, ts1, 'gpu')
    # heating, 250 to 300 K for 10 ns
    # 250-300
    pre_eq = current_milli_time()
    s5 = dynamix(s4, nstep2, prmtop, t2_0, ts1, 'gpu')
    post_eq = current_milli_time()
    total_eq, eq_rate = t_diff(pre_eq, post_eq, ts1, nstep2)
    print('bench: ' + str(eq_rate) + ' hr/ns')

if __name__ == '__main__':
    main()
