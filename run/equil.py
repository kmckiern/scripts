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
from simtk.unit import Quantity, nanometer, angstrom

parser = argparse.ArgumentParser(description='equilibrate structures')
parser.add_argument('--pdb', type=str, help='system pdb preface')
parser.add_argument('--nstep', type=int, help='number of steps for final eq step', default=5000000)
parser.add_argument('--sint', type=int, help='save interval', default=500000)
parser.add_argument('--save', type=str, help='file specification to which we are saving', default='10ns')
args = parser.parse_args()

systm = args.pdb
ns = args.nstep
sint = args.sint
save = args.save
print(systm)
# load initial parameters and geometry
prmtop = app.AmberPrmtopFile(systm + '.prmtop')
inpcrd = app.AmberInpcrdFile(systm + '.inpcrd')
pdb = app.PDBFile(systm + '.pdb')

# benchmarking
current_milli_time = lambda: int(round(time.time() * 1000))
def milli_to_hr(millis):
    return (1.0 * millis) / (3600000.0)
def t_diff(start, finish, ts, steps):
    elapsed = milli_to_hr(finish-start)
    ns = (ts * steps * 1.0) / (1000000.0)
    return elapsed, elapsed / ns
# simulation setup and running
def setup_sim(prmtop, temperature, timestep, coordinates, which_pu, top, box):
    # update topology box vectors so they're save to the pdb correctly
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
    return simulation, system, integrator
def dynamix(systm, simulation, ns, prmtop, temperature, timestep, which_pu, min=False, print_box=False):
    if simulation == None:
        positions = pdb.positions # inpcrd.positions
        box = inpcrd.getBoxVectors()
    else:
        positions = simulation.context.getState(getPositions=True).getPositions()
        box = simulation.context.getState().getPeriodicBoxVectors()
    if min != False:
        top = pdb.topology # prmtop.topology
        bv = np.diag(box)
        qvs = [i.value_in_unit(nanometer) for i in bv]
        top.setUnitCellDimensions(Quantity(value=tuple(qvs), unit=nanometer))
        simulation, system, integrator = setup_sim(prmtop, temperature, timestep, positions, 'cpu', top, box)
        simulation.minimizeEnergy(maxIterations=min)
        positions = simulation.context.getState(getPositions=True).getPositions()
        of = systm + '_min.pdb'
        app.PDBFile.writeFile(top, positions, open(of, 'w'))
    slen = str(ns*timestep/1000.0).replace('.', 'p')
    of = systm + '_eq_' + str(temperature).split('.')[0] + '_' + slen + '.pdb'
    if print_box:
        with open(of + '.box', 'w') as file:
            for item in box:
                file.write('{}\n'.format(item))
    top = pdb.topology # prmtop.topology
    bv = np.diag(box)
    qvs = [i.value_in_unit(nanometer) for i in bv]
    top.setUnitCellDimensions(Quantity(value=tuple(qvs), unit=nanometer))
    simulation, system, integrator = setup_sim(prmtop, temperature, timestep, positions, which_pu, top, box)
    simulation.context.setVelocitiesToTemperature(temperature*unit.kelvin)
    # save trj every sint
    simulation.reporters.append(app.DCDReporter('trajectory.dcd', sint))
    simulation.reporters.append(app.StateDataReporter(stdout, sint, step=True,
        potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
        speed=True, totalSteps=ns, separator='\t'))
    simulation.step(ns)
    # record results
    positions = simulation.context.getState(getPositions=True).getPositions()
    box = simulation.context.getState().getPeriodicBoxVectors()
    top = prmtop.topology
    app.PDBFile.writeFile(top, positions, open(of, 'w'))
    system.setDefaultPeriodicBoxVectors(*box)
    return simulation, system, integrator
# serializing
def serializeObject(obj, objname):
    filename = objname
    objfile = open(filename,'w')
    objfile.write(mm.XmlSerializer.serialize(obj))
    objfile.close()

def main(): 
    # heating temperatures, 50 K increments
    t_min = 1.0
    t0_0, t1_0, t1_1, t1_2, t1_3, t2_0 = np.arange(50.0, 301.0, 50.0)
    # for debugging: t0_0, t1_0, t1_1, t1_2, t1_3, t2_0 = 10.0*np.ones(6)
    # timesteps
    ts0 = 1.0
    ts1 = 2.0
    # number of steps for each process
    minstep = 200
    nstep0 = 200
    nstep1 = 200
    nstep2 = 200 #ns
    
    #### simulations ####
    ## cpu / 1 fs timestep ##
    # minimize and 0 to 1 K for 1 ps
    smin, sysmin, imin = dynamix(systm, None, nstep0, prmtop, t_min, ts0, 'cpu', min=minstep)
    s4, sys4, i4 = dynamix(systm, smin, 2*nstep1, prmtop, t0_0, ts0, 'gpu')
    """
    ## gpu / 1 fs timestep ##
    # heating, 1 to 50 K for 100 ps
    s0, sys0, i0 = dynamix(systm, smin, 2*nstep1, prmtop, t0_0, ts0, 'gpu')
    ## gpu / 2 fs timestep ##
    # heating, 50 to 250 K for 100 ps each
    # 50-100
    s1, sys1, i1 = dynamix(systm, s0, nstep1, prmtop, t1_0, ts1, 'gpu')
    # 100-150
    s2, sys2, i2 = dynamix(systm, s1, nstep1, prmtop, t1_1, ts1, 'gpu')
    # 150-200
    s3, sys3, i3 = dynamix(systm, s2, nstep1, prmtop, t1_2, ts1, 'gpu')
    # 200-250
    s4, sys4, i4 = dynamix(systm, s3, nstep1, prmtop, t1_3, ts1, 'gpu')
    # heating, 250 to 300 K for 10 ns
    # 250-300
    """
    pre_eq = current_milli_time()
    simulation, system, integrator = dynamix(systm, s4, nstep2, prmtop, t0_0, ts1, 'gpu')
    post_eq = current_milli_time()
    total_eq, eq_rate = t_diff(pre_eq, post_eq, ts1, nstep2)
    print('bench: ' + str(eq_rate) + ' hr/ns')

    serializeObject(system, 'system_' + save + '.xml')
    serializeObject(integrator,'integrator_' + save + '.xml')
    state = simulation.context.getState(getPositions=True, getVelocities=True,
        getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
    serializeObject(state, 'state_' + save + '.xml')

if __name__ == '__main__':
    main()
