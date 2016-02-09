from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import os
import numpy as np
import argparse
import IPython
import mdtraj
import datetime
from forcebalance.molecule import Molecule

start_time = str(datetime.datetime.now().time())

parser = argparse.ArgumentParser(description='minimize water')
parser.add_argument('--prmtop', type=str, help='prmtop file')
parser.add_argument('--trj', type=str, help='trajectory file', default='raw_01.xtc')
args = parser.parse_args()

# custom force, harmonic positional constraint
def make_cef(k, hold=['HOH']):
    cef = mm.CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    cef.addPerParticleParameter('k')
    cef.addPerParticleParameter('x0')
    cef.addPerParticleParameter('y0')
    cef.addPerParticleParameter('z0')
    for i in range(M.na):
        if M.resname[i] not in hold:
            cef.addParticle(i, [k] + [M.xyzs[0][i][j]/10 for j in range(3)])
        else:
            cef.addParticle(i, [0.0] + [M.xyzs[0][i][j]/10 for j in range(3)])
    return cef

# run min or eq simulation, possibly with restraint via make_cef
def rest_sim(ts, qs, min_steps, k, dont_hold, eq_steps, platf):
    lf.write('pre: ' + str(datetime.datetime.now().time()))
    lf.write('\n')
    integ = mm.VerletIntegrator(ts*unit.picosecond)
    plat = mm.Platform.getPlatformByName(platf)
    system = prmtop.createSystem(nonbondedMethod=app.PME,
        nonbondedCutoff=1.4*unit.nanometers, rigidWater=True,
        constraints=app.HBonds, ewaldErrorTolerance=0.0005)
    # system.addForce(mm.MonteCarloBarostat(1*unit.atmospheres, 300*unit.kelvin))
    if k > 0:
        cef = make_cef(k, dont_hold)
        system.addForce(cef)
        lf.write('Created System with restraint' + str(k))
        lf.write('\n')
    sim = app.Simulation(pdb.topology, system, integ, plat)
    sim.context.setPositions(qs)
    init_V = sim.context.getState(getEnergy=True).getPotentialEnergy()
    lf.write('Created Simulation, Potential = ' + str(init_V))
    lf.write('\n')
    if min_steps > 0:
        lf.write('min steps: ' +  str(min_steps))
        lf.write('\n')
        sim.minimizeEnergy(maxIterations=min_steps)
        final_V = sim.context.getState(getEnergy=True).getPotentialEnergy()
        lf.write('Energy Minimized, Potential = ' + str(final_V))
        lf.write('\n')
    if eq_steps > 0:
        lf.write('eq steps: ' + str(min_steps))
        lf.write('\n')
        sim.context.setVelocitiesToTemperature(300*unit.kelvin)
        sim.step(eq_steps)
        final_V = sim.context.getState(getEnergy=True).getPotentialEnergy()
        lf.write('Energy Eq, Potential = ' + str(final_V))
        lf.write('\n')
    qs = sim.context.getState(getPositions=True).getPositions()
    lf.write('post: ' + str(datetime.datetime.now().time()))
    lf.write('\n')
    return qs, final_V

# set everything up
pref = args.trj.split('.')[0]
log = pref + '.log'
lf = open(log, 'w')

prmtop = app.AmberPrmtopFile(args.prmtop)
gro_file = pref + '.gro'
gro = app.GromacsGroFile(gro_file)
x = mdtraj.load(gro_file)
pdb_file = pref + '.pdb'
x.save_pdb(pdb_file)
pdb = app.PDBFile(pdb_file)
M = Molecule(pdb_file)

pos = pdb.positions

ts = 0.0001
steps = 4

# minimize cpu
for i in range(25):
    if steps > 200:
        break
    if i == 0:
        ref = -1
    else:
        change = V1/ref
        lf.write('dV = ' + str(change))
        lf.write('\n')
        lf.write('----------------------------')
        lf.write('\n')
        if .9 < change < 1.1:
            steps *= 2
        ref = V1
    pos, V1 = rest_sim(ts, pos, steps, 0, None, 0, 'CPU')
    V1 = V1.value_in_unit(unit.kilojoule/unit.mole)
x.xyz[0] = np.array(pos.value_in_unit(unit.nanometer))
x.save_pdb(pref + '_cpu_min_0.pdb')
out_gro = pref + '_cpu_min_0.gro'
x.save_gro(out_gro)

"""
# minimize gpu
steps = 5000
pos, V1 = rest_sim(ts, pos, steps, 0, None, 0, 'CUDA')
V1 = V1.value_in_unit(unit.kilojoule/unit.mole)
x.xyz[0] = np.array(pos.value_in_unit(unit.nanometer))
x.save_pdb(pref + '_gpu_min_0.pdb')

if V1 > 0:
    steps = 1
# eq gpu
pos, V2 = rest_sim(ts, pos, 0, 0, None, steps, 'CUDA')
x.xyz[0] = np.array(pos.value_in_unit(unit.nanometer))
x.save_pdb(pref + '_gpu_eq_0.pdb')

# minimize gpu
steps = 10000
pos, V1 = rest_sim(ts, pos, steps, 0, None, 0, 'CUDA')
x.xyz[0] = np.array(pos.value_in_unit(unit.nanometer))
x.save_pdb(pref + '_gpu_min_1.pdb')
# eq gpu
pos, V2 = rest_sim(ts, pos, 0, 0, None, steps, 'CUDA')
V2 = V2.value_in_unit(unit.kilojoule/unit.mole)
x.xyz[0] = np.array(pos.value_in_unit(unit.nanometer))
x.save_pdb(pref + '_gpu_eq_1.pdb')

x.save_gro(pref + '_min.gro')
"""

end_time = str(datetime.datetime.now().time())

lf.write('----------------------------')
lf.write('\n')
lf.write('FINAL POTENTIAL: ' + str(V1))
lf.write('\n')
lf.write('start: ' + start_time)
lf.write('\n')
lf.write('end: ' + end_time)
lf.write('\n')
lf.write('----------------------------')
Lf.write('\n')
