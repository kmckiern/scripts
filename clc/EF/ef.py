from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import argparse
import parmed as pmd
import os

parser = argparse.ArgumentParser(description='ef sim')
parser.add_argument('--prmtop', type=str, help='prm', default='3org_l.prmtop')
parser.add_argument('--inpcrd', type=str, help='inp', default='3org_ic_0.inpcrd')
args = parser.parse_args()

parm = pmd.load_file(args.prmtop, args.inpcrd)
parmed_pdb = args.inpcrd.replace('.inpcrd', '_p.pdb')
if not os.path.exists(parmed_pdb):
    parm.write_pdb(parmed_pdb)
pdb = app.PDBFile(parmed_pdb)
system = parm.createSystem()

integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds,
    2.0*unit.femtoseconds)
integrator.setConstraintTolerance(0.00001)

platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}

simulation = app.Simulation(parm.topology, system, integrator, platform,
    properties)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
simulation.step(1)

box = simulation.context.getState().getPeriodicBoxVectors()
V = -.12*unit.volt
L_z = box[2].value_in_unit(unit.angstrom)[-1]*unit.angstrom
force_v = mm.CustomExternalForce("V*q_i*((z + (L_z/2))/L_z)")
force_v.addGlobalParameter("V", V)
force_v.addGlobalParameter("L_z", L_z)
force_v.addPerParticleParameter("q_i")

nonbonded = [f for f in system.getForces()][3]
for i in range(parm.topology.getNumAtoms()):
    charge = nonbonded.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
    force_v.addParticle(i, [charge])
system.addForce(force_v)

# 10 ns
nsteps     = 5000000
# every 200 ps
chkp_intv  = 100000
# every 100 ps
report_int = int(chkp_intv * .5)

ipc = args.inpcrd.split('.')[0] + '_0'
trj_f = ipc + '.dcd'
chk_f = ipc + '.chk'

simulation.reporters.append(app.DCDReporter(trj_f, report_int))
simulation.reporters.append(app.StateDataReporter(stdout, report_int, step=True, 
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True, 
    speed=True, totalSteps=nsteps, separator='\t'))

num_intervals = int(nsteps * (1.0/ (chkp_intv*1.0)))
for i in range(num_intervals):
    simulation.step(chkp_intv)

    with open(chk_f, 'wb') as f:
        f.write(simulation.context.createCheckpoint())
