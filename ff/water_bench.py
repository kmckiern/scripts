from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import argparse

"""
ex usage
    >> python ~/scripts/ff/water_bench.py --pff whatever.xml --wff tip3p.xml 
        --pdb pdb/pdb0.pdb --t 300.0 --od lig/
"""

parser = argparse.ArgumentParser(description='run benchmark sims')
parser.add_argument('--pff', type=str, help='protein FF')
parser.add_argument('--wff', type=str, help='water FF')
parser.add_argument('--pdb', type=str, help='pdb')
parser.add_argument('--t', type=float, help='temperature')
parser.add_argument('--od', type=str, help='output directory')
args = parser.parse_args()

tK = args.t
pdb = app.PDBFile(args.pdb)
forcefield = app.ForceField(args.pff + '.xml', args.wff + '.xml')

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.CutoffPeriodic, 
    nonbondedCutoff=0.6*unit.nanometers, constraints=app.HBonds, rigidWater=True)
integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 
    2.0*unit.femtoseconds)
integrator.setConstraintTolerance(0.00001)
system.addForce(mm.MonteCarloBarostat(1*unit.atmospheres, tK*unit.kelvin, 25))

platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}
simulation = app.Simulation(pdb.topology, system, integrator, platform,
    properties)
simulation.context.setPositions(pdb.positions)

print('Minimizing...')
simulation.minimizeEnergy(maxIterations=5000)

simulation.context.setVelocitiesToTemperature(tK*unit.kelvin)
print('Equilibrating...')
simulation.step(100000)

simulation.reporters.append(app.DCDReporter(args.od + args.pff + '_' + 
    args.wff + '.dcd', 25000))
simulation.reporters.append(app.StateDataReporter(stdout, 25000, step=True,
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
    speed=True, totalSteps=12500000, separator='\t'))

print('Running Production...')
simulation.step(12500000)
print('Done!')
