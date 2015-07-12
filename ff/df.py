import numpy as np
import argparse
from numpy import linalg as LA
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

"""
compare force components between md packages
example usage
    >> python ~/scripts/ff/df.py --omm omm/HarmonicBondForce.dat 
        --gmx gmx/Bond.dat --name diff/bond.dat
"""

parser = argparse.ArgumentParser(description='get force differences')
parser.add_argument('--omm', type=str, help='file of openmm forces')
parser.add_argument('--gmx', type=str, help='file of gmx forces')
parser.add_argument('--name', type=str, help='name of output file')
args = parser.parse_args()

o = np.genfromtxt(args.omm, delimiter=",")
g = np.genfromtxt(args.gmx, delimiter=",")

dF =  LA.norm(o, axis=1)-LA.norm(g, axis=1)
np.savetxt(args.name + '.dat', dF)

hist, bin_edges = np.histogram(dF, bins = 100)
plt.bar(bin_edges[:-1], hist, width = 1)
plt.ylim([0, 20])
plt.title(args.name)
plt.xlabel('|F|_o - |F|_g / kJ/(nm mol)')
plt.savefig(args.name + '.png')
