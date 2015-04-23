import numpy as np
import argparse

"""
read in g_energy files for an original configuration, and a perturbation trajectory, 
    and calculates force components for each atom
example usage
    >> python ~/scripts/ff/gmx_f.py --up zero_ge.xvg --p tip3p.xml --h .001
"""

parser = argparse.ArgumentParser(description='get force components')
parser.add_argument('--up', type=str, help='unperturbed g_energy')
parser.add_argument('--p', type=str, help='perturbed g_energy')
parser.add_argument('--dx', type=float, help='perturbation amount (angstroms)')
parser.add_argument('--od', type=str, help='output directory')
args = parser.parse_args()

compare = ['Bond', 'Angle', 'Torsion']#, 'Nb']
relevant = {'Bond': 0, 'Angle': 1, 'Proper Dih.': 2, 'Improper Dih.': 3, 'LJ-14': 4, 'Coulomb-14': 5, 'LJ (SR)': 6, 'Disper. corr.': 7, 'Coulomb (SR)': 8, 'Coul. recip.': 9, 'Position Rest.': 10}
torsions = ['Proper Dih.', 'Improper Dih.']
nonbond = ['LJ-14', 'Coulomb-14', 'LJ (SR)', 'Coulomb (SR)', 'Coul. recip.', 'Position Rest.']

def FD(z, b, a, h):
    return -1.0*np.gradient(np.array([b, z, a]), h)[1]

upf = np.genfromtxt(args.up, skip_header=32)
pf = np.genfromtxt(args.p, skip_header=32)
h = args.dx
od = args.od

for fg in compare:
    f = open(od + '/' + fg + '.dat', 'w')
    if fg == 'Torsions':
        data = None
        for t in torsions:
            data += pf[:, relevant[t] + 1]
    else:
        data = pf[:, relevant[fg] + 1]
    until = len(data)
    for i in range(until / 6):
        # x, y, and z
        for j in [0, 2, 4]:
            if j < 4:
                f.write(str(FD(upf[relevant[fg] + 1], data[j+(i*6)], data[(j+1)+(i*6)], h)) + ', ')
            else:
                f.write(str(FD(upf[relevant[fg] + 1], data[j+(i*6)], data[(j+1)+(i*6)], h)))
        f.write('\n')
    f.close()
