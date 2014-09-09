#!/bin/env python

"""
Generates psi input file from molecule xyz geometry file.
Example usage:
>> for i in *xyz; do PREF=$(echo $i | awk -F '.' '{print $1}'); python gen_psi.py --df $i --theory scf --basis 6-31G* --of hf/${PREF}.dat; done
"""

from optparse import OptionParser

parser = OptionParser()
parser.add_option('--df', help='Input coordinate (geometry) file')
parser.add_option('--theory', help='mp2, wB97x-D, B97-D3, b3lyp, scf (HF)')
parser.add_option('--basis', help='heavy-aug-cc-pvTZ, heavy-aug-cc-pvDZ, 6-311++G**, 6-31G*')
parser.add_option('--of', help='Output file')
(opts, args) = parser.parse_args()

def gen_specz(b, t):
    if t == 'mp2':
        return 'set {\n    basis    %s\n    guess    sad\n    print    1\n    freeze_core    true\n}\n\nenergy(\'%s\')\n' % (b, t)
    else:
        return 'set {\n    basis    %s\n    guess    sad\n    print    1\n}\n\nenergy(\'%s\')\n' % (b, t)

def main():
    options = opts.__dict__
    
    mol_xyz = options['df']
    theory = options['theory']
    basis = options['basis']
    of = options['of']

    with open(of, 'w') as file:
        file.write('memory 12 gb\n\n')
        file.write('molecule DPPC {\n0 1\n')
        file.writelines(open(mol_xyz, 'r').readlines()[2:])
        file.write('no_reorient\n}\n\n')
        file.write(gen_specz(basis, theory))

if __name__ == '__main__':
    main()
