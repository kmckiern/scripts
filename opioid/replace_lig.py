#!/usr/bin/env python
"""
Reads in a template pdb and ligand pdb, and adds the ligand pdb to the template.
"""
import os
import argparse

parser = argparse.ArgumentParser(description='lig-template pdb generation')
parser.add_argument('--t_pdb', type=str, help='template pdb')
parser.add_argument('--lig_pdb', type=str, help='ligand pdb')
parser.add_argument('--lig_name', type=str, help='ligand residue name')
parser.add_argument('--final_pdb', type=str, help='final pdb name (extension optional)')
args = parser.parse_args()

# Take ligand coordinates, insert into template pdb, and generate filled-in pdb.
def fill_template(template_pdb, lig_coords, resname, generated_pdb):
    # consistency
    if len(generated_pdb.split('.')) > 1:
        out_file = generated_pdb
    else:
        out_file = generated_pdb + '.pdb'
    # fill template
    with open(template_pdb, "r") as template:
        lines = template.readlines()
    with open(out_file, "w") as of:
        for line in lines:
            l = line.split()
            if resname in l:
                l[5], l[6], l[7] = lig_coords[l[2]]
                # hella hacky.  pdb column formatting
                # there is definitely a better way to do this
                if len(l[2]) == 2:
                    space_12 = '   '
                if len(l[2]) == 3:
                    space_12 = '  '
                if l[5].startswith('-'):
                    space_45 = '     '
                else:
                    space_45 = '      '
                if l[6].startswith('-'):
                    space_56 = '  '
                else:
                    space_56 = '   '
                if l[7].startswith('-'):
                    space_67 = '  '
                else:
                    space_67 = '   '
                of.write('%s  %s%s%s %s   %s%s%s%s%s%s%s  %s  %s\n' % (l[0], l[1], space_12, l[2], l[3], l[4], space_45, l[5], space_56, l[6], space_67, l[7], l[8], l[9]))
            else:
                of.write(line)
    return out_file

# Build dictionary of ligand coordinate data
def parse_lig(lig_pdb, resname):
    atm_coords = {}
    for line in open(lig_pdb, 'r'):
        l = line.split()
        if resname in l:
            atm_coords[l[2]] = l[-6:-3]
    return atm_coords

def main():
    res = args.lig_name.upper()

    # read in ligand pdb, store lines
    res_qs = parse_lig(args.lig_pdb, res)

    # parse template pdb, and replace lig coordinates
    of = fill_template(args.t_pdb, res_qs, res, args.final_pdb)

if __name__ == '__main__':
    main()
