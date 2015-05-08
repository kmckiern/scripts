#!/usr/bin/env python

"""
steps:
0. throw out pdb atoms w d_rl > 30 A (prevent seg faulting), save mol2
1. remove hydrogens from trimmed pdb, save mol2
2. use chimera to generate receptor surface file
3. generate docking spheres
4. edit number of spheres (cut those > 10 A from ligands)
5. generate box/grid
6. run dock program
"""

import argparse
import os
import mdtraj
import sys

parser = argparse.ArgumentParser(description='automate clc docking')
parser.add_argument('--rec', type=str, help='receptor file')
parser.add_argument('--lig', type=str, help='ligand pdb')
parser.add_argument('--dist', type=float, help='dist from lig for spheres (/nm)')
parser.add_argument('--sr', type=str, help='script root', default='/home/kmckiern/scripts/')
parser.add_argument('--trim_rec', action='store_true', help='trim pdb?')
args = parser.parse_args()

sr = args.sr
sys.path.insert(0, sr + 'py_general/')
from toolz import call_cl

def trim_rl(r, l, dist):
    rt = ['python', sr + 'clc/dock/trim_pdb.py', '--rec', r, '--lig', l, '--d', str(dist)]
    call_cl(rt)

def call_chimera(args):
    command = ['chimera', '--nogui', '--script', args]
    call_cl(command)

def gen_templ8(tf, fill, dest):
    keys = fill.keys()
    fld = dest + '/' + tf.split('/')[-1]
    with open(tf, "r") as template:
        lines = template.readlines()
    with open(fld, 'w+') as of:
        for line in lines:
            l = line.split('$')
            if len(l) > 1:
                temp = []
                for i in l:
                    if i in keys:
                        temp.append(fill[i])
                    else:
                        temp.append(i)
                line = ''.join(temp)
            of.write(line)

def write_out(name, o, e):
    oe = open(name + '_oe.dat', 'w+')
    oe.write(o)
    oe.write(e)
    oe.close()

def main():
    cwd = os.getcwd()
    rec = args.rec
    lig = args.lig
    d = args.dist
    pref = rec.split('.')[0]
    ligpref = lig.split('.')[0]
    d6bin = '/home/kmckiern/src/dock6/bin/'

    # if pdb is gt 15000 atoms, prob need to trim.
    if args.trim_rec:
        trim_rl(rec, lig, d*3.0)
        rec = pref + '_trim.pdb'
        pref = rec.split('.')[0]

    # generate ligand and receptor (h and no h) mol2 files
    gm2 = sr + 'clc/dock/gen_mol2.py '
    call_chimera(gm2 + lig)
    call_chimera(gm2 + rec)
    call_chimera(gm2 + rec + ' --noH')

    # write receptor surface file
    rnoH = pref + '_noH.mol2'
    pnoH = rnoH.split('.')[0]
    gsurf = sr + 'clc/dock/gen_surf.py '
    call_chimera(gsurf + rnoH)

    # fill in all of the dock templates
    tmap = {'SURFACE': pnoH, 'SPHERE': pref, 'RECEPTOR': pref}
    template_dir = sr + 'clc/dock/templates/'
    tfs = [template_dir + t for t in os.listdir(template_dir)]
    for tf in tfs:
        gen_templ8(tf, tmap, cwd)

    # spheres
    try:
        os.remove(cwd + '/OUTSPH')
    except OSError:
        pass
    gs = [d6bin + 'sphgen', '-i', 'INSPH', '-o', 'OUTSPH']
    call_cl(gs)
    # trim sphere file
    ss = [d6bin + 'sphere_selector', pref + '_big.sph', ligpref + '.mol2', str(d*10.0)]
    call_cl(ss)
    # gen box and grid
    gb = [d6bin + 'showbox']
    gbin = open('showbox.in').read().splitlines() 
    call_cl(gb, gbin)
    ggrid = [d6bin + 'grid', '-i', 'grid.in']
    grid_out, grid_err = call_cl(ggrid)
    write_out('grid', grid_out, grid_err):
    # dock some ligands!!
    dock = [d6bin + 'dock6', '-i', 'dock.in']
    dock_out, dock_err = call_cl(dock)
    write_out('dock', dock_out, dock_err):

if __name__ == '__main__':
    main()
