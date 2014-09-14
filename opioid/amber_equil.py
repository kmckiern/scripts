#!/usr/bin/env python

"""
This script is used for generating all of the files and scripts necessary for running two 
cycles of gpcr equilibration on titan using amber12 and ambertools12.

Steps:
read in list of input pdbs in a given directory
creates simulation directories for each pdb

moves restraint file into job dir
generates submission scripts for each dir via this template: equil_c1.sh

once those crash, it does the full equilibration process again from the restart file.
prmtop file needs to be edited for restart file
new submission script generated
"""

import subprocess as subp
import os
import shutil
import argparse

parser = argparse.ArgumentParser(description='Set up serial AMBER equilibration of GPCR system')
parser.add_argument('--root_dir', type=str, help='Input pdb file directory')
parser.add_argument('--ligand', type=str, help='oxy, nal, or dmt')
parser.add_argument('--rst', action='store_true', help='Generate files for restrainted equilibration')
opts = parser.parse_args()

# this should remove inconsistancies in input path names (eg leading and trailing forward slashes)
# has ability to append to a path preface
def format_path(path_pref, append_ele=None):
    # for formatting consistency, parse paths and extension
    pref_list = filter(None, path_pref.split('/'))
    if append_ele:
        pref_list.append(append_ele)
    # combine formatting dir preface and pdb file preface to create job dir
    return '/' + '/'.join(pref_list) + '/'

def gen_job_dir(job_dir, pdb):
    jd = format_path(job_dir, pdb)
    # create dir and return name
    if not os.path.exists(jd):
        os.makedirs(jd)
    return jd

def fill_template(template_file, script_destination, ext, replacement_map):
    """
    template_file: leaprc or bash submission script template
    replacement_map: map from template variables to simulation variables.
    """
    out_file = script_destination + replacement_map['PDB_PREF'] + ext
    with open(template_file, "r") as template:
        lines = template.readlines()
    with open(out_file, "w+") as of:
        for line in lines:
            # i thought maybe '$$' would be an unambiguous split flag.  hence, my
            # templates are written with the variables as follows: $$VAR$$
            l = line.split('$$')
            # if variables are present, replace them using the replacement_map
            if len(l) > 1:
                for region in l:
                    if region in replacement_map:
                        l[l.index(region)] = replacement_map[region]
                of.write(''.join(l))
            else:
                of.write(line)
    return out_file

def gen_leaprc(template_file, lig, pdb_dir, pdb_pref, job_dir, leap_dir):
    # create replacement_map using function inputs
    rm = {'LIG': lig, 'PDB_DIR': pdb_dir, 'PDB_PREF': pdb_pref, 'JOB_DIR': job_dir}
    # replace template variables with input
    leaprc = fill_template(template_file, leap_dir, '_leaprc', rm)
    return leaprc

def gen_equil(template_file, ext, job_name, leap_dir, leaprc, job_dir, script_dir, min_pref, q_0, pdb_pref, nvt_in, npt_in, equil_pref):
    # populate template variable map
    rm = {'JOBNAME': job_name, 'LEAP_DIR': leap_dir, 'LEAPRC': leaprc, 'JOB_DIR': job_dir, 'IN_DIR': script_dir, 'MIN_PREF': min_pref, 'RESTART0': q_0, 'PDB_PREF': pdb_pref, 'H1_PREF': nvt_in, 'H2_PREF': npt_in, 'EQUIL_PREF': equil_pref}
    # generate script from template
    ec = fill_template(template_file, script_dir, ext, rm)
    return ec

def main():
    # formatting paths
    # is it better to have less option typing by mandating directory name conventions?
    root = opts.root_dir
    pd = format_path(root, 'pdbs')
    rd = format_path(root, 'rst_files')
    ld = format_path(root, 'leap_src')
    sd = format_path(root, 'script_templates')
    jr = format_path(root, 'jerbs')

    # pulling relevant files
    pdb_files = [p for p in os.listdir(pd) if '.pdb' in p]
    rst_files = [r for r in os.listdir(rd)]

    # iterate over all pdbs in pdb_dir
    for pdb in pdb_files:
        # pdb pref for output
        f = pdb.split('/')[-1].split('.pdb')[0]
        # format a bit
        f = f.replace('.', '_').replace('-', '_')

        # create job directory
        jd = gen_job_dir(jr, f)

        # setup directory
        if opts.rst:
            # move restraint file
            for rf in rst_files:
                rf_path = rd + rf
                shutil.copy2(rf_path, jd)

        # generate all submission scripts
        # create leaprc
        leap_template = sd + 'leaprc_c1'
        leaprc = gen_leaprc(leap_template, opts.ligand, pd, f, jd, ld)
        # gen first equilibration cycle submission script
        ec1_template = format_paths(td, 'equil_c1.sh')
        ec1_script = gen_leaprc(ec1_template, '_c1.sh', f + '_c1', ld, leaprc, jd, sd*, 'min', pdb_pref + '_out', pdb_pref, 'h1_nvt', 'h2_npt', 'eq1')
        # gen second equil cycle sub script
        ec2_template = format_paths(td, 'equil_c1.sh)
        ec2_script = gen_leaprc(ec2_template, '_c2.sh', f + '_c2', ld, leaprc, jd, sd*, 'min', 'eq1.rst', pdb_pref, 'h1_nvt', 'h2_npt', 'eq1')

if __name__ == '__main__':
    main()
