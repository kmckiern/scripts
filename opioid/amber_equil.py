#!/usr/bin/env python

"""
This script is for generation of production trajectories.
Borrowed heavily from Morgan.

Steps:
read in list of input pdbs in a given directory
creates simulation directories for each pdb

moves restraint file into job dir
generates submission scripts for each dir via template

once those crash, it does the full equilibration process again from the restart file.

updates prmtop file just before the start of the second minimization.
    via:
    1. ambpdb to convert rst file to pdb
    2. leap pdb to generate new inpcrd and prmtop files

example usage:
    >> ./amber_equil.py --root_dir path/to/root/dir --ligand oxy --rst --e1

ex 2:
    >> for i in {61..70}; do python gen_amb.py --root_dir path/meh --ligand nal --script_name long --template eq_gpu_cont.sh --cyc $i --pdb_sub_dir all; done
"""

import subprocess as subp
import os
import shutil
import argparse
import math
import copy

parser = argparse.ArgumentParser(description='Set up serial AMBER equilibration of GPCR system')
parser.add_argument('--root_dir', type=str, help='Input root file directory')
parser.add_argument('--pdb_sub_dir', type=str, help='Input pdb file directory')
parser.add_argument('--jerb_sub_dir', type=str, help='Jerb sub directory', default=None)
parser.add_argument('--ligand', type=str, help='oxy, nal, or dmt')
parser.add_argument('--script_name', type=str, help='Name of run script')
parser.add_argument('--template', type=str, help='Name of template script')
parser.add_argument('--min', type=str, help='cpu or gpu')
parser.add_argument('--cyc', type=str, help='number of equilibration cycle')
parser.add_argument('--rst', action='store_true', help='Generate files for restrained equilibration')
parser.add_argument('--leap', action='store_true', help='Generate files for restrained equilibration')
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

# create dir and return name
def gen_job_dir(job_dir, pdb):
    jd = format_path(job_dir, pdb)
    if not os.path.exists(jd):
        os.makedirs(jd)
    return jd

def fill_template(template_file, script_destination, ext, replacement_map, multi=None):
    """
    template_file: leaprc or bash submission script template
    replacement_map: map from template variables to simulation variables.
    """
    out_file = script_destination + opts.script_name + ext
    if multi:
        num_pdbs = len(replacement_map['PDB_PREF'])
    else:
        num_pdbs = 1
    replacement_map['NUM_JOBS'] = str(num_pdbs)
    with open(template_file, "r") as template:
        lines = template.readlines()
    with open(out_file, "w") as of:
        for line in lines:
            # i thought maybe '$$' would be an unambiguous split flag.  hence, my
            # templates are written with the variables as follows: $$VAR$$
            l = line.split('$$')
            # if variables are present, replace them using the replacement_map
            if len(l) > 1:
                if l[1] == 'VAR_CMD':
                    # remove flag
                    l.pop(0)
                    l.pop(0)
                    l_bk = copy.deepcopy(l)
                    # write an amber command for all pdbs given in the replacement_map
                    if num_pdbs > 1:
                        for i in range(num_pdbs):
                            for region in l:
                                if region in replacement_map:
                                    mapped = replacement_map[region]
                                    if len(mapped) == 1 or isinstance(mapped, str):
                                        l[l.index(region)] = mapped
                                    else:
                                        l[l.index(region)] = mapped[i]
                            of.write(''.join(l))
                            l = copy.deepcopy(l_bk)
                    else:
                            for region in l:
                                if region in replacement_map:
                                    r_map = replacement_map[region]
                                    if isinstance(r_map, list):
                                        l[l.index(region)] = r_map[0]
                                    else:
                                        l[l.index(region)] = r_map
                            of.write(''.join(l))
                else:
                    for region in l:
                        if region in replacement_map:
                            l[l.index(region)] = replacement_map[region]
                    of.write(''.join(l))
            else:
                of.write(line)
    return out_file

def gen_leaprc(template_file, lig, pdb_dir, pdb_pref, job_dir, leap_dir, rc_suffix, cyc, box_list=None):
    # create replacement_map using function inputs
    c_val = 'c' + cyc
    c_m1 = 'c' + str(int(cyc)-1)
    rm = {'LIG': lig.lower(), 'PDB_DIR': pdb_dir, 'PDB_PREF': pdb_pref, 'JOB_DIR': job_dir, 'CYCLE_NUM': c_val, 'CYCLE_MONE': c_m1}
    if box_list:
        # get rst box dimensions.
        rm['BOXX'], rm['BOXY'], rm['BOXZ'] = box_list
    # replace template variables with input
    leaprc = fill_template(template_file, leap_dir, rc_suffix, rm)
    return leaprc

def gen_equil(template_file, leap_dir, leaprc, job_dir, script_dir, min_pref, pdb_pref, nvt_in, npt_in, equil_pref, gen_dir, cyc, prod_pref):
    c_val = 'c' + cyc
    c_m1 = 'c' + str(int(cyc)-1)
    ext = '_' + c_val + '.sh'
    job_name = opts.script_name
    # populate template variable map
    rm = {'FULL_JOB_NAME': job_name, 'LEAP_DIR': leap_dir, 'LEAPRC': leaprc, 'JOB_DIR': job_dir, 'IN_DIR': script_dir, 'MIN_PREF': min_pref, 'PDB_PREF': pdb_pref, 'H1_PREF': nvt_in, 'H2_PREF': npt_in, 'EQUIL_PREF': equil_pref, 'CYCLE_NUM': c_val, 'CYCLE_MONE': c_m1,  'PROD_PREF': prod_pref}
    # generate script from template
    ec = fill_template(template_file, gen_dir, ext, rm, pdb_pref)
    return ec

def main():
    # formatting paths
    # is it better to have less option typing by mandating directory name conventions?
    root = opts.root_dir
    sd = opts.pdb_sub_dir
    if sd is not None:
        sub_d = sd.strip()
        pd = format_path(root, 'pdbs/' + sub_d)
    else:
        pd = format_path(root, 'pdbs/')
    cycle = opts.cyc
    pcycle = str(int(opts.cyc)-1)

    rd = format_path(root, 'rst_files')
    ld = format_path(root, 'leap_src')
    sd = format_path(root, 'static_scripts')
    od = format_path(root, 'run_scripts')
    if opts.jerb_sub_dir == None:
        jr = format_path(root, 'jerbs')
    else:
        jr = format_path(root, 'jerbs/' + opts.jerb_sub_dir)

    # pulling relevant files
    pdb_files = [p for p in os.listdir(pd) if '.pdb' in p]
    restraint_files = [r for r in os.listdir(rd)]

    num_files = len(pdb_files)
    # iterate over all pdbs in pdb_dir
    for j, pdb in enumerate(pdb_files):
        # pdb pref for output
        f = pdb.split('/')[-1].split('.pdb')[0] 
        # format a bit
        f = f.replace('.', '_').replace('-', '_')

        # setup directory
        if opts.rst:
            # create job directory
            jd = gen_job_dir(jr, f)
            # move restraint file
            for rf in restraint_files:
                rf_path = rd + rf
                shutil.copy2(rf_path, jd)
        else:
            jd = format_path(jr, f)

        if opts.leap:
            # generate all submission scripts
            # gen subsequent equil cycle sub script
            if int(cycle) > 1:
                er = 'equil_c' + pcycle + '.rst'
                h2r = 'h2_npt_c' + pcycle + '.rst'
                h1r = 'h1_nvt_c' + pcycle + '.rst'
                # get the final coordinates of the last restart box.
                # determine last non-empty restart file.
                rst_files = [outf for outf in os.listdir(jd) if '.rst' in outf and os.path.getsize(jd + outf) > 0]
                if er in rst_files:
                    rst_file = er
                elif h2r in rst_files:
                    rst_file = h2r
                elif h1r in rst_files:
                    rst_file = h1r
                else:
                    rst_file = None
                    leaprc = None
                if rst_file:
                    box_dims = subp.Popen(['tail', '-n', '1', jd + rst_file], stdout=subp.PIPE).communicate()[0].split()[0:3]
                    b_d = [str(int(math.ceil(float(i)))) for i in box_dims]
                    # create leaprc
                    leap_template = sd + 'leaprc_c' + cycle
                    leaprc = gen_leaprc(leap_template, opts.ligand, pd, f, jd, ld, '_' + f + '_leaprc_c' + cycle, b_d)
            else:
                # gen first equilibration cycle submission script
                # create leaprc
                leap_template = sd + 'leaprc'
                leaprc = gen_leaprc(leap_template, opts.ligand, pd, f, jd, ld, '_' + f + '_leaprc', cycle)
        else:
            leaprc = None

        # record data
        if j == 0:
            lrc = [leaprc]
            jds = [jd]
            prefs = [f]
        else:
            lrc.append(leaprc)
            jds.append(jd)
            prefs.append(f)

    ec_script = gen_equil(sd + opts.template, ld, lrc, jds, sd, 'min', prefs, 'h1_nvt', 'h2_npt', 'equil', od, cycle, 'prod')

if __name__ == '__main__':
    main()
