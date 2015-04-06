#!/usr/bin/env python

"""
feed in bilayer trj and generate output file with thermodynamic and structural characterization
example usage:
    >> python /path/to/script/dir/analyze_bilayer.py --trj /path/to/trj.trr --o meh.dat --xyz tp/
"""

import argparse
from subprocess import Popen, PIPE, STDOUT
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('font', family = 'serif')
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams.update({'font.size': 12})

parser = argparse.ArgumentParser(description='full analysis of bilayer trajectories')
parser.add_argument('--prod_dir', type=str, help='path to directory with trj and tpr')
parser.add_argument('--prod_pref', type=str, help='file preface for gmx files')
parser.add_argument('--valid_dir', type=str, help='name of directory for gmx output xyz files')
parser.add_argument('--trim', type=float, help='percent to trim off of the front of the trjs', default=0.0)
parser.add_argument('--o', type=str, help='name of output file')
args = parser.parse_args()

# globals.  not sure if i should bother to make these arguments
gmx_temp = 'struc.xvg'
n_lip = 128.0
n_h2o = 3655.0
v_w = {}
v_w['spc'] = {'323.15': 0.03070, '333.15': 0.03096, '338.15': 0.03109, '353.15': 0.03153}
v_w['tip3pfb'] = {'323.15': 0.03043, '333.15': 0.03061, '338.15': 0.03071, '353.15': 0.03104}
sn_ndx = {}
sn_ndx['1'] = ['a C15', 'a C17', 'a C18', 'a C19', 'a C20', 'a C21', 'a C22', 'a C23', \
           'a C24', 'a C25', 'a C26', 'a C27', 'a C28', 'a C29', 'a C30', 'a C31', \
           'del 0-5', 'q', '']
sn_ndx['2'] = ['a C34', 'a C36', 'a C37', 'a C38', 'a C39', 'a C40', 'a C41', 'a C42', \
           'a C43', 'a C44', 'a C45', 'a C46', 'a C47', 'a C48', 'a C49', 'a C50', \
           'del 0-5', 'q', '']
plt_style = '.k'

def plt_ts(time, func, save_name):
    fig = plt.figure(figsize=(6, 4), dpi=250)  
    ax = plt.subplot(111)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.plot(time, func, plt_style)
    plt.savefig(save_name)

def cl_gmx(command_lst, pipe_args=[]):
    p = Popen(command_lst, stdout=PIPE, stdin=PIPE, stderr=STDOUT)  
    p.communicate(input='\n'.join(pipe_args))

# area per lipid
def get_al(data, vdir):
    t = data[:,0]
    apl = data[:,1]*data[:,2]/(n_lip/2.0)
    plt_ts(t, apl, vdir + 'al.png')
    return np.mean(apl)

# volume per lipid
def get_vl(data, vdir):
    t = data[:,0]
    # lazy for now
    vpl = (data[:,3] - (n_h2o * v_w['tip3pfb']['323.15']))/n_lip
    plt_ts(t, vpl, vdir + 'vl.png')
    return np.mean(vpl)

# deuterium order parameter
def get_scd(sn, pdir, ppref, vdir):
    sn_i = ['make_ndx', '-f', pdir + ppref + '.tpr', '-o', vdir + 'sn' + sn + '.ndx']
    cl_gmx(sn_i, sn_ndx[sn])
    go = ['g_order', '-s', pdir + ppref + '.tpr', '-f', pdir + ppref + '.trr', '-n', \
          vdir + 'sn' + sn + '.ndx', '-od', vdir + 'sn' + sn + '.xvg', '-xvg', 'no']
    cl_gmx(go)
    data = np.genfromtxt(vdir + 'sn' + sn + '.xvg')
    return data

"""# isothermal compressibility modulus
def get_kap():

# x-ray structure factor
def get_fq(): 

# diffusion constant
def get_dl():"""

def get_struc(pdir, ppref, vdir, trim):
    props = ['Volume', 'Box-X', 'Box-Y']
    ge = ['g_energy', '-f', pdir + ppref + '.edr', '-xvg', 'no', '-o', vdir + gmx_temp]
    cl_gmx(ge, props)

    data = np.genfromtxt(vdir + gmx_temp)
    if trim > 0.0:
        d = data[trim*len(data):]
    else:
        d = data
    al = get_al(d, vdir)
    vl = get_vl(d, vdir)
    scd1 = get_scd('1', pdir, ppref, vdir)
    scd2 = get_scd('2', pdir, ppref, vdir)

def main():
    pdir = args.prod_dir
    ppref = args.prod_pref
    vdir = args.valid_dir
    trim = args.trim
    out = args.o

    get_struc(pdir, ppref, vdir, trim)

if __name__ == '__main__':
    main()
