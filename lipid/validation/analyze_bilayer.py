#!/usr/bin/env python

"""
feed in bilayer trj and generate output file with thermodynamic and structural characterization
example usage:
    >> python /path/to/script/dir/analyze_bilayer.py --trj /path/to/trj.trr --o meh.dat --xyz tp/
"""

import argparse
from subprocess import Popen, PIPE, STDOUT
import numpy as np
import scipy.integrate as integr8
import matplotlib.pyplot as plt
import matplotlib as mpl

parser = argparse.ArgumentParser(description='full analysis of bilayer trajectories')
parser.add_argument('--prod_dir', type=str, help='path to directory with trj and tpr')
parser.add_argument('--prod_pref', type=str, help='file preface for gmx files')
parser.add_argument('--valid_dir', type=str, help='name of directory for gmx output xyz files')
parser.add_argument('--trim', type=float, help='percent to trim off of the front of the trjs', default=0.0)
parser.add_argument('--temp', type=float, help='simulation temperature', default=323.15)
parser.add_argument('--wm', type=str, help='water model', default='tip3pfb')
parser.add_argument('--o', type=str, help='name of output file')
args = parser.parse_args()

# globals.  not sure if i should bother to make these arguments
gmx_temp = 'struc.xvg'
n_lip = 128.0
n_h2o = 3655.0
k = 1.3806488e-23

temp = args.temp
wm = args.wm
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

mpl.rc('font', family = 'serif')
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams.update({'font.size': 12})
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

def plt_save(t, ts, name):
    plt_ts(t, ts, name + '.png')
    np.savetxt(name + '.dat', zip(t, ts), delimiter='\t')

# area per lipid
def get_al(data, vdir):
    t = data[:,0]
    apl = data[:,1]*data[:,2]/(n_lip/2.0)
    plt_save(t, apl, vdir + 'al')
    return apl

# volume per lipid
def get_vl(data, vdir):
    t = data[:,0]
    # lazy for now
    vpl = (data[:,3] - (n_h2o * v_w[wm][str(temp)]))/n_lip
    plt_save(t, vpl, vdir + 'vl')
    return np.mean(vpl)

# deuterium order parameter
def get_scd(sn, pdir, ppref, vdir):
    sn_i = ['make_ndx', '-f', pdir + ppref + '.tpr', '-o', vdir + 'sn' + sn + '.ndx']
    cl_gmx(sn_i, sn_ndx[sn])
    go = ['g_order', '-s', pdir + ppref + '.tpr', '-f', pdir + ppref + '.trr', '-n', \
          vdir + 'sn' + sn + '.ndx', '-od', vdir + 'sn' + sn + '.xvg', '-xvg', 'no']
    cl_gmx(go)
    data = np.genfromtxt(vdir + 'sn' + sn + '.xvg')
    plt_ts(data[:,0], data[:,1], vdir + 'sn' + sn + '.png')
    return data

# isothermal compressibility modulus
def get_kap(al_ts):
    avg = np.average(al_ts)
    fluc = np.average(al_ts**2) - avg**2
    kap = 1e3 * (2 * k * temp / n_lip) * (avg * 1.0 / fluc)
    return kap
def kap_ts(time, al_ts, vdir):
    kap_ts = [0]
    for t in range(len(al_ts)):
        if t == 0:
            continue
        else:
            kap_ts.append(get_kap(al_ts[:t]))
    kap = np.array(kap_ts)
    # flucs are too crazy to plot beginning
    s0 = len(kap)*.25
    plt_save(time[s0:], kap[s0:], vdir + 'kap')
    return np.mean(kap)

# x-ray structure factor
def calc_fq_real(rho_dehyd, q):
    return integr8.trapz(rho_dehyd[:,1] * np.cos(q * rho_dehyd[:,0]), rho_dehyd[:,0])
def get_fq(): 
    # gen EDP

    # parse file name.
    pref = df.split('.')[0]
    # system EDP
    data = np.genfromtxt(df)
    # convert to angstroms and shift zero to bilayer center
    max_vals = (np.max(data[:, 0]), np.max(data[:, 1]))
    new_max = 10 * max_vals[0] / 2
    z_shift = np.linspace(-new_max, new_max, num = data.shape[0])
    # center EDP data
    rescale_sys = np.zeros(data.shape)
    mid_bin = data.shape[0] / 2
    rescale_sys[0:mid_bin, 1] = data[mid_bin:, 1]
    rescale_sys[mid_bin:, 1] = data[0:mid_bin, 1]
    rescale_sys[:, 0] = z_shift
    # subtract bulk water electron density
    dehy_sys = np.zeros(rescale_sys.shape)
    dehy_sys[:,0] = rescale_sys[:,0]
    dehy_sys[:,1] = (rescale_sys[:,1] - rescale_sys[:,1][-1])/1000 # Angstroms.
    # FT the data
    q_vec = np.linspace(0, 1, num = 1000)
    f_q = np.zeros((len(q_vec),2))
    for q in range(1000):
        f_q[q, 0] = q_vec[q]
        f_q[q, 1] = calc_fq_real(dehy_sys, q_vec[q])
    f_q[:,1] = np.abs(f_q[:,1])
    np.savetxt(pref + '_fq.dat', f_q, delimiter='\t')
    np.savetxt(pref + '_edp.dat', rescale_sys, delimiter='\t')
    return f_q, rescale_sys

"""# diffusion constant
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
    al_ts = get_al(d, vdir)
    al = np.mean(al_ts)
    #vl = get_vl(d, vdir)
    #scd1 = get_scd('1', pdir, ppref, vdir)
    #scd2 = get_scd('2', pdir, ppref, vdir)
    kap_a = kap_ts(d[:,0], al_ts, vdir)

    print al
    #print vl
    #print scd1
    #print scd2

def main():
    pdir = args.prod_dir
    ppref = args.prod_pref
    vdir = args.valid_dir
    trim = args.trim
    out = args.o

    get_struc(pdir, ppref, vdir, trim)

if __name__ == '__main__':
    main()
