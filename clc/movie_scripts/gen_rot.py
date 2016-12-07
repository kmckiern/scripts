from chimera import runCommand as rc # use 'rc' as shorthand for runCommand

clrs = {
    64: '0.260,0.650,0.800', 
    113: '0.165,0.545,0.743', 
    133: '0.363,0.733,0.808', 
    121: '0.479,0.798,0.769', 
    120: '0.601,0.844,0.728', 
    105: '0.705,0.884,0.730', 
    43: '0.074,0.452,0.696', 
    140: '0.798,0.921,0.772', 
    60: '0.031,0.354,0.616'
}

fn = 'rotate_example.pdb'
rc('open ' + fn)
rc('sel #0')
clr = clrs[64]
rc('color ' + clr + ' sel')

ion_fn = 'rotate_ex_ions.pdb'
rc('open ' + fn)

rc('open style_base.cmd')
rc('open align.cmd')
rc('open rotations.cmd')
