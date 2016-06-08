import os
import glob
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages
import numpy as np

# change to folder with data files
out_rd = '/Users/kerimckiernan/Downloads/dock/latest/moar_follow'
pdbs = glob.glob(out_rd + '/mf*.mol2')

rc("windowsize 1050 1000")

# loop through the files, opening and processing each in turn
for ndx, fn in enumerate(pdbs):
    rc("open " + fn)
    rc("sel #" + str(ndx) + ":@/serialNumber=15 or serialNumber=12")
    rc("align sel")
    rc("turn y 90 center view")
    rc("sel #" + str(ndx))
    rc("center sel")

    rc("sel #" + str(ndx) + ":@C")
    rc("color black sel")

rc("preset apply publication 1")
rc("light mode ambient")
rc("rep bs")
rc("light mode ambient")
rc("setattr m stickScale .32")
rc("setattr m ballScale .12")

rc("scale 1.4")

for ndx, fn in enumerate(pdbs):
    rc("sel")
    rc("~display sel")

    fls = fn.split('/')[-1].split('_')
    if len(fls) == 2:
        lbl = fls[0]
    if len(fls) == 3:
        lbl = '_'.join(fls[:2])
    else:
        lbl = fls[0].split('.')[0]

    rc("sel #" + str(ndx))
    rc("display sel")

    rc("~select")

    rc("2dlabels create lab" + str(ndx) + " text '" + lbl + "' style bold color black size 86 xpos .1 ypos .2 visibility show")

    rc("copy file imgs/" + lbl + ".png width 400")

    rc("2dlabels delete lab" + str(ndx))
