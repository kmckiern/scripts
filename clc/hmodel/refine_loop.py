import sys
from copy import copy
from chimera import runCommand as rc
from chimera import openModels, Molecule
from MultAlignViewer.MAViewer import MAViewer
from MultAlignViewer.AddSeqDialog import AddSeqDialog
from MultAlignViewer import ModellerBase
from align_model import findMAVs

def model(tmpl8, loops):
    # open template structure, and generate sequence
    rc('open ' + tmpl8)
    rc('sequence #0; wait')
    
    # MAV objects
    fmav = findMAVs()
    template = fmav[0]
    temp_seq = template.seqs[0]
    seq_name = tmpl8.split('/')[-1]
    
    # run modeller on alignment
    ModellerBase.model(template, temp_seq, openModels.list(modelTypes=[Molecule]), 
        '10', 1, 1, 0, veryFast=0, loopInfo=('', loops), **kw)

def main():
    # read in args
    _, tmpl8, l0_str, l1_str = sys.argv
    loop0 = tuple([int(i) for i in l0_str.split(',')])
    loops = [loop0]
    if l1_str != 'None':
        loop1 = tuple([int(i) for i in l1_str.split(',')])
        loops.append(loop1)

    # model target to template
    model(tmpl8, loops)
    
if __name__ == '__main__':
    main()
