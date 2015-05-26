import sys
from copy import copy
import time
import threading

import Pmw, Tkinter
from chimera import runCommand as rc
from chimera import openModels, Molecule
from chimera.extension import manager
from MultAlignViewer.MAViewer import MAViewer
from MultAlignViewer.AddSeqDialog import AddSeqDialog
from MultAlignViewer.prefs import MATRIX, GAP_OPEN, GAP_EXTEND, \
    USE_SS, SS_MIXTURE, SS_SCORES, HELIX_OPEN, STRAND_OPEN, OTHER_OPEN
from MatchMaker.gui import SSParams
from MultAlignViewer import ModellerBase

# borrowed from 
# http://plato.cgl.ucsf.edu/pipermail/chimera-users/2011-July/006573.html
def findMAVs():
    # locate and return list of MultAlign Viewer instances
    mavs = [inst for inst in manager.instances
        if isinstance(inst, MAViewer)]
    if not mavs:
        raise AssertionError("No MAV instances!")
    return mavs

def model(trgt, tmpl8, seq_name):
    # open target sequence
    rc('open ' + trgt)
    # open template structure, and generate sequence
    rc('open ' + tmpl8)
    rc('sequence #0')
    
    # MAV objects
    fmav = findMAVs()
    target, template = fmav
    tar_seq = copy(target.seqs[0])
    temp_seq = copy(template.seqs[0])
    
    # align sequences
    # get template secondary structure matrix
    nb = Pmw.NoteBook()
    structPage = nb.add("From Structure")
    ssParams = SSParams(structPage, template.prefs)
    kw = {'ssMatrix': ssParams.getMatrix()}
    # generalize these vars later.  tired of hacking rn.
    # match target fasta sequence to template pdb sequence
    target.alignSeq(temp_seq, displayName=seq_name, 
        matrix=template.prefs[MATRIX], gapOpenStrand=-18.0, 
        scoreGap=-1, scoreGapOpen=-12, gapOpenHelix=-18.0, 
        gapOpenOther=-6.0, gapChar='.', guideSeqs=None,
        ssFraction=0.3, **kw)
    
    # run modeller on alignment
    ModellerBase.model(target, tar_seq, openModels.list(modelTypes=[Molecule]), 
        '5', 1, 1, 0, veryFast=0, **kw)

# hack due to the fact that python does not wait for modeller to finish
# before trying to write output files
def wait_and_write(seq_name, od):
    while len(openModels.list(modelTypes=[Molecule])) < 2:
        time.sleep(120)
    else:
        rc('write #1 ' + od + '$number_' + seq_name)
        rc('close all')
        rc('stop now')

def main():
    # read in args
    _, trgt, tmpl8, od = sys.argv
    seq_name = tmpl8.split('/')[-1]
    # model target to template
    threading.Thread(target=model(trgt, tmpl8, seq_name)).start()
    # write results
    # threading.Thread(target=wait_and_write(seq_name, od)).start()
    
if __name__ == '__main__':
    main()
