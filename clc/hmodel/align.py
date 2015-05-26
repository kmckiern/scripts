import sys
from copy import copy
import Pmw, Tkinter
from chimera.extension import manager
from MultAlignViewer.MAViewer import MAViewer
from MultAlignViewer.AddSeqDialog import AddSeqDialog
from MultAlignViewer.prefs import MATRIX, GAP_OPEN, GAP_EXTEND, \
    USE_SS, SS_MIXTURE, SS_SCORES, HELIX_OPEN, STRAND_OPEN, OTHER_OPEN
from MatchMaker.gui import SSParams

# borrowed from 
# http://plato.cgl.ucsf.edu/pipermail/chimera-users/2011-July/006573.html
def findMAVs():
    # locate and return the newest instance of MultAlign Viewer
    mavs = [inst for inst in manager.instances
                    if isinstance(inst, MAViewer)]
    if not mavs:
        raise AssertionError("No MAV instances!")
    return mavs

def get_frame():
    ss = [inst for inst in manager.instances
                    if isinstance(inst, AddSeqDialog)]
    return ss

fmav = findMAVs()
target = fmav[0]
template = fmav[1]
seq_name = sys.argv[2].split('/')[-1]

# get template secondary structure matrix
nb = Pmw.NoteBook()
structPage = nb.add("From Structure")
ssParams = SSParams(structPage, template.prefs)
kw['ssMatrix'] = ssParams.getMatrix()

# align sequences
target.alignSeq(copy(template.seqs[0]), displayName=seq_name, matrix=template.prefs[MATRIX], gapOpenStrand=-18.0, scoreGap=-1, scoreGapOpen=-12, gapOpenHelix=-18.0, gapOpenOther=-6.0, gapChar='.', guideSeqs=None, kw)
