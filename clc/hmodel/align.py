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
    # locate and return list of MultAlign Viewer instances
    mavs = [inst for inst in manager.instances
        if isinstance(inst, MAViewer)]
    if not mavs:
        raise AssertionError("No MAV instances!")
    return mavs

fmav = findMAVs()
target = fmav[0]
template = fmav[1]
seq = copy(template.seqs[0])
seq_name = sys.argv[2].split('/')[-1]

# get template secondary structure matrix
nb = Pmw.NoteBook()
structPage = nb.add("From Structure")
ssParams = SSParams(structPage, template.prefs)
kw = {}
kw['ssMatrix'] = ssParams.getMatrix()

# align sequences
# generalize these vars later.  tired of hacking rn.
target.alignSeq(seq, displayName=seq_name, 
    matrix=template.prefs[MATRIX], gapOpenStrand=-18.0, 
    scoreGap=-1, scoreGapOpen=-12, gapOpenHelix=-18.0, 
    gapOpenOther=-6.0, gapChar='.', guideSeqs=None,
    ssFraction=0.3, **kw)
