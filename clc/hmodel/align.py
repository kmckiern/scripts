# borrowed from 
# http://plato.cgl.ucsf.edu/pipermail/chimera-users/2011-July/006573.html
def findMAVs():
    # locate and return the newest instance of MultAlign Viewer
    from MultAlignViewer.MAViewer import MAViewer
    from chimera.extension import manager
    mavs = [inst for inst in manager.instances
                    if isinstance(inst, MAViewer)]
    if not mavs:
        raise AssertionError("No MAV instances!")
    return mavs

fmav = findMAVs()
target = fmav[0]
template = fmav[1]

# align sequences
template.addSeqs(target.seqs)
