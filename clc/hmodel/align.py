def get_class_members(klass):
    ret = dir(klass)
    if hasattr(klass,'__bases__'):
        for base in klass.__bases__:
            ret = ret + get_class_members(base)
    return ret


def uniq( seq ): 
    """ the 'set()' way ( use dict when there's no set ) """
    return list(set(seq))

def get_object_attrs( obj ):
    # code borrowed from the rlcompleter module ( see the code for Completer::attr_matches() )
    ret = dir( obj )
    ## if "__builtins__" in ret:
    ##    ret.remove("__builtins__")

    if hasattr( obj, '__class__'):
        ret.append('__class__')
        ret.extend( get_class_members(obj.__class__) )

        ret = uniq( ret )

    return ret

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
for mav in fmav:
    print '--------------------------------------------'
    print "mav: ", mav.title
    print "attributes: "
    for at in get_object_attrs(mav):
        print at

target = fmav[0]
template = fmav[1]

# align sequences
template.addSeqs(target.seqs)
