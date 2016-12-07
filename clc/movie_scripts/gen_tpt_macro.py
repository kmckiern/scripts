from chimera import runCommand as rc 
from chimera import openModels, MSMSModel
from Movie.gui import MovieDialog
from chimera.extension import manager

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

tpath = [140, 105, 120, 121, 133,  64, 113,  43,  60]

def findMDMovie():
    mdms = [inst for inst in manager.instances if isinstance(inst,
MovieDialog)]
    if not mdms:
        raise AssertionError('No MD Movie Instances!')
    return mdms[-1]


def do_it(x, y, lbl):
    fn = 'full_trj'
    meta = open('metafile', 'w')
    meta_input = '''pdb\nsingle\n%s.pdb\n''' % fn
    meta.write(meta_input)
    meta.close()
    
    rc('open movie:metafile')
    
    rc('open /home/kmckiern/Downloads/smoothMD.py')
    
    rc('sel #0')
    clr = clrs[tpath[0]]
    rc('color 0.4,.4,.4 sel')

    rc('open axes.bild')
    
    rc('open style_base.cmd')
    rc('open align.cmd')

    x = str(x)
    y = str(y)

    rc('turn x ' + x)
    rc('turn y ' + y)
    
    rc('open idk_surf.cmd')
    
    rc('~sel')

    mdm = findMDMovie()
    
    rc('movie record supersample 4')
    
    for tndx, j in enumerate(tpath):
        # if j == 133:
        #     rc('turn y 2 -90')
        #     rc('wait 45')
        #     rc('freeze')            

        clr = clrs[j]
        rc('sel #0:85,128 #0:Cl-')
        rc('color ' + clr + ' ~sel')
        rc('color byhetero #0')
    
        start = int((tndx*50.0) + 1)
        end = int((tndx*50.0) + 51)
        for x in range(start, end):
            mdm.LoadFrame(x)
            rc('wait 1')
    
    rc('movie stop')
    rc('movie encode roundtrip true format mp4 output movies/' + str(fn) + '_macrocl_' + lbl + '.mp4') 
    
    rc('sel')
    rc('close sel')
    
do_it(-90, 90, 'xy')
do_it(-90, -90, 'xny')
do_it(0, 0, 'zz')
do_it(-90, 0, 'xry')
do_it(-180, 0, 'nz')

rc('stop')
