import re
from subprocess import Popen, PIPE, STDOUT
import zipfile

# from stack overflow: natural sorting of a list (so 11 comes after 7)
def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

# call the command line from within python
def call_cl(command_lst, std=PIPE, pipe_args=[]):
    p = Popen(command_lst, stdout=PIPE, stdin=PIPE, stderr=STDOUT, shell=True)
    out, err = p.communicate(input='\n'.join(pipe_args))
    return out, err

# write the output of a CL call to some file
def write_out(name, o, e=None):
    oe = open(name, 'w+')
    oe.write(str(o))
    if e != None:
        oe.write(str(e))
    oe.close()

# combines call and write functions
def call_write(command_lst, name, std=PIPE, pipe_args=[]):
    o, e = call_cl(command_lst)#, std, pipe_args)
    write_out(name, o)

"""
extracts contents of a zip file to an arbitrary out dir
adapated from:
http://stackoverflow.com/questions/9431918/
    extracting-zip-file-contents-to-specific-directory-in-python-2-7
"""
def xtract(fn, dest=None):
    zfile = zipfile.ZipFile(fn)
    if dest != None:
        zfile.extractall(dest)
    else:
        pth = '/'.join(fn.split('/')[:-1])
        zfile.extractall(pth)
