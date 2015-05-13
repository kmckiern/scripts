import re
from subprocess import Popen, PIPE, STDOUT
import zipfile

def natural_sort(l):
    """ From stack overflow: Natural sorting of a list (so 11 comes after 7) """
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

def call_cl(command_lst, std=PIPE, pipe_args=[]):
    p = Popen(command_lst, stdin=PIPE, stdout=std, stderr=std, shell=True)
    out, err = p.communicate(input='\n'.join(pipe_args))
    return out, err

def write_out(name, o, e=None):
    oe = open(name, 'w+')
    oe.write(str(o))
    if e != None:
        oe.write(str(e))
    oe.close()

def call_write(command_lst, std=PIPE, pipe_args=[], name):
    o, e = call_cl(command_lst, std, pipe_args)
    write_out(name, o)

"""
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
