import re
from subprocess import Popen, PIPE, STDOUT
import zipfile

def natural_sort(l):
    """ From stack overflow: Natural sorting of a list (so 11 comes after 7) """
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

# ex: call_cl(['grep', '-R', 'idk', '.'])
def call_cl(command_lst, pipe_args=[]):
    p = Popen(command_lst, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    out, err = p.communicate(input='\n'.join(pipe_args))
    return out

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
        pth = '/'.join(fn.split('/')[:1])
        zfile.extractall(pth)
