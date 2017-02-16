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
    out, err = p.communicate(input=b'\n'.join(pipe_args))
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

# via http://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays
def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.
    """
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype
    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)
    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out
