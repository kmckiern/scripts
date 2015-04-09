import re
import subprocess as subp

def natural_sort(l):
    """ From stack overflow: Natural sorting of a list (so 11 comes after 7) """
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

def call_cml(cmd):
    cmd_ls = cmd.split()
    return subp.Popen(cmd_ls, stdout=subp.PIPE).communicate()[0]
