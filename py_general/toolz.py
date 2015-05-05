import re
from subprocess import Popen, PIPE, STDOUT

def natural_sort(l):
    """ From stack overflow: Natural sorting of a list (so 11 comes after 7) """
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

# ex: call_cl(['grep', '-R', 'idk', '.'])
def call_cl(command_lst, pipe_args=[]):
    p = Popen(command_lst, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    p.communicate(input='\n'.join(pipe_args))
