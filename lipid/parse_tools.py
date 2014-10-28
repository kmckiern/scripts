#!/usr/bin/env python
"""
file of commonly used parsing tools
"""
def get_ls(list_file):
    # get list of atoms to remove.
    ls = []
    for line in open(list_file, 'r').readlines():
        ls.append(line.strip())
    return ls
