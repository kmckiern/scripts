# this should remove inconsistancies in input path names (eg leading and trailing forward slashes)
# has ability to append to a path preface
def format_path(path_pref, append_ele=None):
    # for formatting consistency, parse paths and extension
    pref_list = filter(None, path_pref.split('/'))
    if append_ele:
        pref_list.append(append_ele)
    # combine formatting dir preface and pdb file preface to create job dir
    return '/' + '/'.join(pref_list) + '/'
