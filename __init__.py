''' 
Handle getting of data objects from repository subfolders.
'''

### import ####################################################################


import os
import imp
import collections


### get method ################################################################


def get(match_string, process=True):
    # get all subfolder paths
    d = os.path.dirname(__file__)
    folders = [os.path.join(d,o) for o in os.listdir(d) if os.path.isdir(os.path.join(d,o)) and not '.git' in o]
    # choose particular folder
    matching_folders = []
    for f in folders:
        if match_string in f:
            matching_folders.append(f)
    if len(matching_folders) > 1:
        raise LookupError('\'' + match_string + '\' matches multile data folders')
    elif len(matching_folders) == 0:
        raise LookupError('\'' + match_string + '\' does not match any data folders')
    else:
        folder = matching_folders[0]
    # get data dictionary
    module = imp.load_source('workup', os.path.join(folder, 'workup.py'))
    if process:
        dictionary = module.processed_dictionary
    else:
        dictionary = module.raw_dictionary
    # finish
    return dictionary
