''' 
Handle getting of data objects from repository subfolders.
'''


### import ####################################################################


import os
import imp
import sys
import collections

import configparser


### define ####################################################################


d = os.path.dirname(__file__)

folders = [os.path.join(d,o) for o in os.listdir(d) if os.path.isdir(os.path.join(d,o)) and not '.git' in o]


### google drive config #######################################################


# create private google drive file if none exists
google_drive_private_path = os.path.join(d, 'google drive private.ini')
if not os.path.isfile(google_drive_private_path):
    with open(google_drive_private_path, 'w') as f:
        ini = configparser.ConfigParser()
        ini.add_section('id')
        ini.write(f)


### keys ######################################################################


def keys():
    '''
    Get a list of the keys in WrightData.
    '''
    return [os.path.basename(p) for p in folders]


### get #######################################################################


def get(match_string, process=True, check_remote=True):
    '''
    Recieve a dictionary of data objects from a particular publication.
    
    Parameters
    ----------
    match_string : string
        A string matching a publications key uniquely. Does not have to be the
        entire key. 'everything' may be used to workup all publications - in
        this case nothing will be returned.
    process : bool (optional)
        Toggle if recieved data objects are processed or raw. Default is True.
    
    Returns
    -------
    dictionary
        Dictionary of WrightTools data objects for matched publication.

    See Also
    --------
    keys
    '''
    if match_string == 'everything':
        for key in keys():
            print('working up', key)
            get(key)
    else:
        # choose folder
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
        key = os.path.basename(folder)
        # download
        if check_remote:
            from . import shared
            directory = os.path.join(d, key)            
            shared.download(key, directory)
        # get data dictionary
        module = imp.load_source('workup', os.path.join(folder, 'workup.py'))
        if process:
            dictionary = module.processed_dictionary
        else:
            dictionary = module.raw_dictionary
        # finish
        return dictionary
