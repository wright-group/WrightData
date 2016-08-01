''' 
Handle getting of data objects from repository subfolders.
'''

### import ####################################################################


import os
import imp
import sys
import collections

try:
    import configparser as ConfigParser  # python 3
except ImportError:
    import ConfigParser as ConfigParser  # python 2

from . import shared


### define ####################################################################


d = os.path.dirname(__file__)

folders = [os.path.join(d,o) for o in os.listdir(d) if os.path.isdir(os.path.join(d,o)) and not '.git' in o]
if sys.version[0] == '3':
    google_drive_ini = ConfigParser.ConfigParser()
else:
    google_drive_ini = ConfigParser.SafeConfigParser()
google_drive_ini.read(os.path.join(d, 'google drive.ini'))


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
            directory = os.path.join(d, key)
            folder_id = google_drive_ini.get('id', key)
            shared.download(folder_id, directory)
        # get data dictionary
        module = imp.load_source('workup', os.path.join(folder, 'workup.py'))
        if process:
            dictionary = module.processed_dictionary
        else:
            dictionary = module.raw_dictionary
        # finish
        return dictionary
