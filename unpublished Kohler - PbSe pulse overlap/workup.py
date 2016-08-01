'''
First Created 2016/08/01 by Blaise Thompson

Last Edited 2016/08/01 by Blaise Thompson

Contributors: Blaise Thompson
'''


### import ####################################################################


import os
import sys
import importlib
import collections
from distutils import util

import numpy as np

try:
    import configparser as ConfigParser  # python 3
except ImportError:
    import ConfigParser as ConfigParser  # python 2


### define ####################################################################


# paths
directory = os.path.dirname(__file__)
key = os.path.basename(directory)
package_folder = os.path.dirname(directory)

# shared module
spec = importlib.util.spec_from_file_location('shared', os.path.join(package_folder, 'shared.py'))
shared_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(shared_module)

# google drive
if sys.version[0] == '2':
    google_drive_ini = ConfigParser.SafeConfigParser()
else:
    google_drive_ini = ConfigParser.ConfigParser()
google_drive_ini.read(os.path.join(package_folder, 'google drive.ini'))

# dictionaries to fill
raw_dictionary = collections.OrderedDict()
processed_dictionary = collections.OrderedDict()


### batch A absorbance ########################################################


# raw and processed are identical in this case
out_path = os.path.join(directory, 'batch A absorbance.p')

force_workup = False

def workup():
    path = os.path.join(directory, 'dilution study', 'fsb 19.1.txt')
    data = wt.data.from_JASCO(path, name='Batch A')
    data.save(out_path)
    return data, data.copy()

# get from pickle or create
if os.path.isfile(out_path) and not force_workup:
    raw_data = wt.data.from_pickle(out_path, verbose=False)
    processed_data = raw_data.copy()
else:
    raw_data, processed_data = workup()

# check version
if raw_data is None:
    pass
elif not raw_data.__version__.split('.')[0] == wt.__version__.split('.')[0]:
    raw_data, processed_data = workup()

# add to dictionaries
dict_name = 'batch A absorbance'
raw_dictionary[dict_name] = raw_data
processed_dictionary[dict_name] = processed_data


### batch B absorbance ########################################################


# raw and processed are identical in this case
out_path = os.path.join(directory, 'batch B absorbance.p')

force_workup = False

def workup():
    path = os.path.join(directory, 'DK-2014.09.17-PbSe 025.4.a.txt')
    data = wt.data.from_JASCO(path, name='Batch B')
    data.save(out_path)
    return data, data.copy()

# get from pickle or create
if os.path.isfile(out_path) and not force_workup:
    raw_data = wt.data.from_pickle(out_path, verbose=False)
    processed_data = raw_data.copy()
else:
    raw_data, processed_data = workup()

# check version
if raw_data is None:
    pass
elif not raw_data.__version__.split('.')[0] == wt.__version__.split('.')[0]:
    raw_data, processed_data = workup()

# add to dictionaries
dict_name = 'batch B absorbance'
raw_dictionary[dict_name] = raw_data
processed_dictionary[dict_name] = processed_data


### batch A 2D delay ##########################################################


raw_out_path = os.path.join(directory, 'batch A 2D delay raw.p')
processed_out_path = os.path.join(directory, 'batch A 2D delay processed.p')

force_workup = False

def workup():
    path = os.path.join(directory, 'TrEE (1) 500 ms wait [21x21] Delay.dat')
    raw = wt.data.from_COLORS(path, name='Batch A')
    raw.save(raw_out_path)
    processed = raw.copy()
    processed.level(0, 'd1', -3)
    processed.level(0, 'd2', -3)
    processed.scale(kind='amplitude')
    processed.normalize()
    processed.save(processed_out_path)
    return raw, processed

# get from pickle or create
if os.path.isfile(os.path.join(directory, raw_out_path)) and not force_workup:
    raw_data = wt.data.from_pickle(os.path.join(directory, raw_out_path), verbose=False)
    processed_data = wt.data.from_pickle(os.path.join(directory, processed_out_path), verbose=False)
else:
    raw_data, processed_data = workup()

# check version
if raw_data is None:
    pass
elif not raw_data.__version__.split('.')[0] == wt.__version__.split('.')[0]:
    raw_data, processed_data = workup()

# add to dictionaries
dict_name = 'batch A 2D delay raw'
raw_dictionary[dict_name] = raw_data
processed_dictionary[dict_name] = processed_data



### 2D frequencies vs OD ######################################################

if False:
    raw_out_path = os.path.join(directory, '2D frequency vs OD raw.p')
    processed_out_path = os.path.join(directory, '2D frequency vs OD processed.p')
    
    force_workup = True
    
    def workup():
        
        return data, data.copy()
    
    # get from pickle or create
    if os.path.isfile(os.path.join(directory, raw_out_path)) and not force_workup:
        raw_data = wt.data.from_pickle(os.path.join(directory, raw_out_path), verbose=False)
        processed_data = wt.data.from_pickle(os.path.join(directory, processed_out_path), verbose=False)
    else:
        raw_data, processed_data = workup()
    
    # check version
    if raw_data is None:
        pass
    elif not raw_data.__version__.split('.')[0] == wt.__version__.split('.')[0]:
        raw_data, processed_data = workup()
    
    # add to dictionaries
    raw_dictionary['2D frequency vs OD'] = raw_data
    processed_dictionary['2D frequency vs OD'] = processed_data
