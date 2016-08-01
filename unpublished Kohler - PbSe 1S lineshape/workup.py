'''
First Created 2016/05/06 by Blaise Thompson

Last Edited 2016/08/01 by Blaise Thompson

Contributors: Dan Kohler, Blaise Thompson
'''


### import ####################################################################


import os
import sys
import importlib
import collections

try:
    import configparser as ConfigParser  # python 3
except ImportError:
    import ConfigParser as ConfigParser  # python 2

import WrightTools as wt


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


### download ##################################################################


bypass_download = False

if __name__ == '__main__' and not bypass_download:
    folder_id = google_drive_ini.get('id', key)
    shared_module.download(folder_id, directory)


### batch A absorbance ########################################################


raw_pickle_path = os.path.join(directory, 'batch A absorbance.p')
processed_pickle_path = raw_pickle_path

def workup():
    path = os.path.join(directory, 'dilution study', 'fsb 19.1.txt')
    data = wt.data.from_JASCO(path, name='Batch A')
    data.save(raw_pickle_path)
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='batch A absorbance',
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)


### batch B absorbance ########################################################


raw_pickle_path = os.path.join(directory, 'batch B absorbance.p')
processed_pickle_path = raw_pickle_path

def workup():
    path = os.path.join(directory, 'DK-2014.09.17-PbSe 025.4.a.txt')
    data = wt.data.from_JASCO(path, name='Batch B')
    data.save(raw_pickle_path)
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='batch B absorbance',
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)


### batch A 2D delay ##########################################################


raw_pickle_path = os.path.join(directory, 'batch A 2D delay raw.p')
processed_pickle_path = os.path.join(directory, 'batch A 2D delay processed.p')

def workup():
    path = os.path.join(directory, 'TrEE (1) 500 ms wait [21x21] Delay.dat')
    raw = wt.data.from_COLORS(path, name='Batch A')
    raw.save(raw_pickle_path)
    processed = raw.copy()
    processed.level(0, 'd1', -3)
    processed.level(0, 'd2', -3)
    processed.scale(kind='amplitude')
    processed.normalize()
    processed.save(processed_pickle_path)
    return raw, processed

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='batch A 2D delay', 
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)
