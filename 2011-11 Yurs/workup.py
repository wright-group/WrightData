'''
First Created 2016/05/05 by Blaise Thompson

Last Edited 2016/08/01 by Blaise Thompson

Contributors: Blaise Thompson
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


bypass_download = True

if __name__ == '__main__' and not bypass_download:
    folder_id = google_drive_ini.get('id', key)
    shared_module.download(folder_id, directory)


### absorbance ################################################################


# TODO:


### 2D frequency progression ##################################################


# TODO:


### on diagonal 2D delay ######################################################


raw_pickle_path = os.path.join(directory, 'on diagonal 2D delay.p')
processed_pickle_path = raw_pickle_path

def workup():
    folder = os.path.join(directory, 'on diagonal 2D delay')
    files = wt.kit.glob_handler('', folder=folder)
    data = wt.data.from_KENT(files, name='PbSe')
    # fine tune constant labels for clarity
    constant = data.constants[1]
    constant.label_seed = ['1', 'm', '2']
    data.constants = [constant]
    # finish
    data.save(raw_pickle_path)
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='on diagonal 2D delay', 
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)


### off diagonal 2D delay #####################################################


raw_pickle_path = os.path.join(directory, 'off diagonal 2D delay.p')
processed_pickle_path = raw_pickle_path

def workup():
    folder = os.path.join(directory, 'off diagonal 2D delay')
    files = wt.kit.glob_handler('', folder=folder)
    data = wt.data.from_KENT(files, name='PbSe')
    data.save(raw_pickle_path)
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='off diagonal 2D delay', 
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)


### high resolution 2D delay a ################################################


raw_pickle_path = os.path.join(directory, 'high resolution 2D delay a.p')
processed_pickle_path = raw_pickle_path

def workup():
    folder = os.path.join(directory, 'high resolution 2D delay a')
    files = wt.kit.glob_handler('', folder=folder)
    data = wt.data.from_KENT(files, name='PbSe', delay_tolerance=0.001)
    # enforce good behavior
    data = data.split('d1', 3, direction='above')[0]
    # finish
    data.save(raw_pickle_path)
    return data, data

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='high resolution 2D delay a', 
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)


### high resolution 2D delay b ################################################


raw_pickle_path = os.path.join(directory, 'high resolution 2D delay b.p')
processed_pickle_path = raw_pickle_path

def workup():
    folder = os.path.join(directory, 'high resolution 2D delay b')
    files = wt.kit.glob_handler('', folder=folder)
    data = wt.data.from_KENT(files, name='PbSe', delay_tolerance=0.001)
    data.save(raw_pickle_path)
    return data, data

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='high resolution 2D delay b',
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)
