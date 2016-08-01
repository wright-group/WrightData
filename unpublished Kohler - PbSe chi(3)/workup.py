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
    path = os.path.join(directory, 'absorbance', 'fsb 19.1.txt')
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
    path = os.path.join(directory, 'absorbance', 'DK-2014.09.17-PbSe 025.4.a.txt')
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


### OPA1 power curve - movie ##################################################


raw_pickle_path = os.path.join(directory, 'OPA1 power curve - movie.p')
processed_pickle_path = raw_pickle_path

def workup():
    path = os.path.join(directory, 'movie calibration', 'power curves', 'opa1 power smoothness (1) 4000 ms wait [1x51] Freq.dat')
    data = wt.data.from_COLORS(path, name='OPA1 power curve')
    data.save(raw_pickle_path)
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='OPA1 power curve movie', 
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)



### OPA2 power curve - movie ##################################################


raw_pickle_path = os.path.join(directory, 'OPA2 power curve - movie.p')
processed_pickle_path = raw_pickle_path

def workup():
    path = os.path.join(directory, 'movie calibration', 'power curves', 'opa2 power smoothness (1) 4000 ms wait [1x51] Freq.dat')
    data = wt.data.from_COLORS(path, name='OPA2 power curve')
    data.save(raw_pickle_path)
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='OPA2 power curve movie', 
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)


### batch A movie #############################################################


raw_pickle_path = os.path.join(directory, 'batch A movie raw.p')
processed_pickle_path = os.path.join(directory, 'batch A movie.p')

def workup():
    # get raw
    filepaths = wt.kit.glob_handler('.dat', folder=os.path.join(directory, 'CMDS movie A'))
    raw = wt.data.from_COLORS(filepaths, name='Batch A Movie')
    raw.save(raw_pickle_path)
    # process
    # TODO:
    processed = raw.copy()
    processed.save(processed_pickle_path)
    # finish
    return (raw, processed)

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='batch A movie', 
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)


### batch B movie #############################################################


# TODO:


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
