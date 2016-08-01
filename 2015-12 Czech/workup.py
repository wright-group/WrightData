'''
First Created 2016/05/05 by Blaise Thompson

Last Edited 2016/08/01 by Blaise Thompson

Contributors: Blaise Thompson, Kyle Czech
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


### movie #####################################################################


raw_pickle_path = os.path.join(directory, 'raw_movie.p')
processed_pickle_path = os.path.join(directory, 'processed_movie.p')

def workup():
    # raw
    data_paths = wt.kit.glob_handler('.dat', folder=os.path.join(directory, 'movie'))
    raw_movie = wt.data.from_COLORS(data_paths, name='MoS2 TrEE Movie')
    raw_movie.save(raw_pickle_path)
    # processed
    processed_movie = raw_movie.copy()
    processed_movie.level('ai0', 'd2', -3)
    processed_movie.smooth([2, 2, 0], channel='ai0')
    processed_movie.scale(channel='ai0', kind='amplitude')
    processed_movie.normalize(channel='ai0')
    processed_movie.save(processed_pickle_path)
    # finish
    return raw_movie, processed_movie

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='movie', 
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)


### absorbance ################################################################


raw_pickle_path = os.path.join(directory, 'absorbance_data.p')
processed_pickle_path = raw_pickle_path

def workup():
    absorbance_path = os.path.join(directory, 'MoS2_TF_III_ebeam_1nm_Mo_onQuartz_T=300K__corrected.txt')
    absorbance_data = wt.data.from_shimadzu(absorbance_path, name='MoS2 thin film absorbance')
    absorbance_data.save(raw_pickle_path)
    return absorbance_data, absorbance_data

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='absorbance', 
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)
