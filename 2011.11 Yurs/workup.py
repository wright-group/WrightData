'''
First Created 2016/05/05 by Blaise Thompson

Last Edited 2016/05/05 by Blaise Thompson

Contributors: Blaise Thompson
'''

### import ####################################################################


import os
import imp
import collections

import ConfigParser

import WrightTools as wt


### define ####################################################################


directory = os.path.dirname(__file__)
key = os.path.basename(directory)
package_folder = os.path.dirname(directory)
shared_module = imp.load_source('shared', os.path.join(package_folder, 'shared.py'))
google_drive_ini = ConfigParser.SafeConfigParser()
google_drive_ini.read(os.path.join(package_folder, 'google drive.ini'))

raw_dictionary = collections.OrderedDict()
processed_dictionary = collections.OrderedDict()


### download ##################################################################


bypass_download = False

if __name__ == '__main__' and not bypass_download:
    folder_id = google_drive_ini.get('id', key)
    shared_module.download(folder_id, directory)


### absorbance ################################################################


# TODO:


### 2D frequency progression ##################################################


# TODO:


### on diagonal 2D delay ######################################################


def workup():
    folder = os.path.join(directory, 'on diagonal 2D delay')
    files = wt.kit.glob_handler('', folder=folder)
    data = wt.data.from_KENT(files, name='PbSe')
    # fine tune constant labels for clarity
    constant = data.constants[1]
    constant.label_seed = ['1', 'm', '2']
    data.constants = [constant]
    # finish
    data.save(os.path.join(directory, 'on diagonal 2D delay.p'))
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
kwargs = {}
kwargs['key'] = 'on diagonal 2D delay'
kwargs['raw_pickle_path'] = os.path.join(directory, 'on diagonal 2D delay.p')
kwargs['processed_pickle_path'] = os.path.join(directory, 'on diagonal 2D delay.p')
shared_module.process(workup_method=workup, raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary, **kwargs)


### off diagonal 2D delay #####################################################


def workup():
    folder = os.path.join(directory, 'off diagonal 2D delay')
    files = wt.kit.glob_handler('', folder=folder)
    data = wt.data.from_KENT(files, name='PbSe')
    data.save(os.path.join(directory, 'off diagonal 2D delay.p'))
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
kwargs = {}
kwargs['key'] = 'off diagonal 2D delay'
kwargs['raw_pickle_path'] = os.path.join(directory, 'off diagonal 2D delay.p')
kwargs['processed_pickle_path'] = os.path.join(directory, 'off diagonal 2D delay.p')
shared_module.process(workup_method=workup, raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary, **kwargs)


### high resolution 2D delay a ################################################


def workup():
    folder = os.path.join(directory, 'high resolution 2D delay a')
    files = wt.kit.glob_handler('', folder=folder)
    data = wt.data.from_KENT(files, name='PbSe', delay_tolerance=0.001)
    # enforce good behavior
    data = data.split('d1', 3, direction='above')[0]
    # finish
    data.save(os.path.join(directory, 'high resolution 2D delay a.p'))
    return data, data

# force workup
if False:
    workup()

# automatically process
kwargs = {}
kwargs['key'] = 'high resolution 2D delay a'
kwargs['raw_pickle_path'] = os.path.join(directory, 'high resolution 2D delay a.p')
kwargs['processed_pickle_path'] = os.path.join(directory, 'high resolution 2D delay a.p')
shared_module.process(workup_method=workup, raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary, **kwargs)


### high resolution 2D delay b ################################################


def workup():
    folder = os.path.join(directory, 'high resolution 2D delay b')
    files = wt.kit.glob_handler('', folder=folder)
    data = wt.data.from_KENT(files, name='PbSe', delay_tolerance=0.001)
    data.save(os.path.join(directory, 'high resolution 2D delay b.p'))
    return data, data

# force workup
if False:
    workup()

# automatically process
kwargs = {}
kwargs['key'] = 'high resolution 2D delay b'
kwargs['raw_pickle_path'] = os.path.join(directory, 'high resolution 2D delay b.p')
kwargs['processed_pickle_path'] = os.path.join(directory, 'high resolution 2D delay b.p')
shared_module.process(workup_method=workup, raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary, **kwargs)
