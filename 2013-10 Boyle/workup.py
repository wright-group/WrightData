'''
First Created 2016/06/14 by Blaise Thompson

Last Edited 2016/08/08 by Blaise Thompson

Contributors: Blaise Thompson
'''


### import ####################################################################


import os
import sys
import importlib
import collections

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

# dictionaries to fill
raw_dictionary = collections.OrderedDict()
processed_dictionary = collections.OrderedDict()


### download ##################################################################


bypass_download = False

if __name__ == '__main__' and not bypass_download:
    shared_module.download(key, directory)


### full 2D TRSF ##############################################################


raw_pickle_path = os.path.join(directory, 'full 2D TRSF.p')
processed_pickle_path = raw_pickle_path

def workup():
    folder = os.path.join(directory, 'full 2D TRSF')
    filepaths = wt.kit.glob_handler('', folder=folder)
    data = wt.data.from_KENT(filepaths, name='TRSF', ignore=['d1', 'd2', 'wm'])
    data.save(raw_pickle_path)
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='full 2D TRSF',
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)
