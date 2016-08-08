'''
First Created 2016/05/05 by Blaise Thompson

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


### absorbance A ##############################################################


raw_pickle_path = os.path.join(directory, 'absorbance_A.p')
processed_pickle_path = raw_pickle_path

def workup():
    path = os.path.join(directory, 'DK-10.13.12-PbSe 011a.txt')
    data = wt.data.from_JASCO(path, name='A')
    data.convert('wn')
    data.save(raw_pickle_path)
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='absorbance A',
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)


### absorbance B ##############################################################


raw_pickle_path = os.path.join(directory, 'absorbance_B.p')
processed_pickle_path = raw_pickle_path

def workup():
    path = os.path.join(directory, 'DK-04.23.12-PbSe 010.a.txt')
    data = wt.data.from_JASCO(path, name='B')
    data.convert('wn')
    data.save(raw_pickle_path)
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='absorbance B',
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)


### absorbance C ##############################################################


raw_pickle_path = os.path.join(directory, 'absorbance_C.p')
processed_pickle_path = raw_pickle_path

def workup():
    path = os.path.join(directory, 'DK-03.20.12-PbSe 005.a.txt')
    data = wt.data.from_JASCO(path, name='C')
    data.convert('wn')
    data.save(raw_pickle_path)
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='absorbance C', 
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)


### 2D frequency progression (sample A) #######################################


raw_pickle_path = os.path.join(directory, '2D_frequencies_A.p')
processed_pickle_path = raw_pickle_path

def workup():
    folder = os.path.join(directory, '2D frequencies')
    files = wt.kit.glob_handler('.dat', folder=folder)
    data = wt.data.from_COLORS(files, name='2D frequencies A', ignore=['num', 'd1'])
    data.save(raw_pickle_path)
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='2D frequencies A', 
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)


### 2D delays (sample B) ######################################################


raw_pickle_path = os.path.join(directory, '2D_delays_B.p')
processed_pickle_path = raw_pickle_path

def workup():
    folder = os.path.join(directory, '2D delays')
    files = wt.kit.glob_handler('.dat', folder=folder)
    data = wt.data.from_COLORS(files, name='2D delays B')
    data.d1.points *= -1  # invert d1s
    data.flip('d1')
    data.save(raw_pickle_path)
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='2D delays B', 
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)


### diagonal wigner (sample C) ################################################


raw_pickle_path = os.path.join(directory, 'diagonal_wigner_C.p')
processed_pickle_path = raw_pickle_path

def workup():
    folder = os.path.join(directory, 'diagonal wigner')
    files = wt.kit.glob_handler('.dat', folder=folder)
    data = wt.data.from_COLORS(files, name='diagonal wigner C', ignore=['num', 'd1', 'wm'])
    data.transpose()
    data = data.split('d2', 200)[0]
    data.save(raw_pickle_path)
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='diagonal wigner C', 
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)


### w2 wigners (sample C) ####################################################


raw_pickle_path = os.path.join(directory, 'w2_wigners_C.p')
processed_pickle_path = raw_pickle_path

def workup():
    superfolder = os.path.join(directory, 'w2 wigners')
    subfolders = ['w1 6600', 'w1 6900', 'w1 7200']
    files = []
    for folder in [os.path.join(superfolder, f) for f in subfolders]:
        dats = wt.kit.glob_handler('.dat', folder=folder)
        files += dats
    data = wt.data.from_COLORS(files, name='w2 wigners C', ignore=['num', 'd1', 'wm'])
    data.transpose()
    data = data.split('d2', 200)[0]    
    data.save(raw_pickle_path)
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='w2 wigners C',
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)
