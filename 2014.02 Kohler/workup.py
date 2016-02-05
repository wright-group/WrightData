### import ####################################################################


import os
import itertools
import collections
from distutils import util

import scipy
from scipy.optimize import leastsq
from scipy.interpolate import griddata, interp1d, interp2d, UnivariateSpline

import numpy as np

import WrightTools as wt


### define ####################################################################


directory = os.path.dirname(__file__)

raw_dictionary = collections.OrderedDict()
processed_dictionary = collections.OrderedDict()


### absorbance A ##############################################################


# raw and processed are identical in this case
out_path = os.path.join(directory, 'absorbance_A.p')

force_workup = False

def workup():
    path = os.path.join(directory, 'DK-10.13.12-PbSe 011a.txt')
    data = wt.data.from_JASCO(path, name='A')
    data.convert('wn')
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
dict_name = 'absorbance A'
raw_dictionary[dict_name] = raw_data
processed_dictionary[dict_name] = processed_data


### absorbance B ##############################################################


# raw and processed are identical in this case
out_path = os.path.join(directory, 'absorbance_B.p')

force_workup = False

def workup():
    path = os.path.join(directory, 'DK-04.23.12-PbSe 010.a.txt')
    data = wt.data.from_JASCO(path, name='B')
    data.convert('wn')
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
dict_name = 'absorbance B'
raw_dictionary[dict_name] = raw_data
processed_dictionary[dict_name] = processed_data


### absorbance C ##############################################################


# raw and processed are identical in this case
out_path = os.path.join(directory, 'absorbance_C.p')

force_workup = False

def workup():
    path = os.path.join(directory, 'DK-03.20.12-PbSe 005.a.txt')
    data = wt.data.from_JASCO(path, name='C')
    data.convert('wn')
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
dict_name = 'absorbance C'
raw_dictionary[dict_name] = raw_data
processed_dictionary[dict_name] = processed_data


### 2D frequency progression (sample A) #######################################


# raw and processed are identical in this case
out_path = os.path.join(directory, '2D_frequencies_A.p')

force_workup = False

def workup():
    folder = os.path.join(directory, '2D frequencies')
    files = wt.kit.glob_handler('.dat', folder=folder)
    data = wt.data.from_COLORS(files, name='2D frequencies A', ignore=['num', 'd1'])
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
dict_name = '2D frequencies A'
raw_dictionary[dict_name] = raw_data
processed_dictionary[dict_name] = processed_data


### 2D delays (sample B) ######################################################

# raw and processed are identical in this case
out_path = os.path.join(directory, '2D_delays_B.p')

force_workup = False

def workup():
    folder = os.path.join(directory, '2D delays')
    files = wt.kit.glob_handler('.dat', folder=folder)
    data = wt.data.from_COLORS(files, name='2D delays B')
    data.d1.points *= -1  # invert d1s
    data.flip('d1')
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
dict_name = '2D delays B'
raw_dictionary[dict_name] = raw_data
processed_dictionary[dict_name] = processed_data


### diagonal wigner (sample C) ################################################


# raw and processed are identical in this case
out_path = os.path.join(directory, 'diagonal_wigner_C.p')

force_workup = False

def workup():
    folder = os.path.join(directory, 'diagonal wigner')
    files = wt.kit.glob_handler('.dat', folder=folder)
    data = wt.data.from_COLORS(files, name='diagonal wigner C', ignore=['num', 'd1', 'wm'])
    data.transpose()
    data = data.split('d2', 200)[0]
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
dict_name = 'diagonal wigner C'
raw_dictionary[dict_name] = raw_data
processed_dictionary[dict_name] = processed_data


### w2 wigners (sample C) ####################################################


# raw and processed are identical in this case
out_path = os.path.join(directory, 'w2_wigners_C.p')

force_workup = False

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
dict_name = 'w2 wigners C'
raw_dictionary[dict_name] = raw_data
processed_dictionary[dict_name] = processed_data
