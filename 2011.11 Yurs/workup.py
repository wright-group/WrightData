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


### download ##################################################################


try:
    drive = wt.google_drive.Drive()
    ids = drive.list_folder('0BzJTClorMBuwcldzTTA3cXpJVkU')
    for fileid in ids:
        drive.download(fileid, directory=directory)
except:
    pass


### absorbance ################################################################


# TODO:


### 2D frequency progression ##################################################


# TODO:


### on diagonal 2D delay ######################################################


# raw and processed are identical in this case
out_path = os.path.join(directory, 'on diagonal 2D delay.p')

force_workup = False

def workup():
    folder = os.path.join(directory, 'on diagonal 2D delay')
    files = wt.kit.glob_handler('', folder=folder)
    data = wt.data.from_KENT(files, name='PbSe')
    # fine tune constants for clarity
    constant = data.constants[1]
    constant.label_seed = ['1', 'm', '2',]
    data.constants = [constant]
    # finish
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
dict_name = 'on diagonal 2D delay'
raw_dictionary[dict_name] = raw_data
processed_dictionary[dict_name] = processed_data


### off diagonal 2D delay #####################################################


# raw and processed are identical in this case
out_path = os.path.join(directory, 'off diagonal 2D delay.p')

force_workup = False

def workup():
    folder = os.path.join(directory, 'off diagonal 2D delay')
    files = wt.kit.glob_handler('', folder=folder)
    data = wt.data.from_KENT(files, name='PbSe')
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
dict_name = 'off diagonal 2D delay'
raw_dictionary[dict_name] = raw_data
processed_dictionary[dict_name] = processed_data


### high resolution 2D delay a ################################################


# raw and processed are identical in this case
out_path = os.path.join(directory, 'high resolution 2D delay a.p')

force_workup = False

def workup():
    folder = os.path.join(directory, 'high resolution 2D delay a')
    files = wt.kit.glob_handler('', folder=folder)
    data = wt.data.from_KENT(files, name='PbSe', delay_tolerance=0.001)
    # enforce good behavior
    data = data.split('d1', 3, direction='above')[0]
    # finish
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
dict_name = 'high resolution 2D delay a'
raw_dictionary[dict_name] = raw_data
processed_dictionary[dict_name] = processed_data


### high resolution 2D delay b ################################################


# raw and processed are identical in this case
out_path = os.path.join(directory, 'high resolution 2D delay b.p')

force_workup = False

def workup():
    folder = os.path.join(directory, 'high resolution 2D delay b')
    files = wt.kit.glob_handler('', folder=folder)
    data = wt.data.from_KENT(files, name='PbSe', delay_tolerance=0.001)
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
dict_name = 'high resolution 2D delay b'
raw_dictionary[dict_name] = raw_data
processed_dictionary[dict_name] = processed_data
