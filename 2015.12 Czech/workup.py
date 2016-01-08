### import ####################################################################


import os
import collections

import WrightTools as wt


### define ####################################################################


directory = os.path.dirname(__file__)

raw_dictionary = collections.OrderedDict()
processed_dictionary = collections.OrderedDict()


### movie #####################################################################


def workup_movie():
    # raw
    data_paths = wt.kit.glob_handler('.dat', folder=os.path.join(directory, 'movie'))
    raw_movie = wt.data.from_COLORS(data_paths, name='MoS2 TrEE Movie')
    raw_movie.save(os.path.join(directory, 'raw_movie.p'))
    # processed
    processed_movie = raw_movie.copy()
    processed_movie.level('ai0', 'd2', -3)
    processed_movie.smooth([2, 2, 0], channel='ai0')
    processed_movie.scale(channel='ai0', kind='amplitude')
    processed_movie.normalize(channel='ai0')
    processed_movie.save(os.path.join(directory, 'processed_movie.p'))
    # finish
    return raw_movie, processed_movie

# get from pickle or create
if os.path.isfile(os.path.join(directory, 'raw_movie.p')):
    raw_movie = wt.data.from_pickle(os.path.join(directory, 'raw_movie.p'), verbose=False)
    processed_movie = wt.data.from_pickle(os.path.join(directory, 'processed_movie.p'), verbose=False)
else:
    raw_movie, processed_movie = workup_movie()

# check version
if not raw_movie.__version__.split('.')[0] == wt.__version__.split('.')[0]:
    raw_movie, processed_movie = workup_movie()

# add to dictionaries
raw_dictionary['movie'] = raw_movie
processed_dictionary['movie'] = processed_movie


### absorbance ################################################################


def workup_absorbance():
    absorbance_path = os.path.join(directory, 'MoS2_TF_III_ebeam_1nm_Mo_onQuartz_T=300K__corrected.txt')
    absorbance_data = wt.data.from_shimadzu(absorbance_path, name='MoS2 thin film absorbance')
    absorbance_data.save(os.path.join(directory, 'absorbance_data.p'))
    return absorbance_data
    
# get from pickle or create
if os.path.isfile(os.path.join(directory, 'absorbance_data.p')):
    absorbance_data = wt.data.from_pickle(os.path.join(directory, 'absorbance_data.p'), verbose=False)
else:
    absorbance_data = workup_absorbance()
    
# check version
if not absorbance_data.__version__.split('.')[0] == wt.__version__.split('.')[0]:
    absorbance_data = workup_absorbance()
    
# add to dictionaries
raw_dictionary['absorbance'] = absorbance_data
processed_dictionary['absorbance'] = absorbance_data
