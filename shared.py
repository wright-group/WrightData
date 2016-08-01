''' 
Methods shared accross WrightData.
'''

### import ####################################################################


import os

import WrightTools as wt


### google drive download #####################################################


def download(folder_id, directory):
    try:
        drive = wt.google_drive.Drive()
        ids = drive.list_folder(folder_id)
        for fileid in ids:
            drive.download(fileid, directory=directory)
    except Exception as inst:
        print(inst)


### process ###################################################################


def process(key, workup_method, raw_dictionary, processed_dictionary,
            raw_pickle_path, processed_pickle_path):
    # get from pickle or create
    is_raw_pickle = os.path.isfile(raw_pickle_path)
    is_processed_pickle = os.path.isfile(processed_pickle_path)
    if is_raw_pickle and is_processed_pickle:
        try:
            raw = wt.data.from_pickle(raw_pickle_path, verbose=False)
            processed = wt.data.from_pickle(processed_pickle_path, verbose=False)
        except UnicodeDecodeError:  # pickles written in wrong version of python, probably
            raw, processed = workup_method()
    else:
        raw, processed = workup_method()
    # check version (should be identical, check only processed)
    if not processed.__version__.split('.')[0] == wt.__version__.split('.')[0]:
        raw_movie, processed_movie = workup_method()
    # add to dictionaries
    raw_dictionary[key] = raw
    processed_dictionary[key] = processed
