### import ####################################################################


import os
import collections

import WrightTools as wt


### define ####################################################################


directory = os.path.dirname(__file__)

raw_dictionary = collections.OrderedDict()
processed_dictionary = collections.OrderedDict()


### download ##################################################################


try:
    drive = wt.google_drive.Drive()
    ids = drive.list_folder('0BzJTClorMBuwT1ptZ2JMalJPdUk')
    for fileid in ids:
        drive.download(fileid, directory=directory)
except:
    pass
