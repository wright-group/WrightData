'''
First Created 2016-05-06 by Blaise Thompson

Last Edited 2016-08-08 by Blaise Thompson

Contributors: Dan Kohler, Blaise Thompson
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
shared_module_path = os.path.join(package_folder, 'shared.py')
spec = importlib.util.spec_from_file_location('shared', shared_module_path)
shared_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(shared_module)

# dictionaries to fill
raw_dictionary = collections.OrderedDict()
processed_dictionary = collections.OrderedDict()


### download ##################################################################


bypass_download = False

if __name__ == '__main__' and not bypass_download:
    shared_module.download(key, directory)
