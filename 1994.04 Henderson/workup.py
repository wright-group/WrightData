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
