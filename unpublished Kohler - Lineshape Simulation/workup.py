'''
First Created 2016/05/06 by Blaise Thompson

Last Edited 2016/05/06 by Blaise Thompson

Contributors: Dan Kohler, Blaise Thompson
'''


### import ####################################################################


import os
import itertools
import collections
from distutils import util

import ConfigParser

import scipy
from scipy.optimize import leastsq
from scipy.interpolate import griddata, interp1d, interp2d, UnivariateSpline

import numpy as np

import NISE
from NISE.lib import pulse
from NISE.lib.misc.__init__ import NISE_path
import numpy as np
# important:  we can't call on scan directly for some reason; use trive to do it
import NISE.lib.measure as m
import NISE.experiments.trive as trive
import NISE.hamiltonians.H0 as H0
import NISE.hamiltonians.params.inhom as inhom

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


### simulation overview 2D frequency (without smear) ##########################


raw_pickle_path = os.path.join(directory, 'simulation_overview.p')
processed_pickle_path = raw_pickle_path

def workup():
    trive.exp.set_coord(trive.d1, 100.)
    trive.exp.set_coord(trive.d2, 200.)
    trive.exp.set_coord(trive.ss, 50.)
    # hamiltonian
    H0.tau_ag = 50.
    H0.tau_2aa = 50.
    H0.tau_2ag = 50.
    H = H0.Omega()
    H.out_group = [[5, 6]]
    # scan
    w1 = trive.w1
    w1.points = np.linspace(5000, 9000, 41)
    w2 = trive.w2
    w2.points = np.linspace(5000, 9000, 41)
    scan = trive.exp.scan(w1, w2, H=H)
    scan.run(mp=False, autosave=False)
    # measure
    mono = m.Mono()
    mono.slitwidth = 120  # wn
    sld = m.SLD()
    measure = m.Measure(scan, mono, sld)
    measure.run()
    # data
    data = wt.data.from_NISE(measure)
    data.transpose()
    data.scale()  # stored as amplitude level
    data.normalize()
    # save
    data.save(raw_pickle_path)
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='simulation overview', 
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)


### simulation overview 2D frequency (with smear) #############################


raw_pickle_path = os.path.join(directory, 'simulation_overview_smeared.p')
processed_pickle_path = raw_pickle_path

def workup():
    trive.exp.set_coord(trive.d1, 100.)
    trive.exp.set_coord(trive.d2, 200.)
    trive.exp.set_coord(trive.ss, 50.)
    # hamiltonian
    H0.tau_ag = 50.
    H0.tau_2aa = 50.
    H0.tau_2ag = 50.
    H = H0.Omega()
    H.out_group = [[5, 6]]
    # scan
    w1 = trive.w1
    w1.points = np.linspace(5000, 9000, 41)
    w2 = trive.w2
    w2.points = np.linspace(5000, 9000, 41)
    scan = trive.exp.scan(w1, w2, H=H)
    scan.run(mp=False, autosave=False)
    # kernal
    # currently no way to NOT do this in multiprocessing
    # additionally this spams my NISE/data folder :-(
    if __name__ == '__main__':
        scan.smear(0, 1, 500, 1e-3, theta=np.pi/4, save=False)
    # measure
    mono = m.Mono()
    mono.slitwidth = 120  # wn
    sld = m.SLD()
    measure = m.Measure(scan, mono, sld)
    measure.run()
    # data
    data = wt.data.from_NISE(measure)
    data.transpose()
    data.scale()  # stored as amplitude level
    data.normalize()
    # save
    data.save(raw_pickle_path)
    return data, data.copy()

# force workup
if False:
    workup()

# automatically process
shared_module.process(key='simulation overview smeared', 
                      workup_method=workup, raw_pickle_path=raw_pickle_path,
                      processed_pickle_path=processed_pickle_path,
                      raw_dictionary=raw_dictionary,
                      processed_dictionary=processed_dictionary)
