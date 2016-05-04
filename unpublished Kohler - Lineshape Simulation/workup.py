### import ####################################################################


import os
import itertools
import collections
from distutils import util

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

raw_dictionary = collections.OrderedDict()
processed_dictionary = collections.OrderedDict()


### download ##################################################################


try:
    drive = wt.google_drive.Drive()
    ids = drive.list_folder('0BzJTClorMBuwWkpRbDlLZHdmSXc')
    for fileid in ids:
        drive.download(fileid, directory=directory)
except Exception as inst:
    print inst


### simulation overview 2D frequency (without smear) ##########################


# raw and processed are identical in this case
out_path = os.path.join(directory, 'simulation_overview.p')

force_workup = False

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
    data.save(out_path)
    return data, data.copy()

# get from pickle or create
if os.path.isfile(os.path.join(directory, out_path)) and not force_workup:
    raw_data = wt.data.from_pickle(os.path.join(directory, out_path), verbose=False)
    processed_data = raw_data.copy()
else:
    raw_data, processed_data = workup()

# check version
if raw_data is None:
    pass
elif not raw_data.__version__.split('.')[0] == wt.__version__.split('.')[0]:
    raw_data, processed_data = workup()

# add to dictionaries
raw_dictionary['simulation overview'] = raw_data
processed_dictionary['simulation overview'] = processed_data


### simulation overview 2D frequency (with smear) #############################


# raw and processed are identical in this case
out_path = os.path.join(directory, 'simulation_overview_smeared.p')

force_workup = False

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
    data.save(out_path)
    return data, data.copy()

# get from pickle or create
if os.path.isfile(os.path.join(directory, out_path)) and not force_workup:
    raw_data = wt.data.from_pickle(os.path.join(directory, out_path), verbose=False)
    processed_data = raw_data.copy()
else:
    raw_data, processed_data = workup()

# check version
if raw_data is None:
    pass
elif not raw_data.__version__.split('.')[0] == wt.__version__.split('.')[0]:
    raw_data, processed_data = workup()

# add to dictionaries
raw_dictionary['simulation overview smeared'] = raw_data
processed_dictionary['simulation overview smeared'] = processed_data

