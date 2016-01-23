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
