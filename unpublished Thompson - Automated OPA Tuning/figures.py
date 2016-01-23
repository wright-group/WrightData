### import ####################################################################


import os
import itertools

import matplotlib.pyplot as plt
import matplotlib.gridspec as grd

import scipy
from scipy.optimize import leastsq
from scipy.interpolate import griddata, interp1d, interp2d, UnivariateSpline

import numpy as np

import WrightTools as wt


### define ####################################################################


directory = os.path.dirname(__file__)


### TOPAS-C amplitude and center (preamp) #####################################

output_path = 'TOPAS_C_amplitude_and_center.png'

force_workup = False

if not os.path.isfile(output_path) or force_workup:
    data = wt.data.from_PyCMDS('TOPAS_C_full_preamp')
    1/0
    
    
    
    
    
    
    
    
    
    c1 = npz['c1']
    d1 = npz['d1']
    outs = npz['outs']
    outs = outs.T
    cen = outs[0].flatten()
    wid = outs[1].flatten()
    amp = outs[2].flatten()
    # grid out c1, d1
    c1_grid, d1_grid = np.meshgrid(c1, d1, indexing='xy')
    c1_list = c1_grid.flatten()
    d1_list = d1_grid.flatten()
    # remove points with amplitudes that are ridiculous
    amp[amp<0.1] = np.nan
    amp[amp>4] = np.nan
    # remove points with centers that are ridiculous
    cen[cen<1150] = np.nan
    cen[cen>1650] = np.nan
    # remove points with widths that are ridiculous
    wid[wid<5] = np.nan
    wid[wid>500] = np.nan
    amp, cen, wid, c1_list, d1_list = wt.kit.remove_nans_1D([amp, cen, wid, c1_list, d1_list])
    # grid data
    xi = tuple(np.meshgrid(c1, d1, indexing='xy'))
    c1_grid, d1_grid = xi
    points = tuple([c1_list, d1_list])
    amp = griddata(points, amp, xi, method='cubic')
    amp /= np.nanmax(amp)
    cen = griddata(points, cen, xi, method='cubic')
    # prepare plot
    fig = plt.figure(figsize=[14, 6])
    gs = grd.GridSpec(1, 5, hspace=0.05, wspace=0.05, width_ratios=[20, 1, 5, 20, 1])
    # intensity
    cmap = wt.artists.colormaps['default']
    cmap.set_under([0.75]*3, 1)
    ax0 = plt.subplot(gs[0])
    X, Y, Z = wt.artists.pcolor_helper(c1, d1, amp)
    cax = ax0.pcolor(X, Y, Z, vmin=0, vmax=np.nanmax(Z), cmap=cmap)
    ax0.set_xlim(c1.min(), c1.max())
    ax0.set_ylim(1.35, 1.8)
    ax0.set_xlabel('C1 (degrees)', fontsize=16)
    ax0.set_ylabel('D1 (mm)', fontsize=16)
    ax0.grid()
    plt.colorbar(cax, plt.subplot(gs[1]))
    wt.artists.corner_text('intensity (a.u.)', ax=ax0, fontsize=16)
    ax0.contour(c1, d1, amp, 5, colors='k')
    # color
    cmap = wt.artists.colormaps['rainbow']
    cmap.set_under([0.75]*3, 1)
    ax1 = plt.subplot(gs[3])
    X, Y, Z = wt.artists.pcolor_helper(c1, d1, cen)
    cax = ax1.pcolor(X, Y, Z, vmin=np.nanmin(Z), vmax=np.nanmax(Z), cmap=cmap)
    ax1.set_xlim(c1.min(), c1.max())
    ax1.set_ylim(1.35, 1.8)
    ax1.set_xlabel('C1 (degrees)', fontsize=16)
    ax1.set_ylabel('D1 (mm)', fontsize=16)
    ax1.grid()
    plt.colorbar(cax, plt.subplot(gs[4]))
    wt.artists.corner_text('color (nm)', ax=ax1, fontsize=16)
    ax1.contour(c1, d1, cen, 25, colors='k')
    # finish
    plt.savefig(output_path, transparent=True, dpi=300)
    plt.close('all')





# not sure if I need this...
curves = [r'TOPAS-C\OPA1 (10743) base - 2015.11.10 09_25_25.crv',
          r'TOPAS-C\OPA1 (10743) mixer1 - 2015.11.09.crv',
          r'TOPAS-C\OPA1 (10743) mixer2 - 2015.10.26 12_39_55.crv',
          r'TOPAS-C\OPA1 (10743) mixer3 - 2013.06.01.crv']


### TOPAS-C amplitude and center (poweramp) ###################################


# TODO:


### OPA800 signal and idler motortune #########################################


# TODO:


### OPA800 DFG Mixer Motortune ################################################


# TODO:
