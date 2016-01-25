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


### TOPAS-C preamp ############################################################


output_path = os.path.join(directory, 'TOPAS_C_preamp_healed.png')

force_workup = False

if not os.path.isfile(output_path) or force_workup:
    # load data from pickle
    data = wt.data.from_pickle('TOPAS_C_full_preamp.p')
    # clip based on intensity
    data.amplitude.clip(zmin=0.1, zmax=5)
    # clip based on width
    data.width.clip(zmin=10, zmax=100)
    # clip based on center
    data.center.clip(zmin=1140, zmax=1620)
    # share NaNs
    data.share_nans()
    # heal copy of data
    healed_data = data.copy()
    healed_data.heal('center')
    healed_data.heal('amplitude')
    datas = [healed_data, data]
    output_paths = [output_path, output_path.replace('healed', 'real')]
    for data, output_path in zip(datas, output_paths):
        # prepare plot
        fig = plt.figure(figsize=[14, 6])
        gs = grd.GridSpec(1, 5, hspace=0.05, wspace=0.05, width_ratios=[20, 1, 5, 20, 1])
        # intensity
        cmap = wt.artists.colormaps['default']
        cmap.set_under([0.75]*3, 1)
        ax0 = plt.subplot(gs[0])
        xi = data.c1.points
        yi = data.d1.points
        zi = data.amplitude.values
        X, Y, Z = wt.artists.pcolor_helper(xi, yi, zi)
        cax = ax0.pcolor(X, Y, Z, vmin=0, vmax=np.nanmax(Z), cmap=cmap)
        ax0.set_xlim(xi.min(), xi.max())
        ax0.set_ylim(1.35, 1.8)
        ax0.set_xlabel('C1 (degrees)', fontsize=16)
        ax0.set_ylabel('D1 (mm)', fontsize=16)
        ax0.grid()
        plt.colorbar(cax, plt.subplot(gs[1]))
        wt.artists.corner_text('intensity (a.u.)', ax=ax0, fontsize=16)
        ax0.contour(xi, yi, zi, 5, colors='k')
        # color
        cmap = wt.artists.colormaps['rainbow']
        cmap.set_under([0.75]*3, 1)
        ax1 = plt.subplot(gs[3])
        xi = data.c1.points
        yi = data.d1.points
        zi = data.center.values
        X, Y, Z = wt.artists.pcolor_helper(xi, yi, zi)
        cax = ax1.pcolor(X, Y, Z, vmin=np.nanmin(Z), vmax=np.nanmax(Z), cmap=cmap)
        ax1.set_xlim(xi.min(), xi.max())
        ax1.set_ylim(1.35, 1.8)
        ax1.set_xlabel('C1 (degrees)', fontsize=16)
        ax1.set_ylabel('D1 (mm)', fontsize=16)
        ax1.grid()
        plt.colorbar(cax, plt.subplot(gs[4]))
        wt.artists.corner_text('color (nm)', ax=ax1, fontsize=16)
        ax1.contour(xi, yi, zi, 25, colors='k')
        # finish
        plt.savefig(output_path, transparent=True, dpi=300)
        plt.close(fig)


### TOPAS-C poweramp ##########################################################


output_path = os.path.join(directory, 'TOPAS_C_poweramp_healed.png')

force_workup = True

if not os.path.isfile(output_path) or force_workup:
    # load data from pickle
    data = wt.data.from_pickle('TOPAS_C_full_poweramp.p')
    # clean up data
    data.amplitude.clip(zmin=0.1, zmax=5)
    
    # prepare plot
    fig = plt.figure(figsize=[14, 14])
    gs = grd.GridSpec(5, 5, hspace=0.05, wspace=0.05)
    
    # intensities -------------------------------------------------------------
    idxs = [idx for idx in np.ndindex(5, 5)]
    for i, idx in zip(range(data.w1.points.size), idxs):
        # get data
        w1 = data.w1.points[i]
        chop = data.chop('c2', 'd2', {'w1': [w1, 'nm']}, verbose=False)[0]
        xi = chop.c2.points
        yi = chop.d2.points
        zi = chop.amplitude.values
        # plot
        cmap = wt.artists.colormaps['default']
        ax = plt.subplot(gs[idx])
        X, Y, Z = wt.artists.pcolor_helper(xi, yi, zi)
        ax.pcolor(X, Y, Z, vmin=0, vmax=4, cmap=cmap)
        ax.set_xlim(xi.min(), xi.max())
        ax.set_ylim(yi.min(), yi.max())
        wt.artists.corner_text(str(w1), fontsize=16)
        
        
    # centers -----------------------------------------------------------------
        
    # finish
    plt.savefig(output_path, transparent=True, dpi=300)
    plt.close(fig)


### OPA800 signal and idler motortune #########################################


# TODO:


### OPA800 DFG Mixer Motortune ################################################


# TODO:
