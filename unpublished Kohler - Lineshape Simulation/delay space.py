### import ####################################################################


import os

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import matplotlib.colors as mplcolors
import matplotlib.patheffects as PathEffects
from mpl_toolkits.axes_grid1 import make_axes_locatable
matplotlib.rcParams['font.size'] = 18
plt.ioff()
plt.close('all')

import WrightTools as wt


### define ####################################################################


data_folder = 'npzs'
template_fname = r'npzs\dpr {0} TO {1} w1 w2 d1 d2 arr.npz'
dprs = [0.5, 1., 2.0]


### delay space figure ########################################################


if True:
    # prepare figure
    plt.figure(figsize=(30, 10))
    gs = grd.GridSpec(1, 3)
    # subplot method
    def plot(sps, dpr):
        ax = plt.subplot(sps)        
        cmap = wt.artists.colormaps['default']
        # get data from zip
        filepath = template_fname.format(dpr, 'all')
        npz = np.load(filepath)
        d1 = -npz['d1']/50.
        d2 = -npz['d2']/50.
        arr = np.sqrt(npz['arr'][20, 20].T)
        TO_all = arr.copy()
        arr /= arr.max()
        # plot amplitude
        levels = np.linspace(0, 1, 200)
        ax.contourf(d1, d2, arr, levels=levels, cmap=cmap)
        ax.set_xlim(d1.min(), d1.max())
        ax.set_ylim(d2.min(), d2.max())    
        ax.grid()
        ax.axhline(0, c='grey', ls='-', lw=3)
        ax.axvline(0, c='grey', ls='-', lw=3)
        wt.artists.diagonal_line(d1, d2, ax=ax, c='grey', ls='-', lw=3)
        def plot_contours(TO):
            filepath = template_fname.format(dpr, TO)
            npz = np.load(filepath)
            d1 = -npz['d1']/50.
            d2 = -npz['d2']/50.
            arr = np.sqrt(npz['arr'][20, 20].T)
            arr /= TO_all
            arr /= arr.max()
            # generate levels
            levels = [0.5, 0.75]
            current = 0.75
            for _ in range(10):
                current = current+((1-current)/1.5)
                levels.append(current)
            CS = ax.contour(d1, d2, arr, colors='k', levels=levels)
            plt.clabel(CS, inline=1, fontsize=12)
        plot_contours(1)
        plot_contours(3)
        plot_contours(5)
        plot_contours(6)        
    plot(gs[0, 0], dprs[0])
    plot(gs[0, 1], dprs[1])
    plot(gs[0, 2], dprs[2])
    
    
    
    
    plt.savefig('delay space.png', dpi=300, transparency=True)
    plt.close('all')