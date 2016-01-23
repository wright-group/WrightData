### import ####################################################################


import os
import itertools

import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
plt.close('all')

import scipy
from scipy.optimize import leastsq
from scipy.interpolate import griddata, interp1d, interp2d, UnivariateSpline

import numpy as np

import WrightTools as wt


### define ####################################################################


directory = os.path.dirname(__file__)

props = dict(boxstyle='square', facecolor='white', alpha=0.8)


### WMELs #####################################################################

output_path = os.path.join(directory, 'WMELs.png')

force_plotting = True

if not os.path.isfile(output_path) or force_plotting:
    # prepare for plot
    fig = plt.figure(figsize=[10, 8])
    gs = grd.GridSpec(3, 8, width_ratios=[1.5, 2, 1, 1, 1, 1, 1, 1])
    # large levels figure
    ax = plt.subplot(gs[:, 0])
    # plot energies
    energies = np.linspace(0, 1, 4)
    for energy in energies:
        ax.axhline(energy, color='k', linewidth=2, ls='-')
    # add arrow
    def add_arrow(between, x_pos, width):
        head_size = 0.05
        color = 'k'
        # calculate arrow length
        arrow_length = energies[between[1]] - energies[between[0]]
        arrow_end = energies[between[1]]
        width_compensation_factor = 1e-3
        if arrow_length > 0:
            direction = 1
            y_poss = [energies[between[0]] + width*width_compensation_factor, energies[between[1]] - head_size]
        elif arrow_length < 0:
            direction = -1
            y_poss = [energies[between[0]] - width*width_compensation_factor, energies[between[1]] + head_size]
        # add line
        length = abs(y_poss[0] - y_poss[1])
        line = ax.plot([x_pos, x_pos], y_poss, linestyle='-', color=color, linewidth=width)
        # add arrow head
        arrow_head = ax.arrow(x_pos, arrow_end - head_size * direction, 
                              0, 0.0001*direction,
                              head_width=head_size*width, 
                              head_length=head_size,
                              fc=color, ec=color, linestyle='solid', linewidth=0)        
    add_arrow([0, 1], 0.25, 14)
    #ax.text(-0.1, 0.15, '0', fontsize=12, horizontalalignment='center', verticalalignment='center')
    add_arrow([1, 2], 0.25, 12)
    #ax.text(-0.1, 0.5, '1', fontsize=12, horizontalalignment='center', verticalalignment='center')
    add_arrow([2, 3], 0.25, 10)
    #ax.text(-0.1, 0.85, '2', fontsize=12, horizontalalignment='center', verticalalignment='center')
    add_arrow([3, 2], 0.75, 8)
    #ax.text(1.1, 0.85, '3', fontsize=12, horizontalalignment='center', verticalalignment='center')
    add_arrow([2, 1], 0.75, 4)
    #ax.text(1.1, 0.5, '4', fontsize=12, horizontalalignment='center', verticalalignment='center')
    add_arrow([1, 0], 0.75, 2)
    #ax.text(1.1, 0.15, '5', fontsize=12, horizontalalignment='center', verticalalignment='center')
    state_text_buffer = 0.25
    state_names = ['N=0', '1', '2', '3']
    for i in range(len(energies)):
        ax.text(-state_text_buffer, energies[i], state_names[i], fontsize=12, verticalalignment='center', horizontalalignment ='right')
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.025, 1.025)
    ax.axis('off')
    # pathway 1 ---------------------------------------------------------------
    energies = [0., 0.5, 1.]
    state_text_buffer = 0.25
    state_names = ['g', 'a', '2a']
    # pathway 1 alpha
    ax = plt.subplot(gs[0, 2])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='I')
    wmel.add_arrow(0, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(1, [0, 1], 'bra', '-2', color='b')
    wmel.add_arrow(2, [1, 0], 'bra', '2\'', color='b')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    for i in range(len(energies)):
        ax.text(-state_text_buffer, energies[i], state_names[i], fontsize=12, verticalalignment='center', horizontalalignment ='right')
    # pathway 1 beta
    ax = plt.subplot(gs[1, 2])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(1, [0, 1], 'bra', '-2', color='b')
    wmel.add_arrow(2, [1, 2], 'bra', '2\'', color='b')
    wmel.add_arrow(3, [2, 1], 'out', '', color='r')
    for i in range(len(energies)):
        ax.text(-state_text_buffer, energies[i], state_names[i], fontsize=12, verticalalignment='center', horizontalalignment ='right')
    # pathway 1 gamma
    ax = plt.subplot(gs[2, 2])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(1, [1, 0], 'ket', '-2', color='b')
    wmel.add_arrow(2, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    for i in range(len(energies)):
        ax.text(-state_text_buffer, energies[i], state_names[i], fontsize=12, verticalalignment='center', horizontalalignment ='right')
    # pathway 2 ---------------------------------------------------------------
    # pathway 2 alpha
    ax = plt.subplot(gs[0, 3])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='II')
    wmel.add_arrow(0, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(1, [1, 2], 'ket', '2\'', color='b')
    wmel.add_arrow(2, [2, 1], 'ket', '-2', color='b')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    # pathway 2 beta
    ax = plt.subplot(gs[1, 3])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(1, [1, 2], 'ket', '2\'', color='b')
    wmel.add_arrow(2, [0, 1], 'bra', '-2', color='b')
    wmel.add_arrow(3, [2, 1], 'out', '', color='r')
    # pathway 3 ---------------------------------------------------------------
    # pathway 3 alpha
    ax = plt.subplot(gs[0, 4])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='III')
    wmel.add_arrow(0, [0, 1], 'bra', '-2', color='b')
    wmel.add_arrow(1, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(2, [1, 0], 'bra', '2\'', color='b')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    # pathway 3 beta
    ax = plt.subplot(gs[1, 4])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'bra', '-2', color='b')
    wmel.add_arrow(1, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(2, [1, 2], 'ket', '2\'', color='b')
    wmel.add_arrow(3, [2, 1], 'out', '', color='r')
    # pathway 3 gamma
    ax = plt.subplot(gs[2, 4])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'bra', '-2', color='b')
    wmel.add_arrow(1, [1, 0], 'bra', '1', color='r')
    wmel.add_arrow(2, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    # pathway 4 ---------------------------------------------------------------
    # pathway 4 alpha
    ax = plt.subplot(gs[0, 5])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='IV')
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(1, [1, 2], 'ket', '1', color='r')
    wmel.add_arrow(2, [2, 1], 'ket', '-2', color='b')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    # pathway 4 beta
    ax = plt.subplot(gs[1, 5])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(1, [1, 2], 'ket', '1', color='r')
    wmel.add_arrow(2, [0, 1], 'bra', '-2', color='b')
    wmel.add_arrow(3, [2, 1], 'out', '', color='r')
    # pathway 5 ---------------------------------------------------------------
    # pathway 5 alpha
    ax = plt.subplot(gs[0, 6])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='V')
    wmel.add_arrow(0, [0, 1], 'bra', '-2', color='b')
    wmel.add_arrow(1, [1, 0], 'bra', '2\'', color='b')
    wmel.add_arrow(2, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    # pathway 5 beta
    ax = plt.subplot(gs[1, 6])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'bra', '-2', color='b')
    wmel.add_arrow(1, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(2, [1, 2], 'ket', '1', color='r')
    wmel.add_arrow(3, [2, 1], 'out', '', color='r')
    # pathway 5 gamma
    ax = plt.subplot(gs[2, 6])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'bra', '-2', color='b')
    wmel.add_arrow(1, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(2, [1, 0], 'bra', '1', color='r')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    # pathway 6 ---------------------------------------------------------------
    text_buffer = 1.3
    # pathway 6 alpha
    ax = plt.subplot(gs[0, 7])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='VI')
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(1, [1, 0], 'ket', '-2', color='b')
    wmel.add_arrow(2, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    ax.text(text_buffer, 0.5, r'$\mathrm{\alpha}$', fontsize=25, verticalalignment='center', horizontalalignment='left')
    # pathway 6 beta
    ax = plt.subplot(gs[1, 7])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(1, [0, 1], 'bra', '-2', color='b')
    wmel.add_arrow(2, [1, 2], 'ket', '1', color='r')
    wmel.add_arrow(3, [2, 1], 'out', '', color='r')
    ax.text(text_buffer, 0.5, r'$\mathrm{\beta}$', fontsize=25, verticalalignment='center', horizontalalignment='left')
    # pathway 6 gamma
    ax = plt.subplot(gs[2, 7])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(1, [0, 1], 'bra', '-2', color='b')
    wmel.add_arrow(2, [1, 0], 'bra', '1', color='r')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    ax.text(text_buffer, 0.5, r'$\mathrm{\gamma}$', fontsize=25, verticalalignment='center', horizontalalignment='left')
    # labels
    dist = 0.05
    fig.text(0.075, 1-dist, 'A', bbox=props, horizontalalignment='left', verticalalignment='top', fontsize=16)
    fig.text(0.35, 1-dist, 'B', bbox=props, horizontalalignment='left', verticalalignment='top', fontsize=16)    
    # save
    plt.savefig(output_path, dpi=300, transparent=True)
    