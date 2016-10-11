### import ####################################################################


import os
import itertools

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.patheffects as PathEffects
from matplotlib.patches import Ellipse
plt.close('all')
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage[cm]{sfmath}'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.serif'] = 'cm'
matplotlib.rcParams['font.sans-serif'] = 'cm'
matplotlib.rcParams['font.size'] = 14

import scipy
from scipy import constants
from scipy import ndimage
from scipy.optimize import leastsq
from scipy.interpolate import griddata, interp1d, interp2d, UnivariateSpline

import numpy as np

import WrightTools as wt


### define ####################################################################


directory = os.path.dirname(__file__)

# colormaps
cmap = wt.artists.colormaps['default']
coolwarm = plt.get_cmap('coolwarm_r')

props = dict(boxstyle='square', facecolor='white', alpha=0.8)

def adjust_spines(ax, color, lw=3):
    for key in ['bottom', 'top', 'right', 'left']:
        ax.spines[key].set_color(color)
        ax.spines[key].set_linewidth(lw)
        ax.spines[key].zorder = 1000

def zoom_arrs(*args):
    return [ndimage.interpolation.zoom(arr, 5) for arr in args]

delta_omega = 2 * np.log(2) / (50*np.pi * 3e10) * 1e15  # wn
def normalize_frequency(arr):
    arr -= 7000.
    arr /= delta_omega
    return arr

delta_t = 50.  # fs
def normalize_delay(arr):
    arr /= delta_t
    return arr

def normalized_gauss(t, FWHM):
    '''
    amplitude level
    '''
    sigma = FWHM / (2.*np.sqrt(np.log(2.)))
    out = np.exp(-0.5*(t/sigma)**2)
    out /= sigma * np.sqrt(2*np.pi)
    return out        

def delay_label(kind):
    if kind == 1:
        return r'$\mathsf{\tau_{22^{\prime}}/\Delta_t}$'
    elif kind == 2:
        return r'$\mathsf{\tau_{21}/\Delta_t}$'
    else:
        return None

def time_label():
    return r'$\mathsf{t/\Delta_t}$'

def frequency_label(sub=''):
    return r'$\mathsf{(\omega_{' + str(sub) + r'}-\omega_{10})/\Delta_{\omega}}$'

freq_ticks = [-1, 0, 1]
delay_ticks = [-2, 0, 2]

contour_levels = np.linspace(0, 1, 11)[1:-1]


### WMELs #####################################################################


output_path = os.path.join(directory, 'WMELs.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # prepare for plot
    fig = plt.figure(figsize=[6.5, 6])
    gs = grd.GridSpec(3, 8, width_ratios=[1.5, 1, 1, 1, 1, 1, 1, 1])
    # large levels figure
    ax = plt.subplot(gs[:, 0])
    # plot energies
    energies = [0., 0.175, 0.35, 0.525, 0.825, 1.]
    for energy in energies:
        ax.axhline(energy, color='k', linewidth=2, ls='-')
    # add arrow method for subplot a
    def add_arrow(between, x_pos, width):
        head_size = 0.05
        color = 'k'
        # calculate arrow length
        arrow_length = energies[between[1]] - energies[between[0]]
        arrow_end = energies[between[1]]
        width_compensation_factor = 2e-3
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
    add_arrow([1, 2], 0.25, 12)
    add_arrow([2, 3], 0.25, 10)
    add_arrow([4, 5], 0.25, 2)
    add_arrow([3, 2], 0.75, 8)
    add_arrow([2, 1], 0.75, 4)
    add_arrow([1, 0], 0.75, 2)
    add_arrow([5, 4], 0.75, 14)
    state_text_buffer = 0.25
    state_names = ['n=0', 'n=1', 'n=2', 'n=3', 'n=N-1', 'n=N']
    for i in range(len(energies)):
        ax.text(-state_text_buffer, energies[i], state_names[i], fontsize=18, verticalalignment='center', horizontalalignment ='right')
    ax.text(-state_text_buffer, 0.675, '...', rotation=-90, fontsize=18, verticalalignment='center', horizontalalignment ='right')
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.025, 1.025)
    ax.axis('off')
    # pathway 1 ---------------------------------------------------------------
    energies = [0., 0.5, 1.]
    state_text_buffer = 0.25
    state_names = ['0', '1', '2']
    # pathway 1 alpha
    ax = plt.subplot(gs[0, 2])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='I')
    wmel.add_arrow(0, [0, 1], 'ket', '1', color='orange')
    wmel.add_arrow(1, [0, 1], 'bra', '2', color='m')
    wmel.add_arrow(2, [1, 0], 'bra', '2\'', color='m')
    wmel.add_arrow(3, [1, 0], 'out', '', color='grey')
    for i in range(len(energies)):
        ax.text(-state_text_buffer, energies[i], state_names[i], fontsize=12, verticalalignment='center', horizontalalignment ='right')
    # pathway 1 beta
    ax = plt.subplot(gs[1, 2])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '1', color='orange')
    wmel.add_arrow(1, [0, 1], 'bra', '2', color='m')
    wmel.add_arrow(2, [1, 2], 'bra', '2\'', color='m')
    wmel.add_arrow(3, [2, 1], 'out', '', color='grey')
    for i in range(len(energies)):
        ax.text(-state_text_buffer, energies[i], state_names[i], fontsize=12, verticalalignment='center', horizontalalignment ='right')
    # pathway 1 gamma
    ax = plt.subplot(gs[2, 2])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '1', color='orange')
    wmel.add_arrow(1, [1, 0], 'ket', '2', color='m')
    wmel.add_arrow(2, [0, 1], 'ket', '2\'', color='m')
    wmel.add_arrow(3, [1, 0], 'out', '', color='grey')
    for i in range(len(energies)):
        ax.text(-state_text_buffer, energies[i], state_names[i], fontsize=12, verticalalignment='center', horizontalalignment ='right')
    # pathway 2 ---------------------------------------------------------------
    # pathway 2 alpha
    ax = plt.subplot(gs[0, 3])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='II')
    wmel.add_arrow(0, [0, 1], 'ket', '1', color='orange')
    wmel.add_arrow(1, [1, 2], 'ket', '2\'', color='m')
    wmel.add_arrow(2, [2, 1], 'ket', '2', color='m')
    wmel.add_arrow(3, [1, 0], 'out', '', color='grey')
    # pathway 2 beta
    ax = plt.subplot(gs[1, 3])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '1', color='orange')
    wmel.add_arrow(1, [1, 2], 'ket', '2\'', color='m')
    wmel.add_arrow(2, [0, 1], 'bra', '2', color='m')
    wmel.add_arrow(3, [2, 1], 'out', '', color='grey')
    # pathway 3 ---------------------------------------------------------------
    # pathway 3 alpha
    ax = plt.subplot(gs[0, 4])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='III')
    wmel.add_arrow(0, [0, 1], 'bra', '2', color='m')
    wmel.add_arrow(1, [0, 1], 'ket', '1', color='orange')
    wmel.add_arrow(2, [1, 0], 'bra', '2\'', color='m')
    wmel.add_arrow(3, [1, 0], 'out', '', color='grey')
    # pathway 3 beta
    ax = plt.subplot(gs[1, 4])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'bra', '2', color='m')
    wmel.add_arrow(1, [0, 1], 'ket', '1', color='orange')
    wmel.add_arrow(2, [1, 2], 'ket', '2\'', color='m')
    wmel.add_arrow(3, [2, 1], 'out', '', color='grey')
    # pathway 3 gamma
    ax = plt.subplot(gs[2, 4])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'bra', '2', color='m')
    wmel.add_arrow(1, [1, 0], 'bra', '1', color='orange')
    wmel.add_arrow(2, [0, 1], 'ket', '2\'', color='m')
    wmel.add_arrow(3, [1, 0], 'out', '', color='grey')
    # pathway 4 ---------------------------------------------------------------
    # pathway 4 alpha
    ax = plt.subplot(gs[0, 5])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='IV')
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='m')
    wmel.add_arrow(1, [1, 2], 'ket', '1', color='orange')
    wmel.add_arrow(2, [2, 1], 'ket', '2', color='m')
    wmel.add_arrow(3, [1, 0], 'out', '', color='grey')
    # pathway 4 beta
    ax = plt.subplot(gs[1, 5])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='m')
    wmel.add_arrow(1, [1, 2], 'ket', '1', color='orange')
    wmel.add_arrow(2, [0, 1], 'bra', '2', color='m')
    wmel.add_arrow(3, [2, 1], 'out', '', color='grey')
    # pathway 5 ---------------------------------------------------------------
    # pathway 5 alpha
    ax = plt.subplot(gs[0, 6])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='V')
    wmel.add_arrow(0, [0, 1], 'bra', '2', color='m')
    wmel.add_arrow(1, [1, 0], 'bra', '2\'', color='m')
    wmel.add_arrow(2, [0, 1], 'ket', '1', color='orange')
    wmel.add_arrow(3, [1, 0], 'out', '', color='grey')
    # pathway 5 beta
    ax = plt.subplot(gs[1, 6])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'bra', '2', color='m')
    wmel.add_arrow(1, [0, 1], 'ket', '2\'', color='m')
    wmel.add_arrow(2, [1, 2], 'ket', '1', color='orange')
    wmel.add_arrow(3, [2, 1], 'out', '', color='grey')
    # pathway 5 gamma
    ax = plt.subplot(gs[2, 6])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'bra', '2', color='m')
    wmel.add_arrow(1, [0, 1], 'ket', '2\'', color='m')
    wmel.add_arrow(2, [1, 0], 'bra', '1', color='orange')
    wmel.add_arrow(3, [1, 0], 'out', '', color='grey')
    # pathway 6 ---------------------------------------------------------------
    text_buffer = 1.3
    # pathway 6 alpha
    ax = plt.subplot(gs[0, 7])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='VI')
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='m')
    wmel.add_arrow(1, [1, 0], 'ket', '2', color='m')
    wmel.add_arrow(2, [0, 1], 'ket', '1', color='orange')
    wmel.add_arrow(3, [1, 0], 'out', '', color='grey')
    ax.text(text_buffer, 0.5, r'$\mathsf{\alpha}$', fontsize=25, verticalalignment='center', horizontalalignment='left')
    # pathway 6 beta
    ax = plt.subplot(gs[1, 7])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='m')
    wmel.add_arrow(1, [0, 1], 'bra', '2', color='m')
    wmel.add_arrow(2, [1, 2], 'ket', '1', color='orange')
    wmel.add_arrow(3, [2, 1], 'out', '', color='grey')
    ax.text(text_buffer, 0.5, r'$\mathsf{\beta}$', fontsize=25, verticalalignment='center', horizontalalignment='left')
    # pathway 6 gamma
    ax = plt.subplot(gs[2, 7])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='m')
    wmel.add_arrow(1, [0, 1], 'bra', '2', color='m')
    wmel.add_arrow(2, [1, 0], 'bra', '1', color='orange')
    wmel.add_arrow(3, [1, 0], 'out', '', color='grey')
    ax.text(text_buffer, 0.5, r'$\mathsf{\gamma}$', fontsize=25, verticalalignment='center', horizontalalignment='left')
    # labels
    dist = 0.05
    fig.text(0.075, 1-dist, 'a', bbox=props, horizontalalignment='left', verticalalignment='top', fontsize=18, alpha=1)
    fig.text(0.33, 1-dist, 'b', bbox=props, horizontalalignment='left', verticalalignment='top', fontsize=18, alpha=1)    
    # finish
    wt.artists.subplots_adjust(fig)
    plt.savefig(output_path, dpi=300, transparent=True,  pad_inches=1.)
    plt.close(fig)


### simulation overview #######################################################


output_path = os.path.join(directory, 'simulation overview.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # prepare for plot
    fig, gs = wt.artists.create_figure(width='double', nrows=2, cols=[1, 0, 0, 0, 1.5, 0, 0, 0, 1, 'cbar'],
                                       hspace=0.75)
    # set of SQCs vs lab time -------------------------------------------------
    ax = plt.subplot(gs[0, 0])                         
    p = os.path.join(directory, 'precalculated', 'simulation overview a.hdf5')
    d = wt.kit.read_h5(p)
    xi = d['xi']
    # pulse
    yi = d['pulse']
    ax.fill_between(xi, 0, yi, facecolor='grey', alpha=0.25)
    # density matrix terms
    cs = ['b', 'g', 'r']
    labels = ['$\mathsf{2.0}$', '$\mathsf{1.0}$', '$\mathsf{0.5}$']
    for i, c in enumerate(cs):
        yi = d['rho%d'%i]
        ax.plot(xi, yi, c=c, lw=2, label=labels[i])
    # finish
    ax.legend(fontsize=18, loc='upper right')
    ax.set_ylim(0, 1.15)
    ax.set_xlim(-3, 11)
    ax.set_xlabel(time_label(), fontsize=18)
    ax.grid()
    wt.artists.corner_text('a', ax=ax, fontsize=18, background_alpha=1)
    ax.set_ylabel('amplitude', fontsize=18)
    adjust_spines(ax, 'grey')
    # evolution of density matrix terms in pw5 --------------------------------
    ax = plt.subplot(gs[0, 1:5])
    p = os.path.join(directory, 'precalculated', 'simulation overview b.hdf5')
    d = wt.kit.read_h5(p)
    xi = d['xi']
    # -2
    yi = d['-2']
    ax.fill_between(xi, 0, yi, facecolor='b', alpha=0.25)
    ax.text(-4, 0.1, '2', fontsize=18, ha='center', va='center')
    ax.text(-4, 0.3, 'x', fontsize=18, ha='center', va='center')
    # 2'
    yi = d['2prime']
    ax.fill_between(xi, 0, yi, facecolor='b', alpha=0.25)
    ax.text(-2, 0.1, '2\'', fontsize=18, ha='center', va='center')
    ax.text(-2, 0.3, 'y', fontsize=18, ha='center', va='center')
    # 1
    yi = d['1']
    ax.fill_between(xi, 0, yi, facecolor='r', alpha=0.25)
    ax.text(0, 0.1, '1', fontsize=18, ha='center', va='center')
    ax.text(0   , 0.3, 'z', fontsize=18, ha='center', va='center')
    # plot density matrix terms
    # ga
    yi = d['ga']
    ax.plot(xi, yi, lw=2, c='g', label=r'$\mathsf{\rho_{01}}$')
    # aa
    yi = d['aa']
    ax.plot(xi, yi, lw=2, c='m', label=r'$\mathsf{\rho_{11}}$')
    # ag2
    yi = d['ag']
    ax.plot(xi, yi, lw=2, c='c', label=r'$\mathsf{\rho_{01}}$')
    # arrows
    ax.annotate('', xy=(-2, 1.05), xytext=(-4, 1.05), arrowprops=dict(arrowstyle="->", ec='k'))
    ax.annotate('', xy=(0, 1.1), xytext=(-4, 1.1), arrowprops=dict(arrowstyle="->"))
    ax.text(-4.1, 1.05, r'$\mathsf{\tau_{22^{\prime}}}$', ha='right', va='center', fontsize=18)
    ax.text(0.1, 1.1, r'$\mathsf{\tau_{21}}$', ha='left', va='center', fontsize=18)
    # finish
    ax.legend(fontsize=18, loc='right')
    ax.set_ylim(0, 1.15)
    ax.set_xlim(-7, 7)
    plt.grid()    
    wt.artists.corner_text('b', ax=ax, fontsize=18, background_alpha=1)
    ax.set_xlabel(time_label(), fontsize=18)
    plt.setp(ax.get_yticklabels(), visible=False)
    adjust_spines(ax, 'g')
    # FT of FID above and mono pass function ----------------------------------
    ax = plt.subplot(gs[0, 5:])
    p = os.path.join(directory, 'precalculated', 'simulation overview c.hdf5')
    d = wt.kit.read_h5(p)
    xi = d['xi']
    yi = d['yi']
    ax.plot(xi, yi, c='c', lw=2)
    # mono window
    window = 150.  # wn
    limit = (window/2) / (2 * np.log(2) / (50*np.pi * 3e10) * 1e15)
    ax.axvline(-limit, c='k', ls='--')
    ax.axvline(limit, c='k', ls='--')
    low_index = np.argmin(np.abs(xi+limit))
    high_index = np.argmin(np.abs(xi-limit)) + 1
    ax.fill_between(xi[low_index:high_index], 0, yi[low_index:high_index], facecolor='c', alpha=0.25, edgecolor='none')
    # finish
    ax.set_xlim(-2, 2)
    ax.set_ylim(0, 1.15)
    ax.set_xticks(freq_ticks)
    ax.grid()
    wt.artists.corner_text('c', ax=ax, fontsize=18, background_alpha=1)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.set_xlabel(r'$\mathsf{(\omega_{out}-\omega_{1})/\Delta_{\omega}}$', fontsize=18)
    adjust_spines(ax, 'g')
    # 2D delay ----------------------------------------------------------------
    ax = plt.subplot(gs[1, 0])
    p = os.path.join(directory, 'measured', 'smear 0.0', 'dpr 1.0 TO all w1 w2 d1 d2.hdf5')
    d = wt.kit.read_h5(p)
    xi = normalize_delay(-d['d1'])
    yi = normalize_delay(-d['d2'])
    zi = np.sqrt(d['arr'][20, 20, :, :]).T
    zi /= zi.max()
    levels= np.linspace(0, 1, 200)
    xi, yi, zi = zoom_arrs(xi, yi, zi)
    ax.contourf(xi, yi, zi, cmap=cmap, levels=levels)
    ax.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
    ax.grid()
    ax.axhline(0, c='k', lw=2)
    ax.axvline(0, c='k', lw=2)
    ticks = [-2, 0, 2]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xlabel(delay_label(1), fontsize=18)
    ax.set_ylabel(delay_label(2), fontsize=18)
    wt.artists.diagonal_line(xi, yi, ls='-', lw=2)
    adjust_spines(ax, 'g')
    wt.artists.corner_text('d', ax=ax, fontsize=18, background_alpha=1)
    factor = 4
    ax.text(-0.5*factor, 0.5*factor, 'I', fontsize = 35, verticalalignment = 'center', horizontalalignment = 'center', path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")])
    ax.text(0.25*factor, 0.6*factor, 'II', fontsize = 35, verticalalignment = 'center', horizontalalignment = 'center', path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")])
    ax.text(-0.6*factor, -0.25*factor, 'III', fontsize = 35, verticalalignment = 'center', horizontalalignment = 'center', path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")])
    ax.text(0.6*factor, 0.25*factor, 'IV', fontsize = 35, verticalalignment = 'center', horizontalalignment = 'center', path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")])
    ax.text(-0.25*factor, -0.6*factor, 'V', fontsize = 35, verticalalignment = 'center', horizontalalignment = 'center', path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")])
    ax.text(0.5*factor, -0.5*factor, 'VI', fontsize = 35, verticalalignment = 'center', horizontalalignment = 'center', path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")])
    # wigner ------------------------------------------------------------------
    ax = plt.subplot(gs[1, 4])
    p = os.path.join(directory, 'measured', 'smear 0.0', 'dpr 1.0 TO all w1 w2 d1 d2.hdf5')
    d = wt.kit.read_h5(p)
    xi = normalize_frequency(d['w1'])
    yi = normalize_delay(-d['d2'])
    zi = np.sqrt(d['arr'][:, 20, 10, :]).T
    zi /= zi.max()
    levels= np.linspace(0, 1, 200)
    xi, yi, zi = zoom_arrs(xi, yi, zi)
    ax.contourf(xi, yi, zi, cmap=cmap, levels=levels)
    ax.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
    ax.grid()
    ax.axhline(0, c='k', lw=2)
    ax.axvline(0, lw=4, alpha=0.25, c='k')
    ax.set_xlim([-2, 2])
    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([-2, 0, 2])
    ax.set_xlabel(frequency_label('1'), fontsize=18)
    ax.set_ylabel(delay_label(2), fontsize=18)
    adjust_spines(ax, 'g')
    wt.artists.corner_text('e', ax=ax, fontsize=18, background_alpha=1)
    # 2D frequency ------------------------------------------------------------
    ax = plt.subplot(gs[1, 8])
    p = os.path.join(directory, 'measured', 'smear 0.0', 'dpr 1.0 TO all w1 w2 d1 d2.hdf5')
    d = wt.kit.read_h5(p)    
    xi = normalize_frequency(d['w1'])
    yi = normalize_frequency(d['w2'])
    zi = np.sqrt(d['arr'][:, :, 10, 10]).T
    zi /= zi.max()
    levels= np.linspace(0, 1, 200)
    xi, yi, zi = zoom_arrs(xi, yi, zi)
    ax.contourf(xi, yi, zi, cmap=cmap, levels=levels)
    ax.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
    ax.grid()
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_xticks(freq_ticks)
    ax.set_yticks(freq_ticks)
    wt.artists.diagonal_line(xi, yi)
    wt.artists.corner_text('f', ax=ax, fontsize=18, background_alpha=1)
    ax.set_xlabel(frequency_label('1'), fontsize=18)
    ax.set_ylabel(frequency_label('2'), fontsize=18)
    adjust_spines(ax, 'g')
    # finish ------------------------------------------------------------------
    # colorbar    
    cax = plt.subplot(gs[1, -1])
    cbar = wt.artists.plot_colorbar(cax=cax, label='amplitude')
    # finish
    wt.artists.savefig(output_path)
        

### fid vs dpr ################################################################


output_path = os.path.join(directory, 'fid vs dpr.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # prepare figure
    fig, gs = wt.artists.create_figure(nrows=5, cols=[1, 'cbar'],
                                       default_aspect=0.3, hspace=0.1)
    def add_data(ax, index, label, show_x=False):
       p = os.path.join(directory, 'precalculated', 'fid vs dpr.hdf5')
       d = wt.kit.read_h5(p)
       xi = d['xi']
       # color
       yi = d['rho%d'%index]
       freq = normalize_frequency(d['freq%d'%index])
       cNorm  = colors.Normalize(vmin=-0.1, vmax=0.1)
       scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=coolwarm)
       for i in range(len(xi[1:])):
           # calculate the color to use, given the instantaneous frequency
           colorVal = scalarMap.to_rgba(freq[i])
           ax.plot([xi[i+1], xi[i+1]], [0, yi[i]], color=colorVal)
       # envelope
       yi = d['rho%d'%index]
       plt.plot(xi, yi, color='k', lw=2)
       # pulse
       yi = d['pulse']
       plt.plot(xi, yi, color='grey', lw=2, zorder=15)
       # finish
       if show_x:
           ax.set_xlabel(r'$\mathsf{t / \Delta_t}$', fontsize=18)
       else:
           plt.setp(ax.get_xticklabels(), visible=False)
       ax.set_xlim(-1.5, 2.25)
       ax.set_ylim(0, 1.1)
       ax.axvline(0, c='k', lw=2)
       plt.setp(ax.get_yticklabels(), visible=False)
       wt.artists.corner_text(label, fontsize=18, background_alpha=1.)
       ax.grid()
    # fill out axes
    # 0
    ax = plt.subplot(gs[0, 0])
    add_data(ax, 0, label='$\infty$')
    adjust_spines(ax, 'grey')
    # 1
    ax = plt.subplot(gs[1, 0])
    add_data(ax, 1, label='2.0')
    adjust_spines(ax, 'b')
    # 2
    ax = plt.subplot(gs[2, 0])
    add_data(ax, 2, label='1.0')
    ax.set_ylabel('amplitude', fontsize=18, labelpad=20)
    adjust_spines(ax, 'g')
    # 3
    ax = plt.subplot(gs[3, 0])
    add_data(ax, 3, label='0.5')
    adjust_spines(ax, 'r')
    # 4
    ax = plt.subplot(gs[4, 0])
    add_data(ax, 4, label='0.0', show_x=True)
    adjust_spines(ax, 'grey')
    # colorbar    
    cax = plt.subplot(gs[:, 1])
    label = r'$\mathsf{((d\phi / dt)-\omega_{10})/\Delta_{\omega}}$'
    cbar = wt.artists.plot_colorbar(cax=cax, cmap=coolwarm, 
                                    ticks=[0., 0.05, 0.1], vlim=[0, 0.1],
                                    clim=[-0.1, 0.1], label=label)
    # finish
    plt.savefig(output_path, dpi=300, transparent=True, pad_inches=1)
    plt.close('fig')


### fid vs detuning ###########################################################


output_path = os.path.join(directory, 'fid vs detuning.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # prepare figure
    fig, gs = wt.artists.create_figure(nrows=5, cols=[1, 'cbar'],
                                       default_aspect=0.3, hspace=0.1)
    def add_data(ax, index, label, show_x=False):
       p = os.path.join(directory, 'precalculated', 'fid vs detuning.hdf5')
       d = wt.kit.read_h5(p)
       xi = d['xi']
       # color
       yi = d['rho%d'%index]
       freq = normalize_frequency(d['freq%d'%index])
       cNorm  = colors.Normalize(vmin=-2, vmax=2)
       scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=coolwarm)
       for i in range(len(xi[1:])):
           # calculate the color to use, given the instantaneous frequency
           colorVal = scalarMap.to_rgba(freq[i])
           ax.plot([xi[i+1], xi[i+1]], [0, yi[i]], color=colorVal)
       # envelope
       yi = d['rho%d'%index]
       plt.plot(xi, yi, color='k', lw=2)
       # pulse
       yi = d['pulse']
       plt.plot(xi, yi, color='grey', lw=2, zorder=15)
       # finish
       if show_x:
           ax.set_xlabel(r'$\mathsf{t / \Delta_t}$', fontsize=18)
       else:
           plt.setp(ax.get_xticklabels(), visible=False)
       ax.set_xlim(-1.5, 2.25)
       ax.set_ylim(0, 1.1)
       ax.axvline(0, c='k', lw=2)
       plt.setp(ax.get_yticklabels(), visible=False)
       wt.artists.corner_text(label, fontsize=18, background_alpha=1.)
       ax.grid()
       adjust_spines(ax, 'g')
    # fill out axes
    # 0
    ax = plt.subplot(gs[0, 0])
    add_data(ax, 0, label='2')
    # 1
    ax = plt.subplot(gs[1, 0])
    add_data(ax, 1, label='1')
    # 2
    ax = plt.subplot(gs[2, 0])
    add_data(ax, 2, label='0')
    ax.set_ylabel('amplitude', fontsize=18, labelpad=20)
    # 3
    ax = plt.subplot(gs[3, 0])
    add_data(ax, 3, label='-1')
    # 4
    ax = plt.subplot(gs[4, 0])
    add_data(ax, 4, label='-2', show_x=True)
    # colorbar
    cax = plt.subplot(gs[:, 1])
    label = r'$\mathsf{((d\phi / dt)-\omega_{10})/\Delta_{\omega}}$'
    wt.artists.plot_colorbar(cax=cax, cmap=coolwarm, ticks=[-2, -1, 0, 1, 2],
                             label=label)
    # finish
    plt.savefig(output_path, dpi=300, transparent=True, pad_inches=1)
    plt.close('fig')


### pathway 1 lineshapes ######################################################


output_path = os.path.join(directory, 'pw1 lineshapes.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # preapre figure
    fig, gs = wt.artists.create_figure(width='double', nrows=2, 
                                       cols=[1, 1, 0.25, 1, 1, 'cbar'])
    # data
    title = None        
    filepath = os.path.join(directory, 'measured', 'smear 0.0', 'dpr 1.0 TO 1 w1 w2 d1 d2.hdf5')
    d = wt.kit.read_h5(filepath)
    arr = np.sqrt(d['arr'])
    # 2D delay
    ax = plt.subplot(gs[0:2, 0:2])
    xi = -d['d1']/50.
    yi = -d['d2']/50.
    zi = arr[20, 20].copy().T
    xi, yi, zi = zoom_arrs(xi, yi, zi)
    zi /= zi.max()
    levels = np.linspace(0, 1, 200)
    ax.contourf(xi, yi, zi, cmap=cmap, levels=levels)
    ax.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
    ax.grid()
    ax.axhline(0, c='k', ls='-', lw=2)
    ax.axvline(0, c='k', ls='-', lw=2)
    wt.artists.diagonal_line(xi, yi, ax=ax, c='k', ls='-', lw=2)
    ax.set_xticks(delay_ticks)
    ax.set_yticks(delay_ticks)
    ax.set_xlabel(delay_label(1), fontsize=18)
    ax.set_ylabel(delay_label(2), fontsize=18)
    ax.scatter([0], [2], s=500, facecolor='m', edgecolor='w', lw=2, zorder=10)
    ax.scatter([-2], [2], s=500, facecolor='cyan', edgecolor='w', lw=2, zorder=10)
    ax.scatter([0], [0], s=500, facecolor='orange', edgecolor='w', lw=2, zorder=10)
    ax.scatter([-2], [0], s=500, facecolor='pink', edgecolor='w', lw=2, zorder=10)
    ax.set_xlim(xi.min(), xi.max())
    ax.set_ylim(yi.min(), yi.max())
    ax.set_title('$\mathsf{\omega_1 = \omega_2 = \omega_{10}}$')
    adjust_spines(ax, 'g')
    # preapre for frequency plots
    xi = normalize_frequency(d['w1'])
    yi = normalize_frequency(d['w2'])
    xi = scipy.ndimage.interpolation.zoom(xi, 5)
    yi = scipy.ndimage.interpolation.zoom(yi, 5)
    freq_ticks = [-2, 0, 2]
    def add_axs(r, c, bg_c='g'):
        ax = plt.subplot(gs[r, c], frameon=False)
        for key in ['bottom', 'top', 'right', 'left']:
            ax.spines[key].set_color(bg_c)
            ax.spines[key].set_linewidth(6)
            ax.spines[key].zorder = 1000
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        newax = fig.add_axes(ax.get_position(), frameon=True)
        return newax
    # UL
    ax = add_axs(0, 3)
    adjust_spines(ax, 'cyan', lw=3)
    zi = arr[:, :, 15, 5].T
    zi /= zi.max()
    zi = scipy.ndimage.interpolation.zoom(zi, 5)
    levels = np.linspace(0, 1, 200)
    ax.contourf(xi, yi, zi, cmap=cmap, levels=levels)
    ax.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    ax.set_xticks(freq_ticks)
    ax.set_yticks(freq_ticks)
    ax.grid()
    wt.artists.diagonal_line(xi, yi)
    ax.set_ylabel(frequency_label('2'))
    plt.setp(ax.get_xticklabels(), visible=False)
    # UR
    ax = add_axs(0, 4)
    adjust_spines(ax, 'm', lw=3)
    zi = arr[:, :, 10, 5].T
    zi /= zi.max()
    zi = scipy.ndimage.interpolation.zoom(zi, 5)
    ax.contourf(xi, yi, zi, cmap=cmap, levels=levels)
    ax.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    ax.set_xticks(freq_ticks)
    ax.set_yticks(freq_ticks)
    ax.grid()
    wt.artists.diagonal_line(xi, yi)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    # LL
    ax = add_axs(1, 3)
    adjust_spines(ax, 'pink', lw=3)
    zi = arr[:, :, 15, 10].T
    zi /= zi.max()
    zi = scipy.ndimage.interpolation.zoom(zi, 5)
    levels = np.linspace(0, 1, 200)
    ax.contourf(xi, yi, zi, cmap=cmap, levels=levels)
    ax.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    ax.set_xticks(freq_ticks)
    ax.set_yticks(freq_ticks)
    ax.grid()
    wt.artists.diagonal_line(xi, yi)
    ax.set_ylabel(frequency_label('2'))
    ax.set_xlabel(frequency_label('1'))
    # LR
    ax = add_axs(1, 4)
    adjust_spines(ax, 'orange', lw=3)
    zi = arr[:, :, 10, 10].T
    zi /= zi.max()
    zi = scipy.ndimage.interpolation.zoom(zi, 5)
    levels = np.linspace(0, 1, 200)
    ax.contourf(xi, yi, zi, cmap=cmap, levels=levels)
    ax.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    ax.set_xticks(freq_ticks)
    ax.set_yticks(freq_ticks)
    ax.grid()
    wt.artists.diagonal_line(xi, yi)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.set_xlabel(frequency_label('1'))
    # colorbar
    ticks = np.linspace(0, 1, 11)
    cax = plt.subplot(gs[:, -1])
    matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, ticks=ticks)
    cax.set_ylabel('amplitude', fontsize=18)
    # finish
    if title:
        fig.suptitle(title, fontsize=20)
    plt.savefig(output_path, dpi=300, transparent=True, pad_inches=1)
    plt.close('fig')


### delay space ratios ########################################################


output_path = os.path.join(directory, 'delay space ratios.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # prepare figure
    aspects = [[[0, 0], 0.25]]
    fig, gs = wt.artists.create_figure(width='double', nrows=2, 
                                       cols=[1, 1, 1, 'cbar'], aspects=aspects,
                                       hspace=0.1)
    # subplot method
    def plot(gs_col, dpr, manuals, yticks, title='', c='grey'):
        sps = gs[1, gs_col]
        ax = plt.subplot(sps)
        if not yticks:
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            ax.set_ylabel(delay_label(2), fontsize=18)
        cmap = wt.artists.colormaps['default']
        # get data from zip
        template_fname = os.path.join('measured', 'smear 0.0', 'dpr {0} TO {1} w1 w2 d1 d2.hdf5')
        filepath = template_fname.format(dpr, 'all')
        npz = wt.kit.read_h5(filepath)
        d1 = -npz['d1']/50.
        d2 = -npz['d2']/50.
        arr = np.sqrt(npz['arr'][20, 20].T)
        TO_all = arr.copy()
        amps_max = arr.max()
        arr /= arr.max()
        amps = arr.copy()
        # plot amplitude
        levels = np.linspace(0, 1, 200)
        ax.contourf(d1, d2, arr, levels=levels, cmap=cmap)
        ax.set_xlim(d1.min(), d1.max())
        ax.set_ylim(d2.min(), d2.max())    
        ax.grid()
        ax.axhline(0, c='k', ls='-', lw=2)
        ax.axvline(0, c='k', ls='-', lw=2)
        wt.artists.diagonal_line(d1, d2, ax=ax, c='k', ls='-', lw=2)
        def plot_contours(TO, manual):
            filepath = template_fname.format(dpr, TO)
            npz = wt.kit.read_h5(filepath)
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
            levels = [0.5, 0.75, 0.9, 0.99]
            CS = ax.contour(d1, d2, arr, colors='k', levels=levels, alpha=0.5)
            plt.clabel(CS, inline=True, fontsize=12, fmt='%.2f', manual=manual, c='k')
        plot_contours(1, manuals[0])
        plot_contours(3, manuals[1])
        plot_contours(5, manuals[2])
        plot_contours(6, manuals[3])
        ax.set_xticks(delay_ticks)
        ax.set_yticks(delay_ticks)
        ax.set_xlabel(delay_label(1), fontsize=18)
        adjust_spines(ax, c)
        # plot slice
        ax = plt.subplot(gs[0, gs_col])
        ax.axvline(0, c='k', ls='-', lw=2, zorder=50)
        xi = d1
        yi = amps[10]
        ax.plot(xi, yi, c='grey', lw=5)
        ax.set_yticks([0, 0.5, 1.])
        ax.set_ylim(0, 1.1)
        plt.setp(ax.get_xticklabels(), visible=False)
        if not yticks:
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            ax.set_ylabel('amp.', fontsize=18)
        ax.grid()
        def plot_slice(TO, c):
            filepath = template_fname.format(dpr, TO)
            npz = wt.kit.read_h5(filepath)
            arr = np.sqrt(npz['arr'][20, 20].T)
            arr /= amps_max
            ax.plot(xi, arr[10], lw=2, c=c)
        plot_slice(5, 'm')
        plot_slice(6, 'orange')
        plot_slice(1, 'c')
        plot_slice(3, 'c')
        ax.set_title(title, fontsize=18, y=1.08)
        adjust_spines(ax, c)
    manuals = [None, None, None, None]
    manuals[0] = [[(-0.31, 0), (-1, 1), (-1.2, 1.5), (-2, 3)],
                  [(-0.2, -0.7), (-2.5, -1.3), (-3, -1.9), (-4, -2.5)],
                  [(-0.7, -0.02), (-1.3, -2.5), (-1.9, -3), (-2.5, -4)],
                  [(0, -0.31), (1, -1), (1.5, -1.2), (3, -2)]]
    manuals[1] = [[(-0.31, 0), (-1, 1), (-1.2, 1.5), (-2, 3)],
                  [(-0.2, -0.7), (-2, -1.3), (-3, -1.7), (-4, -2.25)],
                  [(-0.7, -0.02), (-1.3, -2), (-1.7, -3), (-2.25, -4)],
                  [(0, -0.31), (1, -1), (1.5, -1.2), (3, -2)]]
    manuals[2] = [[(-0.31, 0), (-1, 1), (-1.2, 1.5), (-2, 3)],
                  [(-0.2, -0.7), (-2, -1), (-3, -1.5), (-4, -2.05)],
                  [(-0.7, -0.02), (-1, -2), (-1.5, -3), (-2.05, -4)],
                  [(0, -0.31), (1, -1), (1.5, -1.2), (3, -2)]]    
    dprs = ['0.5', '1.0', '2.0']
    ax = plot(0, dprs[0], manuals=manuals[0], yticks=True, title=r'$\mathsf{\Gamma_{10}\Delta_t = 2.0}$', c='b')
    ax = plot(1, dprs[1], manuals=manuals[1], yticks=False, title=r'$\mathsf{\Gamma_{10}\Delta_t = 1.0}$', c='g')
    ax = plot(2, dprs[2], manuals=manuals[2], yticks=False, title=r'$\mathsf{\Gamma_{10}\Delta_t = 0.5}$', c='r')
    # colorbar
    ticks = np.linspace(0, 1, 11)
    cax = plt.subplot(gs[:, -1])
    matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, ticks=ticks)
    cax.set_ylabel('amplitude', fontsize=18)
    # finish
    plt.savefig(output_path, dpi=300, transparent=True, pad_inches=1)
    plt.close('fig')


### spectral evolution ########################################################


output_path = os.path.join(directory, 'spectral evolution.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # prepare figure
    aspects = [[[1, 0], 0.3]]
    fig, gs = wt.artists.create_figure(width='double', aspects=aspects,
                                       cols=[1, 1, 1, 1, 1], nrows=2,
                                       hspace=0.1)
    # plot method
    p_template = os.path.join(directory, 'measured', 'smear 0.0', 'dpr {} TO {} w1 w2 d1 d2.hdf5')
    def plot(col, d2_index, d2_text):
        ax = plt.subplot(gs[0, col])
        dprs = ['2.0', '1.0', '0.5']
        cs = {'0.5': 'b',
              '1.0': 'g',
              '2.0': 'r'}
        maxima = {}
        # 2D freq
        for dpr in dprs:
            p = p_template.format(dpr, 'all')
            d = wt.kit.read_h5(p)
            xi = normalize_frequency(d['w1'])
            yi = normalize_frequency(d['w2'])
            zi = np.sqrt(d['arr'][:, :, 10, d2_index].T)
            maxima[dpr] = zi.max()
            zi /= zi.max()
            xi, yi, zi = zoom_arrs(xi, yi, zi)
            levels = [0, 0.5, 1]
            ax.contour(xi, yi, zi, levels=levels, colors=cs[dpr],
                       linewidths=2, alpha=0.8)
            CS = plt.contourf(xi, yi, zi, levels=levels, colors='k', alpha=0.0)
            zc = CS.collections[1]
            plt.setp(zc, alpha=0.1)
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)
        ticks = [-1, 0, 1]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top') 
        ax.set_xlabel(frequency_label('1'), fontsize=18)
        if col == 0:
            ax.set_ylabel(frequency_label('2'), fontsize=18)
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        ax.grid()
        label = r'$\mathsf{\tau_{21}=' + d2_text + r'\Delta_t}$'
        wt.artists.corner_text(label, fontsize=18, background_alpha=1)
        wt.artists.diagonal_line(xi, yi)
        # histogram
        ax = plt.subplot(gs[1, col])
        TOs = ['1', '3', '5']
        for i, TO in enumerate(TOs):
            for j, dpr in enumerate(dprs):
                p = p_template.format(dpr, TO)
                d = wt.kit.read_h5(p)
                fraction = np.sqrt(d['arr'][20, 20, 10, d2_index])
                fraction /= maxima[dpr]
                position = i - (j-1)/4 - 1/8
                if TO == '5':
                    fraction *= 2
                plt.bar(position, fraction, color=cs[dpr], width=0.25)
        ax.set_ylim(0, 1.1)
        ax.set_yticks([0, 0.5, 1])
        if col == 0:
            ax.set_ylabel('amp.', fontsize=18)
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        ax.grid()
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(['I', 'III', 'V,VI'], fontsize=18)
    # fill out plot
    plot(0, d2_index=5, d2_text=r'2.0')
    plot(1, d2_index=9, d2_text=r'0.4')
    plot(2, d2_index=10, d2_text=r'0.0')
    plot(3, d2_index=11, d2_text=r'-0.4')
    plot(4, d2_index=15, d2_text=r'-2.0')
    # finish
    wt.artists.savefig(output_path, fig=fig, close=True)    


### wigners ###################################################################


output_path = os.path.join(directory, 'wigners.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    cmap = wt.artists.colormaps['default']
    def plot(ax, p, w2_position):
        d = wt.kit.read_h5(p)
        w2_d = {'rr': 12, 'r':16, 'c': 20, 'b':24, 'bb': 28}    
        w2_index = w2_d[w2_position]
        arr = d['arr']
        arr = arr [:, :, ::-1, ::-1]  # reverse delay axes
        d2 = d['d2']
        w2 = d['w2']
        w1 = d['w1']
        d2 = d['d1']
        xi = normalize_frequency(w1)
        yi = normalize_delay(d2)
        zi = np.sqrt(arr[:, w2_index, 11, :].T)
        xi, yi, zi = zoom_arrs(xi, yi, zi)
        levels = np.linspace(0, zi.max(), 200)
        ax.contourf(xi, yi, zi, levels=levels, cmap=cmap)
        levels = np.linspace(0, zi.max(), 11)[1:-1]
        ax.contour(xi, yi, zi, levels=levels, colors='k', alpha=0.5)
        ax.axhline(0, c='k', lw=2)
        ax.axvline(normalize_frequency(w2[w2_index]), lw=4, alpha=0.25, c='k')
        ax.grid()
        ax.set_xticks([-2, 0, 2])
        ax.set_xlim(-3, 3)
        ax.set_yticks(delay_ticks)
        adjust_spines(ax, 'g')
    fig, gs = wt.artists.create_figure(nrows=1, cols=[1, 1, 1, 1, 1, 'cbar'], width='double')
    p = os.path.join('measured', 'smear 0.0', 'dpr 1.0 TO all w1 w2 d1 d2.hdf5')
    # red
    ax = plt.subplot(gs[0, 0])
    plot(ax, p, 'rr')
    ax.set_xlabel(frequency_label(1), fontsize=18)
    ax.set_ylabel(delay_label(2), fontsize=18)
    # red
    ax = plt.subplot(gs[0, 1])
    plot(ax, p, 'r')
    ax.set_xlabel(frequency_label(1), fontsize=18)
    plt.setp(ax.get_yticklabels(), visible=False)
    # center
    ax = plt.subplot(gs[0, 2])
    plot(ax, p, 'c')
    ax.set_xlabel(frequency_label(1), fontsize=18)
    plt.setp(ax.get_yticklabels(), visible=False)
    # blue
    ax = plt.subplot(gs[0, 3])
    plot(ax, p, 'b')
    ax.set_xlabel(frequency_label(1), fontsize=18)
    plt.setp(ax.get_yticklabels(), visible=False)
    # blue
    ax = plt.subplot(gs[0, 4])
    plot(ax, p, 'bb')
    ax.set_xlabel(frequency_label(1), fontsize=18)
    plt.setp(ax.get_yticklabels(), visible=False)
    # finish ------------------------------------------------------------------
    # cmap
    cax = plt.subplot(gs[:, -1])
    wt.artists.plot_colorbar(cax=cax, label='amplitude')
    # finish
    plt.savefig('wigners.png', dpi=300, transparent=True, pad_inches=1)
    plt.close(fig)


### inhom delay space ratios ##################################################


output_path = os.path.join(directory, 'inhom delay space ratios.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # prepare figure
    aspects = [[[0, 0], 0.25]]
    fig, gs = wt.artists.create_figure(width='single', nrows=2, 
                                       cols=[1, 'cbar'], aspects=aspects,
                                       hspace=0.1)
    # subplot method
    def plot(gs_col, dpr, manuals, yticks):
        sps = gs[1, gs_col]
        ax = plt.subplot(sps)
        if not yticks:
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            ax.set_ylabel(delay_label(2), fontsize=18)
        cmap = wt.artists.colormaps['default']
        # get data from zip
        template_fname = os.path.join('measured', 'smear 1.0', 'dpr {0} TO {1} w1 w2 d1 d2.hdf5')
        filepath = template_fname.format(dpr, 'all')
        d = wt.kit.read_h5(filepath)
        d1 = -d['d1']/50.
        d2 = -d['d2']/50.
        arr = np.sqrt(d['arr'][20, 20].T)
        TO_all = arr.copy()
        amps_max = arr.max()
        arr /= arr.max()
        amps = arr.copy()
        # plot amplitude
        levels = np.linspace(0, 1, 200)
        ax.contourf(d1, d2, arr, levels=levels, cmap=cmap)
        ax.set_xlim(d1.min(), d1.max())
        ax.set_ylim(d2.min(), d2.max())    
        ax.grid()
        ax.axhline(0, c='k', ls='-', lw=2)
        ax.axvline(0, c='k', ls='-', lw=2)
        wt.artists.diagonal_line(d1, d2, ax=ax, c='k', ls='-', lw=2)
        def plot_contours(TO, manual):
            filepath = template_fname.format(dpr, TO)
            npz = wt.kit.read_h5(filepath)
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
            levels = [0.5, 0.75, 0.9, 0.99]
            CS = ax.contour(d1, d2, arr, colors='k', levels=levels, alpha=0.5)
            plt.clabel(CS, inline=True, fontsize=12, fmt='%.2f', manual=manual, c='k')
        plot_contours(1, manuals[0])
        plot_contours(3, manuals[1])
        plot_contours(5, manuals[2])
        plot_contours(6, manuals[3])
        ax.set_xticks(delay_ticks)
        ax.set_yticks(delay_ticks)
        ax.set_xlabel(delay_label(1), fontsize=18)
        adjust_spines(ax, 'g')
        # plot slice
        ax = plt.subplot(gs[0, gs_col])
        ax.axvline(0, c='k', ls='-', lw=2, zorder=50)
        xi = d1
        yi = amps[10]
        ax.plot(xi, yi, c='grey', lw=5)
        ax.set_yticks([0, 0.5, 1.])
        ax.set_ylim(0, 1.1)
        plt.setp(ax.get_xticklabels(), visible=False)
        if not yticks:
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            ax.set_ylabel('amp.', fontsize=18)
        ax.grid()
        def plot_slice(TO, c, ls='-'):
            filepath = template_fname.format(dpr, TO)
            npz = wt.kit.read_h5(filepath)
            arr = np.sqrt(npz['arr'][20, 20].T)
            arr /= amps_max
            ax.plot(xi, arr[10], lw=2, c=c, ls=ls)
        plot_slice(5, 'm')
        plot_slice(6, 'orange')
        plot_slice(1, 'c')
        plot_slice(3, 'c', ls='--')
        adjust_spines(ax, 'g')
    manual = [[(-0.31, 0), (-1, 1), (-1.2, 2.1), (-2, 3)],
              [(-0.2, -0.7), (-2, -1), (-3, -1.5), (-4, -2.05)],
              [(-0.7, -0.02), (-1, -2), (-1.5, -3), (-2.05, -4)],
              [(0, -0.31), (1, -1), (2.1, -1.2), (3, -2)]]    
    plot(0, '1.0', manuals=manual, yticks=True)
    # colorbar
    ticks = np.linspace(0, 1, 11)
    cax = plt.subplot(gs[:, -1])
    matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, ticks=ticks)
    cax.set_ylabel('amplitude', fontsize=18)
    # finish
    plt.savefig(output_path, dpi=300, transparent=True, pad_inches=1)
    plt.close('fig')


### inhom spectral evolution ##################################################


output_path = os.path.join(directory, 'inhom spectral evolution.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # prepare figure
    aspects = [[[1, 0], 0.3]]
    fig, gs = wt.artists.create_figure(width='double', aspects=aspects,
                                       cols=[1, 1, 1, 1, 1], nrows=2,
                                       hspace=0.1)
    # plot method
    p_template = os.path.join(directory, 'measured', 'smear {}', 'dpr {} TO {} w1 w2 d1 d2.hdf5')
    def plot(col, d2_index, d2_text):
        ax = plt.subplot(gs[0, col])
        dprs = ['2.0', '1.0', '0.5']
        cs = {'0.5': 'b',
              '1.0': 'g',
              '2.0': 'r'}
        smear = {'0.5': '2.0',
                 '1.0': '1.0',
                 '2.0': '0.5'}
        maxima = {}
        # 2D freq
        for dpr in dprs:
            p = p_template.format(smear[dpr], dpr, 'all')  # smear = dpr
            d = wt.kit.read_h5(p)
            xi = normalize_frequency(d['w1'])
            yi = normalize_frequency(d['w2'])
            zi = np.sqrt(d['arr'][:, :, 10, d2_index].T)
            maxima[dpr] = zi.max()
            zi /= zi.max()
            xi, yi, zi = zoom_arrs(xi, yi, zi)
            levels = [0, 0.5, 1]
            ax.contour(xi, yi, zi, levels=levels, colors=cs[dpr],
                       linewidths=2, alpha=0.8)
            CS = plt.contourf(xi, yi, zi, levels=levels, colors='k', alpha=0.0)
            zc = CS.collections[1]
            plt.setp(zc, alpha=0.1)
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)
        ticks = [-1, 0, 1]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top') 
        ax.set_xlabel(frequency_label('1'), fontsize=18)
        if col == 0:
            ax.set_ylabel(frequency_label('2'), fontsize=18)
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        ax.grid()
        label = r'$\mathsf{\tau_{21}=' + d2_text + r'\Delta_t}$'
        wt.artists.corner_text(label, fontsize=18, background_alpha=1)
        wt.artists.diagonal_line(xi, yi)
        # histogram
        ax = plt.subplot(gs[1, col])
        TOs = ['1', '3', '5', '6']
        for i, TO in enumerate(TOs):
            for j, dpr in enumerate(dprs):
                p = p_template.format(smear[dpr], dpr, TO)
                d = wt.kit.read_h5(p)
                fraction = np.sqrt(d['arr'][20, 20, 10, d2_index])
                fraction /= maxima[dpr]
                position = i - (j-1)/4 - 1/8
                plt.bar(position, fraction, color=cs[dpr], width=0.25)
        ax.set_ylim(0, 1.1)
        ax.set_yticks([0, 0.5, 1])
        if col == 0:
            ax.set_ylabel('amp.', fontsize=18)
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        ax.grid()
        ax.set_xticks([0, 1, 2, 3])
        ax.set_xticklabels(['I', 'III', 'V', 'VI'], fontsize=18)
    # fill out plot
    plot(0, d2_index=5, d2_text=r'2.0')
    plot(1, d2_index=9, d2_text=r'0.4')
    plot(2, d2_index=10, d2_text=r'0.0')
    plot(3, d2_index=11, d2_text=r'-0.4')
    plot(4, d2_index=15, d2_text=r'-2.0')
    # finish
    wt.artists.savefig(output_path, fig=fig, close=True)


### steady state ##############################################################


output_path = os.path.join(directory, 'steady state.png')

force_plotting = True

if not os.path.isfile(output_path) or force_plotting:
    # create figure
    fig, gs = wt.artists.create_figure(width='double', cols=[1, 1, 1, 1, 1, 1, 'cbar', 0.25, 2], nrows=4)
    def plot(ax, xi, yi, zi, ax_diag, ax_anti, c):
        '''
        accepts data on the amplitude level
        '''
        # 2D
        ax.contourf(xi, yi, zi)
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)
        levels = np.linspace(0, 1, 200)
        plt.contourf(xi, yi, zi, levels=levels, cmap=cmap)
        plt.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
        ax.grid()
        ax.set_xticks([-2, -1, 0, 1, 2])
        ax.set_yticks([-2, -1, 0, 1, 2])
        ax.plot([-2, 2], [-2, 2], c=c, lw=4, ls='--', alpha=0.5)
        ax.plot([-2, 2], [2, -2], c=c, lw=4, ls='-', alpha=0.5)
        adjust_spines(ax, 'g')
        if ax.rowNum == 2:
            ax.set_xlabel(frequency_label('1'), fontsize=18)
        else:
            plt.setp(ax.get_xticklabels(), visible=False)
        if ax.is_first_col():
            ax.set_ylabel(frequency_label('2'), fontsize=18)
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        # diagonal
        diag = np.diag(zi)
        ax_diag.plot(xi, diag, c=c, lw=2, ls='--', alpha=0.5)
        ax_diag.grid()
        ax_diag.set_xlim(-2, 2)
        ax_diag.set_ylim(0, 1.1)
        ax_diag.set_xticks([-2, -1, 0, 1, 2])
        ax_diag.set_yticks([0, 0.5, 1])
        plt.setp(ax_diag.get_xticklabels(), visible=False)
        plt.setp(ax_diag.get_yticklabels(), visible=False)
        adjust_spines(ax_diag, 'g')
        # antidiagonal
        zi = np.flipud(zi)
        diag = np.diag(zi)
        ax_anti.plot(xi, diag, c=c, lw=2, alpha=0.5)
        ax_anti.grid()
        ax_anti.set_xlim(-2, 2)
        ax_anti.set_ylim(0, 1.1)
        ax_anti.set_xticks([-2, -1, 0, 1, 2])
        ax_anti.set_yticks([0, 0.5, 1])
        plt.setp(ax_anti.get_yticklabels(), visible=False)
        if not ax_anti.is_last_row():
            plt.setp(ax_anti.get_xticklabels(), visible=False)
        else:
            ax_anti.set_xlabel(frequency_label('2'), fontsize=18)
        adjust_spines(ax_anti, 'g')
    for delay_index, delay in enumerate(['zero delay', 'pw56']):
        first_row = delay_index*2
        # 1D plots
        ax_diag = plt.subplot(gs[first_row, -1])
        ax_anti = plt.subplot(gs[first_row+1, -1])
        # infinite lifetime
        ax = plt.subplot(gs[first_row:first_row+2, 0:2])
        p = os.path.join('precalculated', 'cw {} infinite lifetime.hdf5'.format(delay))
        d = wt.kit.read_h5(p)
        xi = d['w1']
        yi = d['w2']
        zi = d['arr']
        plot(ax, xi, yi, zi, ax_diag, ax_anti, 'orange')
        if delay_index == 0:
            ax.set_title(r'Steady State, $\Gamma_{11} \rightarrow \infty$', fontsize=18)
        # finite lifetime
        ax = plt.subplot(gs[first_row:first_row+2, 2:4])
        p = os.path.join('precalculated', 'cw {} finite lifetime.hdf5'.format(delay))
        d = wt.kit.read_h5(p)
        xi = d['w1']
        yi = d['w2']
        zi = d['arr']
        plot(ax, xi, yi, zi, ax_diag, ax_anti, 'm')
        if delay_index == 0:
            ax.set_title('Steady State, $\Gamma_{11} \Delta_t = 1$', fontsize=18)
        # actual
        ax = plt.subplot(gs[first_row:first_row+2, 4:6])
        p = os.path.join(directory, 'measured', 'smear 0.0', 'dpr 1.0 TO all w1 w2 d1 d2.hdf5')
        d = wt.kit.read_h5(p)
        xi = normalize_frequency(d['w1'])
        yi = normalize_frequency(d['w2'])
        if delay == 'zero delay':
            d2i = 10
        elif delay == 'pw56':
            d2i = -1
        zi = d['arr'][:, :, 10, d2i].T
        zi **= 0.5
        zi /= zi.max()
        xi, yi, zi = zoom_arrs(xi, yi, zi)
        plot(ax, xi, yi, zi, ax_diag, ax_anti, 'k')
        if delay_index == 0:
            ax.set_title('Actual', fontsize=18)
            ax.text(2.5, 0.5, r'$\mathsf{T = -\tau_{21} = 0}$', rotation=-90, fontsize=18, va='center', transform=ax.transAxes)
        else:
            ax.text(2.5, 0.5, r'$\mathsf{T = -\tau_{21} = 4 \Delta_t}$', rotation=-90, fontsize=18, va='center', transform=ax.transAxes)
    # colorbar
    ticks = np.linspace(0, 1, 11)
    cax = plt.subplot(gs[:, -3])
    matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, ticks=ticks)
    cax.set_ylabel('amplitude', fontsize=18)
    # finish
    wt.artists.savefig(output_path, fig=fig, close=True)


### metrics of correlation ####################################################


output_path = os.path.join(directory, 'metrics.png')

force_plotting = True

if not os.path.isfile(output_path) or force_plotting:
    # create figure
    aspects = [[[0, 0], 0.6], [[1, 0], 0.6], [[2, 0], 0.1]]
    fig, gs = wt.artists.create_figure(width='single', nrows=4, cols=[1, 1, 0], aspects=aspects, hspace=0.33, wspace=0.33)
    for delay_index, delay in enumerate([0., -4.]):
        p = os.path.join('precalculated', 'metrics.txt')
        arr = np.genfromtxt(p)
        # 0 dpr
        # 1 smear
        # 2 delay
        # 3 ellipticity
        # 4 shift
        # organize ------------------------------------------------------------
        delta_gamma = {}
        for key, c in zip(['2.0', '1.0', '0.5'], ['b', 'g', 'r']):
            l = []
            for measurement in arr:
                if measurement[2] == delay and 1/measurement[0] == float(key):
                    l.append(measurement)
            delta_gamma[key] = (l, c)
        delta_inhom_over_gamma = {}
        for key, ls in zip(['0.0', '0.5', '1.0', '2.0'], ['-', ':', '-.', '--']):
            l = []
            for measurement in arr:
                c = measurement[1] * measurement[0]
                if measurement[2] == delay and c == float(key):
                    l.append(measurement)
            delta_inhom_over_gamma[key] = (l, ls)
        # 3PEPS ---------------------------------------------------------------
        ax = plt.subplot(gs[0, delay_index])
        # points
        for item in delta_gamma.values():
            for measurement in item[0]:
                c = measurement[1] * measurement[0]
                if c == 0.:
                    marker = 'x'
                    s = 20
                else:
                    marker = 'o'
                    s = 100*float(c)                
                ax.scatter(measurement[1], measurement[4], s=s, marker=marker, alpha=0.5, color='k', edgecolor='none')
        # delta gamma
        for key, item in delta_gamma.items():
            measurements, c = item
            arr = np.array(measurements)
            if arr.size > 0:
                ax.plot(arr[:, 1], arr[:, 4], c=c, lw=2, label=key)
        # delta_inhom_over_gamma
        for key, item in delta_inhom_over_gamma.items():
            measurements, ls = item
            arr = np.array(measurements)
            if arr.size > 0:
                lw = (float(key)*10) + 1
                ax.plot(arr[:, 1], arr[:, 4], c='k', ls='-', alpha=0.25, lw=2)
        ax.grid()
        ax.set_yticks([0, 0.5, 1])
        ax.set_xlim(0, 2.2)
        ax.set_ylim(0, 1.1)
        plt.setp(ax.get_xticklabels(), visible=False)
        if delay_index == 1:
            plt.setp(ax.get_yticklabels(), visible=False)
            ax.set_title(r'$\mathsf{T = -\tau_{21} = 4 \Delta_t}$')
            legend = plt.legend(bbox_to_anchor=(0.975, 1), loc=2, borderaxespad=0., frameon=False)
            legend.set_title(r'$\mathsf{\Gamma_{10}\Delta_t}$')
        else:
            ax.set_ylabel(r'$\mathsf{\Delta_{\mathrm{3PEPS}}/\Delta_t}$', fontsize=18)
            ax.set_title(r'$\mathsf{T = -\tau_{21} = 0}$')
        # ellipticity ---------------------------------------------------------
        ax = plt.subplot(gs[1, delay_index])
        # points
        for item in delta_gamma.values():
            for measurement in item[0]:
                c = measurement[1] * measurement[0]
                if c == 0.:
                    marker = 'x'
                    s = 20
                else:
                    marker = 'o'
                    s = 100*float(c)                
                ax.scatter(measurement[1], measurement[3], s=s, marker=marker, alpha=0.5, color='k', edgecolor='none')
        # delta gamma
        for key, item in delta_gamma.items():
            measurements, c = item
            arr = np.array(measurements)
            if arr.size > 0:
                ax.plot(arr[:, 1], arr[:, 3], c=c, lw=2)
        # delta_inhom_over_gamma
        for key, item in delta_inhom_over_gamma.items():
            measurements, ls = item
            arr = np.array(measurements)
            if arr.size > 0:
                ax.plot(arr[:, 1], arr[:, 3], c='k', ls='-', alpha=0.25, lw=2)
        ax.grid()
        ax.set_xlabel(r'$\mathsf{\Delta_{\mathrm{inhom}}/\Delta_\omega}$', fontsize=18)
        ax.set_yticks([0, 0.5, 1])
        ax.set_xlim(0, 2.2)
        ax.set_ylim(0, 1.1)
        if delay_index == 1:
            plt.setp(ax.get_yticklabels(), visible=False)
            # fake a legend
            for c in ['0.0', '0.25', '0.5', '1.0', '2.0', '4.0']:
                if c == '0.0':
                    marker = 'x'
                    s = 20
                else:
                    marker = 'o'
                    s = 100*float(c)
                ax.scatter(-1, -1, s=s, marker=marker, alpha=0.5, color='k', edgecolor='none', label=c)
            legend = plt.legend(bbox_to_anchor=(0.975, 1), loc=2, borderaxespad=0., frameon=False, scatterpoints=1)
            legend.set_title(r'$\mathsf{\Delta_{\mathrm{inhom}}/\Gamma_{10}}$')
        else:
            ax.set_ylabel(r'ellipticity', fontsize=18)
        # parametric ----------------------------------------------------------
        ax = plt.subplot(gs[3, delay_index])
        # points
        for item in delta_gamma.values():
            for measurement in item[0]:
                c = measurement[1] * measurement[0]
                if c == 0.:
                    marker = 'x'
                    s = 20
                else:
                    marker = 'o'
                    s = 100*float(c) 
                ax.scatter(measurement[3], measurement[4], s=s, marker=marker, alpha=0.5, color='k', edgecolor='none')
        # delta gamma
        for key, item in delta_gamma.items():
            measurements, c = item
            arr = np.array(measurements)
            if arr.size > 0:
                ax.plot(arr[:, 3], arr[:, 4], c=c, lw=2)
        # delta_inhom_over_gamma
        for key, item in delta_inhom_over_gamma.items():
            measurements, ls = item
            arr = np.array(measurements)
            if arr.size > 0:
                ax.plot(arr[:, 3], arr[:, 4], c='k', ls='-', alpha=0.25, lw=2)
        ax.grid()
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel(r'ellipticity', fontsize=18)
        if delay_index == 0:
            ax.set_ylabel(r'$\mathsf{\Delta_{\mathrm{3PEPS}}/\Delta_t}$', fontsize=18)
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
    # finish
    wt.artists.savefig(output_path, fig=fig, close=True)


### convolve ##################################################################


output_path = os.path.join(directory, 'convolve.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # create figure -----------------------------------------------------------
    fig, gs = wt.artists.create_figure(width='double', cols=[1, 1.5, 1, 'cbar'])
    # measured 2D frequency with homo lineshape -------------------------------
    ax = plt.subplot(gs[0, 0])
    p = os.path.join(directory, 'measured', 'smear 0.0', 'dpr 1.0 TO all w1 w2 d1 d2.hdf5')
    d = wt.kit.read_h5(p)    
    xi = normalize_frequency(d['w1'])
    yi = normalize_frequency(d['w2'])
    zi = np.sqrt(d['arr'][:, :, -6, -1]).T
    zi /= zi.max()
    levels= np.linspace(0, 1, 200)
    xi, yi, zi = zoom_arrs(xi, yi, zi)
    mappable = ax.contourf(xi, yi, zi, cmap=cmap, levels=levels)
    ax.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
    ax.grid()
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_xticks(freq_ticks)
    ax.set_yticks(freq_ticks)
    wt.artists.diagonal_line(xi, yi)
    wt.artists.corner_text('a', ax=ax, fontsize=18, background_alpha=1)
    ax.set_xlabel(frequency_label('1'), fontsize=18)
    ax.set_ylabel(frequency_label('2'), fontsize=18)
    adjust_spines(ax, 'g')
    # representation of kernal ------------------------------------------------
    ax = plt.subplot(gs[0, 1], projection='3d')
    ax.view_init(35, 75 )
    plt.setp(ax.get_zticklabels(), visible=False)
    # box
    ax.plot([-2, -2], [-2, 2], [0, 0], c='g')
    ax.plot([-2, 2], [2, 2], [0, 0], c='g')
    ax.plot([2, 2], [2, -2], [0, 0], c='g')
    ax.plot([-2, 2], [-2, -2], [0, 0], c='g')
    # gaussian
    xi = np.linspace(3, -3, 1000)
    yi = xi
    function = wt.fit.Gaussian()
    FWHM = 2
    sigma = FWHM/2.3548
    zi = function.evaluate([0, sigma, 1, 0], xi)
    verts = [list(zip(xi, yi, zi))]
    poly = Poly3DCollection(verts, facecolors=[[0., 0., 0., 0.1]])
    ax.add_collection3d(poly)
    ax.set_xlabel(r'$\mathsf{\omega_1}$', fontsize=18)
    ax.set_xlim3d(xi.min(), xi.max())
    ax.yaxis.set_rotate_label(False)  # disable automatic rotation
    ax.set_ylabel(r'$\mathsf{\omega_2}$', fontsize=18, rotation=115)
    ax.set_ylim3d(xi.min(), xi.max())
    ax.set_zlim3d(0, 1)
    wt.artists.corner_text('b', ax=ax, fontsize=18, background_alpha=1)
    labels = ['', 2, '', 0,  '', -2, '']
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    # lines
    ax.plot([1, 1], [1, 1], [0, 0.5], c='k', ls=':')
    ax.plot([-1, -1], [-1, -1], [0, 0.5], c='k', ls=':')
    ax.plot([1, -1], [1, -1], [0.5, 0.5], c='k', ls=':')
    # inhomogenious -----------------------------------------------------------
    ax = plt.subplot(gs[0, 2])
    p = os.path.join(directory, 'measured', 'smear 2.0', 'dpr 1.0 TO all w1 w2 d1 d2.hdf5')
    d = wt.kit.read_h5(p)    
    xi = normalize_frequency(d['w1'])
    yi = normalize_frequency(d['w2'])
    zi = np.sqrt(d['arr'][:, :, -6, -1]).T
    zi /= zi.max()
    levels= np.linspace(0, 1, 200)
    xi, yi, zi = zoom_arrs(xi, yi, zi)
    mappable = ax.contourf(xi, yi, zi, cmap=cmap, levels=levels)
    ax.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
    ax.grid()
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_xticks(freq_ticks)
    ax.set_yticks(freq_ticks)
    wt.artists.diagonal_line(xi, yi)
    wt.artists.corner_text('c', ax=ax, fontsize=18, background_alpha=1)
    # show inhomo linewidth applied
    center = 7000.
    width = 2*delta_omega
    l = np.linspace(center-(width/2), center+(width/2), 100)
    l = normalize_frequency(l)
    ax.plot(l, l, c='k', lw=5)
    ax.set_xlabel(frequency_label('1'), fontsize=18)
    ax.set_ylabel(frequency_label('2'), fontsize=18)
    adjust_spines(ax, 'g')
    # colorbar ----------------------------------------------------------------
    ticks = np.linspace(0, 1, 11)
    cax = plt.subplot(gs[:, -1])
    matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, ticks=ticks)
    cax.set_ylabel('amplitude', fontsize=18)
    # finish
    plt.savefig(output_path, dpi=300, transparent=True, pad_inches=1)
    plt.close('fig')


### 2D delays #################################################################


output_path = os.path.join(directory, '2D delays.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # set up figure
    fig, gs = wt.artists.create_figure(width='double', nrows=3,
                                       cols=[1, 1, 1, 1, 0, 'cbar'])
    # fit method
    def measure_shifts(xi, yi, zi):
        '''
        3PEPS shift defined using sig on the intensity level
        '''
        # shape arrays
        xi = xi.copy()
        yi = yi.copy()
        zi = zi.copy()
        zi_shifted = np.full(zi.shape, np.nan)
        for i in range(xi.size):
            if i <= 10:
                zi_shifted[:, i] = zi[:, i]
            else:
                for j in range(yi.size):
                    try:
                        zi_shifted[j, i] = zi[j+(i-10), i]
                    except IndexError:
                        pass
        # go by fitting
        def gauss(p, x):
            A, mu, w = p
            z = (mu - x) / (np.sqrt(2) * w)
            out = A * np.exp(-z**2)
            return out
        def erf1(p, x, y):
            return y - gauss(p,x)
        shifts = np.full((11, 2), np.nan)
        for i in range(10, 21):
            sig = zi_shifted[i]
            p0 = []
            p0.append(np.nanmax(sig))  # amplitude
            p0.append(0.)  # center
            p0.append(1.)  # width
            args = wt.kit.remove_nans_1D([xi, sig])
            params = leastsq(erf1, p0, args=tuple(args), full_output=False)[0]
            shifts[i-10] = [yi[i], params[1]]
        for i in range(7, 11):
            shifts[i, 1] = shifts[7, 1]  # edge is corrupted, esentially
        return shifts
    # worker method
    p_template = os.path.join(directory, 'measured', 'smear {}', 'dpr {} TO all w1 w2 d1 d2.hdf5')
    cs = {'0.5': 'b',
          '1.0': 'g',
          '2.0': 'r'}    
    def plot(ax, dpr, smear):
        p = p_template.format(smear, dpr)
        d = wt.kit.read_h5(p)
        xi = normalize_delay(-d['d1'])
        yi = normalize_delay(-d['d2'])
        zi = d['arr'][20, 20].T
        shifts = measure_shifts(xi, yi, zi)  # intensity level
        zi = np.sqrt(zi)
        xi, yi, zi = zoom_arrs(xi, yi, zi)
        zi /= zi.max()
        levels = np.linspace(0, 1, 200)
        plt.contourf(xi, yi, zi, levels=levels, cmap=cmap)
        plt.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
        adjust_spines(ax, cs[dpr])
        ax.grid()
        wt.artists.diagonal_line(xi, yi, lw=2, ls='-')
        ax.axhline(0, c='k', lw=2)
        ax.axvline(0, c='k', lw=2)
        if ax.is_first_col():
            ax.set_ylabel(delay_label(2), fontsize=18)
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        if ax.is_last_row():
            ax.set_xlabel(delay_label(1), fontsize=18)
        else:
            plt.setp(ax.get_xticklabels(), visible=False)        
        ax.set_xlim(xi.min(), xi.max())
        ax.set_ylim(yi.min(), yi.max())   
        ticks = [-3, -2, -1, 0, 1, 2, 3]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        # shifts
        tau = shifts[:, 1]  # negative
        T = shifts[:, 0] + tau
        tau[tau>0] = -1e-6  # hack for display
        c = 'y'
        ax.plot(tau, T, lw=5, c=c, zorder=1000)
        ax.text(tau[0], T[0], '{:.2f}'.format(abs(tau[0])) , color=c, ha='right', va='bottom', fontsize=18, path_effects=[PathEffects.withStroke(linewidth=2, foreground="k")], zorder=1001)
        ax.text(tau[-1]-0.2, -3.5, '{:.2f}'.format(abs(tau[-1])) , color=c, ha='right', va='center', fontsize=18, path_effects=[PathEffects.withStroke(linewidth=2, foreground="k")], zorder=1001)
    # fill out
    # col 1
    ax = plt.subplot(gs[0, 0])
    plot(ax, '0.5', '0.0')
    ax.set_title(r'$\mathsf{\Delta_{inhom} = 0.0}$', fontsize=18)
    ax = plt.subplot(gs[1, 0])
    plot(ax, '1.0', '0.0')
    ax = plt.subplot(gs[2, 0])
    plot(ax, '2.0', '0.0')
    # col 2
    ax = plt.subplot(gs[0, 1])
    plot(ax, '0.5', '0.5')
    ax.set_title(r'$\mathsf{\Delta_{inhom} = 0.5 \Delta_\omega}$', fontsize=18)
    ax = plt.subplot(gs[1, 1])
    plot(ax, '1.0', '0.5')
    ax = plt.subplot(gs[2, 1])
    plot(ax, '2.0', '0.5')
    # col 3
    ax = plt.subplot(gs[0, 2])
    plot(ax, '0.5', '1.0')
    ax.set_title(r'$\mathsf{\Delta_{inhom} = 1.0 \Delta_\omega}$', fontsize=18)
    ax = plt.subplot(gs[1, 2])
    plot(ax, '1.0', '1.0')
    ax = plt.subplot(gs[2, 2])
    plot(ax, '2.0', '1.0')
    # col 4
    ax = plt.subplot(gs[0, 3])
    plot(ax, '0.5', '2.0')
    ax.set_title(r'$\mathsf{\Delta_{inhom} = 2.0 \Delta_\omega}$', fontsize=18)
    ax.text(1.01, 0.5, r'$\mathsf{\Gamma_{10}\Delta_t = 2.0}$', rotation=-90, fontsize=18, va='center', transform=ax.transAxes)
    ax = plt.subplot(gs[1, 3])
    plot(ax, '1.0', '2.0')
    ax.text(1.01, 0.5, r'$\mathsf{\Gamma_{10}\Delta_t = 1.0}$', rotation=-90, fontsize=18, va='center', transform=ax.transAxes)
    ax = plt.subplot(gs[2, 3])
    plot(ax, '2.0', '2.0')
    ax.text(1.01, 0.5, r'$\mathsf{\Gamma_{10}\Delta_t = 0.5}$', rotation=-90, fontsize=18, va='center', transform=ax.transAxes)
    # colorbar
    cax = plt.subplot(gs[:, -1])
    wt.artists.plot_colorbar(cax=cax, cmap=cmap, label='amplitude', label_fontsize=18)
    # finish
    wt.artists.savefig(output_path, fig=fig, close=True)


### spectral evolution full ###################################################


output_path = os.path.join(directory, 'spectral evolution full.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # set up figure
    fig, gs = wt.artists.create_figure(width='double', nrows=3,
                                       cols=[1, 1, 1, 1, 1, 0, 'cbar'])
    # worker method
    p_template = os.path.join(directory, 'measured', 'smear {}', 'dpr {} TO all w1 w2 d1 d2.hdf5')
    cs = {'0.5': 'b',
          '1.0': 'g',
          '2.0': 'r'}
    labels = {'0.5': r'$\mathsf{\Gamma_{10}\Delta_t = 2.0}$',
              '1.0': r'$\mathsf{\Gamma_{10}\Delta_t = 1.0}$',
              '2.0': r'$\mathsf{\Gamma_{10}\Delta_t = 0.5}$'}
    def plot(ax, dpr, smear, d2_index):
        p = p_template.format(smear, dpr)
        d = wt.kit.read_h5(p)
        xi = normalize_frequency(d['w1'])
        yi = normalize_frequency(d['w2'])
        zi = d['arr'][:, :, 10, d2_index].T
        zi **= 0.5
        xi, yi, zi = zoom_arrs(xi, yi, zi)
        zi /= zi.max()
        levels = np.linspace(0, 1, 200)
        ax.contourf(xi, yi, zi, levels=levels, cmap=cmap)
        ax.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
        ax.grid()
        wt.artists.diagonal_line(xi, yi)
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)
        ticks=[-1, 0, 1]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        if ax.is_first_col():
            ax.set_ylabel(frequency_label(2), fontsize=18)
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        if ax.is_last_row():
            ax.set_xlabel(frequency_label(1), fontsize=18)
        else:
            plt.setp(ax.get_xticklabels(), visible=False)
        if ax.is_first_row():
            d2 = normalize_delay(-d['d2'])[d2_index]
            if d2 == 0.:
                d2 = abs(d2)
            d2_text = '%.1f'%d2
            label = r'$\mathsf{\tau_{21}=' + d2_text + r'\Delta_t}$'
            ax.set_title(label, fontsize=18)
        if d2_index == 15:
            label = labels[dpr]
            ax.text(1.01, 0.5, label, fontsize=18, rotation=-90, ha='left', va='center', transform=ax.transAxes)
        adjust_spines(ax, cs[dpr])
        levels = [0, 0.5, 1]
        ax.contour(xi, yi, zi, levels=levels, colors='k',
                   linewidths=2, alpha=1)
    # fill
    # row 1
    ax = plt.subplot(gs[0, 0])
    plot(ax, '0.5', '0.0', 5)
    ax = plt.subplot(gs[0, 1])
    plot(ax, '0.5', '0.0', 9)
    ax = plt.subplot(gs[0, 2])
    plot(ax, '0.5', '0.0', 10)
    ax = plt.subplot(gs[0, 3])
    plot(ax, '0.5', '0.0', 11)
    ax = plt.subplot(gs[0, 4])
    plot(ax, '0.5', '0.0', 15)
    # row 2
    ax = plt.subplot(gs[1, 0])
    plot(ax, '1.0', '0.0', 5)
    ax = plt.subplot(gs[1, 1])
    plot(ax, '1.0', '0.0', 9)
    ax = plt.subplot(gs[1, 2])
    plot(ax, '1.0', '0.0', 10)
    ax = plt.subplot(gs[1, 3])
    plot(ax, '1.0', '0.0', 11)
    ax = plt.subplot(gs[1, 4])
    plot(ax, '1.0', '0.0', 15)
    # row 3
    ax = plt.subplot(gs[2, 0])
    plot(ax, '2.0', '0.0', 5)
    ax = plt.subplot(gs[2, 1])
    plot(ax, '2.0', '0.0', 9)
    ax = plt.subplot(gs[2, 2])
    plot(ax, '2.0', '0.0', 10)
    ax = plt.subplot(gs[2, 3])
    plot(ax, '2.0', '0.0', 11)
    ax = plt.subplot(gs[2, 4])
    plot(ax, '2.0', '0.0', 15)
    # colorbar
    cax = plt.subplot(gs[:, -1])
    wt.artists.plot_colorbar(cax=cax, cmap=cmap, label='amplitude', label_fontsize=18)

    # finish
    wt.artists.savefig(output_path, fig=fig, close=True)


### wigners full ##############################################################


output_path = os.path.join(directory, 'wigners full.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # set up figure
    fig, gs = wt.artists.create_figure(width='double', nrows=3,
                                       cols=[1, 1, 1, 1, 1, 0, 'cbar'])
    # worker method
    p_template = os.path.join(directory, 'measured', 'smear {}', 'dpr {} TO all w1 w2 d1 d2.hdf5')
    cs = {'0.5': 'b',
          '1.0': 'g',
          '2.0': 'r'}
    labels = {'0.5': r'$\mathsf{\Gamma_{10}\Delta_t = 2.0}$',
              '1.0': r'$\mathsf{\Gamma_{10}\Delta_t = 1.0}$',
              '2.0': r'$\mathsf{\Gamma_{10}\Delta_t = 0.5}$'}
    def plot(ax, dpr, smear, w2_index):
        p = p_template.format(smear, dpr)
        d = wt.kit.read_h5(p)
        xi = normalize_frequency(d['w1'])
        yi = normalize_delay(-d['d2'])
        w2 = normalize_frequency(d['w2'])
        zi = d['arr'][:, w2_index, 10, :].T
        zi **= 0.5
        xi, yi, zi = zoom_arrs(xi, yi, zi)
        zi /= zi.max()
        levels = np.linspace(0, 1, 200)
        ax.contourf(xi, yi, zi, levels=levels, cmap=cmap)
        ax.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
        ax.grid()
        ax.axhline(0, c='k', lw=2)
        ax.axvline(w2[w2_index], lw=4, alpha=0.25, c='k')
        ax.set_xlim(-3, 3)
        ax.set_ylim(-4, 4)
        ticks=[-2, 0, 2]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        if ax.is_first_col():
            ax.set_ylabel(delay_label(2), fontsize=18)
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        if ax.is_last_row():
            ax.set_xlabel(frequency_label(1), fontsize=18)
        else:
            plt.setp(ax.get_xticklabels(), visible=False)
        if w2_index == 28:
            label = labels[dpr]
            ax.text(1.01, 0.5, label, fontsize=18, rotation=-90, ha='left', va='center', transform=ax.transAxes)
        adjust_spines(ax, cs[dpr])
    # fill
    # row 1
    ax = plt.subplot(gs[0, 0])
    plot(ax, '0.5', '0.0', 12)
    ax = plt.subplot(gs[0, 1])
    plot(ax, '0.5', '0.0', 16)
    ax = plt.subplot(gs[0, 2])
    plot(ax, '0.5', '0.0', 20)
    ax = plt.subplot(gs[0, 3])
    plot(ax, '0.5', '0.0', 24)
    ax = plt.subplot(gs[0, 4])
    plot(ax, '0.5', '0.0', 28)
    # row 2
    ax = plt.subplot(gs[1, 0])
    plot(ax, '1.0', '0.0', 12)
    ax = plt.subplot(gs[1, 1])
    plot(ax, '1.0', '0.0', 16)
    ax = plt.subplot(gs[1, 2])
    plot(ax, '1.0', '0.0', 20)
    ax = plt.subplot(gs[1, 3])
    plot(ax, '1.0', '0.0', 24)
    ax = plt.subplot(gs[1, 4])
    plot(ax, '1.0', '0.0', 28)
    # row 3
    ax = plt.subplot(gs[2, 0])
    plot(ax, '2.0', '0.0', 12)
    ax = plt.subplot(gs[2, 1])
    plot(ax, '2.0', '0.0', 16)
    ax = plt.subplot(gs[2, 2])
    plot(ax, '2.0', '0.0', 20)
    ax = plt.subplot(gs[2, 3])
    plot(ax, '2.0', '0.0', 24)
    ax = plt.subplot(gs[2, 4])
    plot(ax, '2.0', '0.0', 28)
    # colorbar
    cax = plt.subplot(gs[:, -1])
    wt.artists.plot_colorbar(cax=cax, cmap=cmap, label='amplitude', label_fontsize=18)

    # finish
    wt.artists.savefig(output_path, fig=fig, close=True)


### inhom spectral evolution full #############################################


output_path = os.path.join(directory, 'inhom spectral evolution full.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # set up figure
    fig, gs = wt.artists.create_figure(width='double', nrows=3,
                                       cols=[1, 1, 1, 1, 1, 0, 'cbar'])
    # worker method
    p_template = os.path.join(directory, 'measured', 'smear {}', 'dpr {} TO all w1 w2 d1 d2.hdf5')
    cs = {'0.5': 'b',
          '1.0': 'g',
          '2.0': 'r'}
    labels = {'0.5': r'$\mathsf{\Gamma_{10}\Delta_t = 2.0}$',
              '1.0': r'$\mathsf{\Gamma_{10}\Delta_t = 1.0}$',
              '2.0': r'$\mathsf{\Gamma_{10}\Delta_t = 0.5}$'}
    def plot(ax, dpr, smear, d2_index):
        p = p_template.format(smear, dpr)
        d = wt.kit.read_h5(p)
        xi = normalize_frequency(d['w1'])
        yi = normalize_frequency(d['w2'])
        zi = d['arr'][:, :, 10, d2_index].T
        zi **= 0.5
        xi, yi, zi = zoom_arrs(xi, yi, zi)
        zi /= zi.max()
        levels = np.linspace(0, 1, 200)
        ax.contourf(xi, yi, zi, levels=levels, cmap=cmap)
        ax.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
        ax.grid()
        wt.artists.diagonal_line(xi, yi)
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)
        ticks=[-1, 0, 1]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        if ax.is_first_col():
            ax.set_ylabel(frequency_label(2), fontsize=18)
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        if ax.is_last_row():
            ax.set_xlabel(frequency_label(1), fontsize=18)
        else:
            plt.setp(ax.get_xticklabels(), visible=False)
        if ax.is_first_row():
            d2 = normalize_delay(-d['d2'])[d2_index]
            if d2 == 0.:
                d2 = abs(d2)
            d2_text = '%.1f'%d2
            label = r'$\mathsf{\tau_{21}=' + d2_text + r'\Delta_t}$'
            ax.set_title(label, fontsize=18)
        if d2_index == 15:
            label = labels[dpr]
            ax.text(1.01, 0.5, label, fontsize=18, rotation=-90, ha='left', va='center', transform=ax.transAxes)
        adjust_spines(ax, cs[dpr])
        levels = [0, 0.5, 1]
        ax.contour(xi, yi, zi, levels=levels, colors='k',
                   linewidths=2, alpha=1)
    # fill
    # row 1
    ax = plt.subplot(gs[0, 0])
    plot(ax, '0.5', '2.0', 5)
    ax = plt.subplot(gs[0, 1])
    plot(ax, '0.5', '2.0', 9)
    ax = plt.subplot(gs[0, 2])
    plot(ax, '0.5', '2.0', 10)
    ax = plt.subplot(gs[0, 3])
    plot(ax, '0.5', '2.0', 11)
    ax = plt.subplot(gs[0, 4])
    plot(ax, '0.5', '2.0', 15)
    # row 2
    ax = plt.subplot(gs[1, 0])
    plot(ax, '1.0', '1.0', 5)
    ax = plt.subplot(gs[1, 1])
    plot(ax, '1.0', '1.0', 9)
    ax = plt.subplot(gs[1, 2])
    plot(ax, '1.0', '1.0', 10)
    ax = plt.subplot(gs[1, 3])
    plot(ax, '1.0', '1.0', 11)
    ax = plt.subplot(gs[1, 4])
    plot(ax, '1.0', '1.0', 15)
    # row 3
    ax = plt.subplot(gs[2, 0])
    plot(ax, '2.0', '0.5', 5)
    ax = plt.subplot(gs[2, 1])
    plot(ax, '2.0', '0.5', 9)
    ax = plt.subplot(gs[2, 2])
    plot(ax, '2.0', '0.5', 10)
    ax = plt.subplot(gs[2, 3])
    plot(ax, '2.0', '0.5', 11)
    ax = plt.subplot(gs[2, 4])
    plot(ax, '2.0', '0.5', 15)
    # colorbar
    cax = plt.subplot(gs[:, -1])
    wt.artists.plot_colorbar(cax=cax, cmap=cmap, label='amplitude', label_fontsize=18)

    # finish
    wt.artists.savefig(output_path, fig=fig, close=True)


### 2D frequencies at zero ####################################################


output_path = os.path.join(directory, '2D frequences at zero.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # set up figure
    fig, gs = wt.artists.create_figure(width='double', nrows=3,
                                       cols=[1, 1, 1, 1, 0, 'cbar'])
    # fit method
    def get_widths(xi, yi, zi):
        '''
        defined on the amplitude level
        '''
        # go by fitting
        def gauss(p, x):
            A, mu, w = p
            z = (mu - x) / (np.sqrt(2) * w)
            out = A * np.exp(-z**2)
            return out
        def erf1(p, x, y):
            return y - gauss(p,x)
        outs = []
        # diagonal
        y = zi.diagonal()
        p0 = []
        p0.append(np.nanmax(y))  # amplitude
        p0.append(0.)  # center
        p0.append(1.)  # width
        out = leastsq(erf1, p0, args=(xi, y), full_output=False)[0]
        outs.append(2.35482*out[2])
        # antidiagonal
        y = zi[...,::-1].diagonal()
        p0 = []
        p0.append(np.nanmax(y))  # amplitude
        p0.append(0.)  # center
        p0.append(1.)  # width
        out = leastsq(erf1, p0, args=(xi, y), full_output=False)[0]
        outs.append(2.35482*out[2])
        return outs
    # worker method
    p_template = os.path.join(directory, 'measured', 'smear {}', 'dpr {} TO all w1 w2 d1 d2.hdf5')
    cs = {'0.5': 'b',
          '1.0': 'g',
          '2.0': 'r'}
    def plot(ax, dpr, smear):
        p = p_template.format(smear, dpr)
        d = wt.kit.read_h5(p)
        xi = normalize_frequency(d['w1'])
        yi = normalize_frequency(d['w2'])
        zi = d['arr'][:, :, 10, 10].T
        zi = np.sqrt(zi)
        xi, yi, zi = zoom_arrs(xi, yi, zi)
        zi /= zi.max()
        levels = np.linspace(0, 1, 200)
        plt.contourf(xi, yi, zi, levels=levels, cmap=cmap)
        plt.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
        adjust_spines(ax, cs[dpr])
        ax.grid()
        wt.artists.diagonal_line(xi, yi)
        if ax.is_first_col():
            ax.set_ylabel(frequency_label(2), fontsize=18)
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        if ax.is_last_row():
            ax.set_xlabel(frequency_label(1), fontsize=18)
        else:
            plt.setp(ax.get_xticklabels(), visible=False)        
        ax.set_xlim(-4, 4)
        ax.set_ylim(-4, 4)   
        ticks = [-3, -2, -1, 0, 1, 2, 3]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        # ellipse
        a, b = get_widths(xi, yi, zi)
        e = Ellipse(xy=[0, 0], width=a, height=b, angle=45, color='y', lw=5, fill=False)
        e.set_alpha(1)
        e.set_zorder(500)
        ax.add_artist(e)
        ellipticity = ((a**2)-(b**2))/((a**2)+(b**2))
        label = r'${\mathcal{E}=' + '{:.2f}'.format(ellipticity) + r'}$'
        wt.artists.corner_text(label, corner='LR', fontsize=18)
    # fill out
    # col 1
    ax = plt.subplot(gs[0, 0])
    plot(ax, '0.5', '0.0')
    ax.set_title(r'$\mathsf{\Delta_{inhom} = 0.0}$', fontsize=18)
    ax = plt.subplot(gs[1, 0])
    plot(ax, '1.0', '0.0')
    ax = plt.subplot(gs[2, 0])
    plot(ax, '2.0', '0.0')
    # col 2
    ax = plt.subplot(gs[0, 1])
    plot(ax, '0.5', '0.5')
    ax.set_title(r'$\mathsf{\Delta_{inhom} = 0.5 \Delta_\omega}$', fontsize=18)
    ax = plt.subplot(gs[1, 1])
    plot(ax, '1.0', '0.5')
    ax = plt.subplot(gs[2, 1])
    plot(ax, '2.0', '0.5')
    # col 3
    ax = plt.subplot(gs[0, 2])
    plot(ax, '0.5', '1.0')
    ax.set_title(r'$\mathsf{\Delta_{inhom} = 1.0 \Delta_\omega}$', fontsize=18)
    ax = plt.subplot(gs[1, 2])
    plot(ax, '1.0', '1.0')
    ax = plt.subplot(gs[2, 2])
    plot(ax, '2.0', '1.0')
    # col 4
    ax = plt.subplot(gs[0, 3])
    plot(ax, '0.5', '2.0')
    ax.set_title(r'$\mathsf{\Delta_{inhom} = 2.0 \Delta_\omega}$', fontsize=18)
    ax.text(1.01, 0.5, r'$\mathsf{\Gamma_{10}\Delta_t = 2.0}$', rotation=-90, fontsize=18, va='center', transform=ax.transAxes)
    ax = plt.subplot(gs[1, 3])
    plot(ax, '1.0', '2.0')
    ax.text(1.01, 0.5, r'$\mathsf{\Gamma_{10}\Delta_t = 1.0}$', rotation=-90, fontsize=18, va='center', transform=ax.transAxes)
    ax = plt.subplot(gs[2, 3])
    plot(ax, '2.0', '2.0')
    ax.text(1.01, 0.5, r'$\mathsf{\Gamma_{10}\Delta_t = 0.5}$', rotation=-90, fontsize=18, va='center', transform=ax.transAxes)
    # colorbar
    cax = plt.subplot(gs[:, -1])
    wt.artists.plot_colorbar(cax=cax, cmap=cmap, label='amplitude', label_fontsize=18)
    # finish
    wt.artists.savefig(output_path, fig=fig, close=True)


### 2D frequencies at -4 ######################################################


output_path = os.path.join(directory, '2D frequences at -4.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # set up figure
    fig, gs = wt.artists.create_figure(width='double', nrows=3,
                                       cols=[1, 1, 1, 1, 0, 'cbar'])
    # fit method
    def get_widths(xi, yi, zi):
        '''
        defined on the amplitude level
        '''
        # go by fitting
        def gauss(p, x):
            A, mu, w = p
            z = (mu - x) / (np.sqrt(2) * w)
            out = A * np.exp(-z**2)
            return out
        def erf1(p, x, y):
            return y - gauss(p,x)
        outs = []
        # diagonal
        y = zi.diagonal()
        p0 = []
        p0.append(np.nanmax(y))  # amplitude
        p0.append(0.)  # center
        p0.append(1.)  # width
        out = leastsq(erf1, p0, args=(xi, y), full_output=False)[0]
        outs.append(2.35482*out[2])
        # antidiagonal
        y = zi[...,::-1].diagonal()
        p0 = []
        p0.append(np.nanmax(y))  # amplitude
        p0.append(0.)  # center
        p0.append(1.)  # width
        out = leastsq(erf1, p0, args=(xi, y), full_output=False)[0]
        outs.append(2.35482*out[2])
        return outs
    # worker method
    p_template = os.path.join(directory, 'measured', 'smear {}', 'dpr {} TO all w1 w2 d1 d2.hdf5')
    cs = {'0.5': 'b',
          '1.0': 'g',
          '2.0': 'r'}
    def plot(ax, dpr, smear):
        p = p_template.format(smear, dpr)
        d = wt.kit.read_h5(p)
        xi = normalize_frequency(d['w1'])
        yi = normalize_frequency(d['w2'])
        zi = d['arr'][:, :, 10, 20].T
        zi = np.sqrt(zi)
        xi, yi, zi = zoom_arrs(xi, yi, zi)
        zi /= zi.max()
        levels = np.linspace(0, 1, 200)
        plt.contourf(xi, yi, zi, levels=levels, cmap=cmap)
        plt.contour(xi, yi, zi, levels=contour_levels, colors='k', alpha=0.5)
        adjust_spines(ax, cs[dpr])
        ax.grid()
        wt.artists.diagonal_line(xi, yi)
        if ax.is_first_col():
            ax.set_ylabel(frequency_label(2), fontsize=18)
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        if ax.is_last_row():
            ax.set_xlabel(frequency_label(1), fontsize=18)
        else:
            plt.setp(ax.get_xticklabels(), visible=False)        
        ax.set_xlim(-4, 4)
        ax.set_ylim(-4, 4)   
        ticks = [-3, -2, -1, 0, 1, 2, 3]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        # ellipse
        a, b = get_widths(xi, yi, zi)
        e = Ellipse(xy=[0, 0], width=a, height=b, angle=45, color='y', lw=5, fill=False)
        e.set_alpha(1)
        e.set_zorder(500)
        ax.add_artist(e)
        ellipticity = ((a**2)-(b**2))/((a**2)+(b**2))
        label = r'${\mathcal{E}=' + '{:.2f}'.format(ellipticity) + r'}$'
        wt.artists.corner_text(label, corner='LR', fontsize=18)
    # fill out
    # col 1
    ax = plt.subplot(gs[0, 0])
    plot(ax, '0.5', '0.0')
    ax.set_title(r'$\mathsf{\Delta_{inhom} = 0.0}$', fontsize=18)
    ax = plt.subplot(gs[1, 0])
    plot(ax, '1.0', '0.0')
    ax = plt.subplot(gs[2, 0])
    plot(ax, '2.0', '0.0')
    # col 2
    ax = plt.subplot(gs[0, 1])
    plot(ax, '0.5', '0.5')
    ax.set_title(r'$\mathsf{\Delta_{inhom} = 0.5 \Delta_\omega}$', fontsize=18)
    ax = plt.subplot(gs[1, 1])
    plot(ax, '1.0', '0.5')
    ax = plt.subplot(gs[2, 1])
    plot(ax, '2.0', '0.5')
    # col 3
    ax = plt.subplot(gs[0, 2])
    plot(ax, '0.5', '1.0')
    ax.set_title(r'$\mathsf{\Delta_{inhom} = 1.0 \Delta_\omega}$', fontsize=18)
    ax = plt.subplot(gs[1, 2])
    plot(ax, '1.0', '1.0')
    ax = plt.subplot(gs[2, 2])
    plot(ax, '2.0', '1.0')
    # col 4
    ax = plt.subplot(gs[0, 3])
    plot(ax, '0.5', '2.0')
    ax.set_title(r'$\mathsf{\Delta_{inhom} = 2.0 \Delta_\omega}$', fontsize=18)
    ax.text(1.01, 0.5, r'$\mathsf{\Gamma_{10}\Delta_t = 2.0}$', rotation=-90, fontsize=18, va='center', transform=ax.transAxes)
    ax = plt.subplot(gs[1, 3])
    plot(ax, '1.0', '2.0')
    ax.text(1.01, 0.5, r'$\mathsf{\Gamma_{10}\Delta_t = 1.0}$', rotation=-90, fontsize=18, va='center', transform=ax.transAxes)
    ax = plt.subplot(gs[2, 3])
    plot(ax, '2.0', '2.0')
    ax.text(1.01, 0.5, r'$\mathsf{\Gamma_{10}\Delta_t = 0.5}$', rotation=-90, fontsize=18, va='center', transform=ax.transAxes)
    # colorbar
    cax = plt.subplot(gs[:, -1])
    wt.artists.plot_colorbar(cax=cax, cmap=cmap, label='amplitude', label_fontsize=18)
    # finish
    wt.artists.savefig(output_path, fig=fig, close=True)
