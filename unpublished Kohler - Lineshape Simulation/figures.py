### import ####################################################################


import os
import itertools

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
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
        return r'$\mathsf{\tau_{22^{\prime}}/w_t}$'
    elif kind == 2:
        return r'$\mathsf{\tau_{21}/w_t}$'
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
    wmel.add_arrow(0, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(1, [0, 1], 'bra', '2', color='b')
    wmel.add_arrow(2, [1, 0], 'bra', '2\'', color='b')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    for i in range(len(energies)):
        ax.text(-state_text_buffer, energies[i], state_names[i], fontsize=12, verticalalignment='center', horizontalalignment ='right')
    # pathway 1 beta
    ax = plt.subplot(gs[1, 2])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(1, [0, 1], 'bra', '2', color='b')
    wmel.add_arrow(2, [1, 2], 'bra', '2\'', color='b')
    wmel.add_arrow(3, [2, 1], 'out', '', color='r')
    for i in range(len(energies)):
        ax.text(-state_text_buffer, energies[i], state_names[i], fontsize=12, verticalalignment='center', horizontalalignment ='right')
    # pathway 1 gamma
    ax = plt.subplot(gs[2, 2])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(1, [1, 0], 'ket', '2', color='b')
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
    wmel.add_arrow(2, [2, 1], 'ket', '2', color='b')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    # pathway 2 beta
    ax = plt.subplot(gs[1, 3])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(1, [1, 2], 'ket', '2\'', color='b')
    wmel.add_arrow(2, [0, 1], 'bra', '2', color='b')
    wmel.add_arrow(3, [2, 1], 'out', '', color='r')
    # pathway 3 ---------------------------------------------------------------
    # pathway 3 alpha
    ax = plt.subplot(gs[0, 4])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='III')
    wmel.add_arrow(0, [0, 1], 'bra', '2', color='b')
    wmel.add_arrow(1, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(2, [1, 0], 'bra', '2\'', color='b')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    # pathway 3 beta
    ax = plt.subplot(gs[1, 4])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'bra', '2', color='b')
    wmel.add_arrow(1, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(2, [1, 2], 'ket', '2\'', color='b')
    wmel.add_arrow(3, [2, 1], 'out', '', color='r')
    # pathway 3 gamma
    ax = plt.subplot(gs[2, 4])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'bra', '2', color='b')
    wmel.add_arrow(1, [1, 0], 'bra', '1', color='r')
    wmel.add_arrow(2, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    # pathway 4 ---------------------------------------------------------------
    # pathway 4 alpha
    ax = plt.subplot(gs[0, 5])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='IV')
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(1, [1, 2], 'ket', '1', color='r')
    wmel.add_arrow(2, [2, 1], 'ket', '2', color='b')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    # pathway 4 beta
    ax = plt.subplot(gs[1, 5])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(1, [1, 2], 'ket', '1', color='r')
    wmel.add_arrow(2, [0, 1], 'bra', '2', color='b')
    wmel.add_arrow(3, [2, 1], 'out', '', color='r')
    # pathway 5 ---------------------------------------------------------------
    # pathway 5 alpha
    ax = plt.subplot(gs[0, 6])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='V')
    wmel.add_arrow(0, [0, 1], 'bra', '2', color='b')
    wmel.add_arrow(1, [1, 0], 'bra', '2\'', color='b')
    wmel.add_arrow(2, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    # pathway 5 beta
    ax = plt.subplot(gs[1, 6])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'bra', '2', color='b')
    wmel.add_arrow(1, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(2, [1, 2], 'ket', '1', color='r')
    wmel.add_arrow(3, [2, 1], 'out', '', color='r')
    # pathway 5 gamma
    ax = plt.subplot(gs[2, 6])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'bra', '2', color='b')
    wmel.add_arrow(1, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(2, [1, 0], 'bra', '1', color='r')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    # pathway 6 ---------------------------------------------------------------
    text_buffer = 1.3
    # pathway 6 alpha
    ax = plt.subplot(gs[0, 7])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='VI')
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(1, [1, 0], 'ket', '2', color='b')
    wmel.add_arrow(2, [0, 1], 'ket', '1', color='r')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    ax.text(text_buffer, 0.5, r'$\mathsf{\alpha}$', fontsize=25, verticalalignment='center', horizontalalignment='left')
    # pathway 6 beta
    ax = plt.subplot(gs[1, 7])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(1, [0, 1], 'bra', '2', color='b')
    wmel.add_arrow(2, [1, 2], 'ket', '1', color='r')
    wmel.add_arrow(3, [2, 1], 'out', '', color='r')
    ax.text(text_buffer, 0.5, r'$\mathsf{\beta}$', fontsize=25, verticalalignment='center', horizontalalignment='left')
    # pathway 6 gamma
    ax = plt.subplot(gs[2, 7])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(1, [0, 1], 'bra', '2', color='b')
    wmel.add_arrow(2, [1, 0], 'bra', '1', color='r')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
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
    fig, gs = wt.artists.create_figure(width='double', nrows=2, cols=[1, 1.5, 1, 'cbar'],
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
    # evolution of density matrix terms in pw5 --------------------------------
    ax = plt.subplot(gs[0, 1])
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
    ax.plot(xi, yi, lw=2, c='g', label=r'$\mathsf{\rho_{ga}}$')
    # aa
    yi = d['aa']
    ax.plot(xi, yi, lw=2, c='m', label=r'$\mathsf{\rho_{aa}}$')
    # ag2
    yi = d['ag']
    ax.plot(xi, yi, lw=2, c='c', label=r'$\mathsf{\rho_{ag}}$')
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
    # FT of FID above and mono pass function ----------------------------------
    ax = plt.subplot(gs[0, 2:])
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
    ax.set_xlabel(frequency_label('out'), fontsize=18)
    # measured 2D frequency with homo lineshape -------------------------------
    ax = plt.subplot(gs[1, 0])
    p = os.path.join(directory, 'measured', 'dpr 1.0 TO all w1 w2 d1 d2.hdf5')
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
    wt.artists.corner_text('d', ax=ax, fontsize=18, background_alpha=1)
    ax.set_xlabel(frequency_label('1'), fontsize=18)
    ax.set_ylabel(frequency_label('2'), fontsize=18)
    # representation of kernal ------------------------------------------------
    ax = plt.subplot(gs[1, 1], projection='3d')
    ax.view_init(35, 75 )
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    plt.setp(ax.get_zticklabels(), visible=False)
    xi = np.linspace(-2, 2, 1000)
    yi = xi
    function = wt.fit.Gaussian()
    zi = function.evaluate([0, 1/3., 1, 0], xi)
    verts = [list(zip(xi, yi, zi))]
    poly = Poly3DCollection(verts, facecolors=[[0., 0., 0., 0.1]])
    ax.add_collection3d(poly)
    ax.set_xlabel(r'$\mathsf{\omega_1}$', fontsize=18)
    ax.set_xlim3d(-2, 2)
    ax.yaxis.set_rotate_label(False)  # disable automatic rotation
    ax.set_ylabel(r'$\mathsf{\omega_2}$', fontsize=18, rotation=115)
    ax.set_ylim3d(-2, 2)
    ax.set_zlim3d(0, 1)
    wt.artists.corner_text('e', ax=ax, fontsize=18, background_alpha=1)
    # inhomogenious -----------------------------------------------------------
    ax = plt.subplot(gs[1, 2])
    p = os.path.join(directory, 'measured with inhomogeneity', 'smear 2.0', 
                     'dpr 1.0 TO all w1 w2 d1 d2.hdf5')
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
    wt.artists.corner_text('f', ax=ax, fontsize=18, background_alpha=1)
    # show inhomo linewidth applied
    # TODO: ensure that this matches the actual linewidth (FWHM) used...
    center = 7000.
    width = 2*delta_omega
    l = np.linspace(center-(width/2), center+(width/2), 100)
    l = normalize_frequency(l)
    ax.plot(l, l, c='k', lw=5)
    ax.set_xlabel(frequency_label('1'), fontsize=18)
    ax.set_ylabel(frequency_label('2'), fontsize=18)
    # colorbar ----------------------------------------------------------------
    ticks = np.linspace(0, 1, 11)
    cax = plt.subplot(gs[1, -1])
    matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, ticks=ticks)
    cax.set_ylabel('amplitude', fontsize=18)
    # finish
    plt.savefig(output_path, dpi=300, transparent=True, pad_inches=1)
    plt.close('fig')
        

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
    # 1
    ax = plt.subplot(gs[1, 0])
    add_data(ax, 1, label='2.0')
    # 2
    ax = plt.subplot(gs[2, 0])
    add_data(ax, 2, label='1.0')
    ax.set_ylabel('amplitude', fontsize=18, labelpad=20)
    # 3
    ax = plt.subplot(gs[3, 0])
    add_data(ax, 3, label='0.5')
    # 4
    ax = plt.subplot(gs[4, 0])
    add_data(ax, 4, label='0.0', show_x=True)
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
# TODO: consider adding diagrams showing pulse envelopes as in simulation overview


output_path = os.path.join(directory, 'pw1 lineshapes.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # preapre figure
    fig, gs = wt.artists.create_figure(width='double', nrows=2, 
                                       cols=[1, 1, 0.25, 1, 1, 'cbar'])
    # data
    title = None        
    filepath = os.path.join(directory, 'measured', 'dpr 1.0 TO all w1 w2 d1 d2.hdf5')
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
    # preapre for frequency plots
    xi = normalize_frequency(d['w1'])
    yi = normalize_frequency(d['w2'])
    xi = scipy.ndimage.interpolation.zoom(xi, 5)
    yi = scipy.ndimage.interpolation.zoom(yi, 5)
    def adjust_spines(ax, color):
        for key in ['bottom', 'top', 'right', 'left']:
            ax.spines[key].set_color(color)
            ax.spines[key].set_linewidth(5)
    freq_ticks = [-2, 0, 2]
    # UL
    ax = plt.subplot(gs[0, 3])
    adjust_spines(ax, 'cyan')
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
    ax = plt.subplot(gs[0, 4])
    adjust_spines(ax, 'm')
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
    ax = plt.subplot(gs[1, 3])
    adjust_spines(ax, 'pink')
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
    ax = plt.subplot(gs[1, 4])
    adjust_spines(ax, 'orange')
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
    def plot(gs_col, dpr, manuals, yticks, title=''):
        sps = gs[1, gs_col]
        ax = plt.subplot(sps)
        if not yticks:
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            ax.set_ylabel(r'$\mathsf{\tau_{21}/w_t}$', fontsize=18)
        cmap = wt.artists.colormaps['default']
        # get data from zip
        template_fname = os.path.join('measured', 'dpr {0} TO {1} w1 w2 d1 d2.hdf5')
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
        ax.set_xlabel(r'$\mathsf{\tau_{22^{\prime}}/w_t}$', fontsize=18)
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
        def plot_slice(TO, c=c):
            filepath = template_fname.format(dpr, TO)
            npz = wt.kit.read_h5(filepath)
            arr = np.sqrt(npz['arr'][20, 20].T)
            arr /= amps_max
            ax.plot(xi, arr[10], lw=2, c=c)
        plot_slice(5, 'r')
        plot_slice(6, 'b')
        plot_slice(1, 'g')
        plot_slice(3, 'g')
        ax.set_title(title, fontsize=18, y=1.08)
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
    plot(0, dprs[0], manuals=manuals[0], yticks=True, title=r'$\mathsf{\Gamma_{10}\Delta_t = 2.0}$')
    plot(1, dprs[1], manuals=manuals[1], yticks=False, title=r'$\mathsf{\Gamma_{10}\Delta_t = 1.0}$')
    plot(2, dprs[2], manuals=manuals[2], yticks=False, title=r'$\mathsf{\Gamma_{10}\Delta_t = 0.5}$')
    # colorbar
    ticks = np.linspace(0, 1, 11)
    cax = plt.subplot(gs[:, -1])
    matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, ticks=ticks)
    cax.set_ylabel('amplitude', fontsize=18)
    # finish
    plt.savefig(output_path, dpi=300, transparent=True, pad_inches=1)
    plt.close('fig')


### wigners ###################################################################


output_path = os.path.join(directory, 'wigners.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    cmap = wt.artists.colormaps['default']
    def plot(ax, p, w2_position):
        d = wt.kit.read_h5(p)
        # TODO (?) put version of this for other dprs in SI
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
    fig, gs = wt.artists.create_figure(nrows=1, cols=[1, 1, 1, 1, 1, 'cbar'], width='double')
    p = os.path.join('measured', 'dpr 1.0 TO all w1 w2 d1 d2.hdf5')
    # red
    ax = plt.subplot(gs[0, 0])
    plot(ax, p, 'rr')
    ax.set_xlabel(frequency_label(1), fontsize=18)
    ax.set_ylabel(time_label(), fontsize=18)
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


### test ######################################################################


output_path = os.path.join(directory, 'test.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:    
    # create figure
    cols = ['cbar', 0.5, 1, 0.5, 1, 'cbar']
    aspects = [[[0, 1], 0.3],
               [[1, 1], 1.0],
               [[2, 1], 0.25],
               [[3, 1], 0.3],
               [[4, 1], 1.0]]
    fig, gs = wt.artists.create_figure(cols=cols, nrows=5, width='double',
                                       aspects=aspects, hspace=0.5)
    # time traces -------------------------------------------------------------
    def plot(ax, d1, d2):
        # TODO: precalculate the things we need here...
        import NISE
        from NISE.lib import pulse
        from NISE.lib.misc.__init__ import NISE_path
        import numpy as np
        # important:  we can't call on scan directly for some reason; use trive to do it
        import NISE.lib.measure as m
        import NISE.experiments.trive as trive
        import NISE.hamiltonians.H0 as H0
        import NISE.hamiltonians.params.inhom as inhom
        trive.exp.set_coord(trive.d2, -d2)
        trive.exp.set_coord(trive.ws, 7000.)  # TODO: change
        trive.exp.set_coord(trive.ss, 50.)
        trive.exp.timestep = 0.1
        trive.exp.early_buffer = 500
        trive.exp.late_buffer = 500
        # hamiltonian
        H0.tau_ag = 50.
        H0.tau_2aa = 50.
        H0.tau_2ag = 50.
        H = H0.Omega()
        H.out_group = [[i] for i in range(7)]
        H.TOs = [1]  # other pathways complicate interpretation...
        # scan
        axis_d1 = trive.d1  # I seem to need an axis argument...
        axis_d1.points = np.array([-d1])  # this is where you set D1, really
        scan = trive.exp.scan(trive.d1, H=H)
        scan.run(mp=False, autosave=False)
        # plot pulses
        arr = scan.efields(windowed=True)
        xi = np.linspace(-scan.early_buffer, scan.late_buffer, arr.shape[-1])
        xi /= 50.
        # -2
        yi = np.abs(arr[0, 1, :])**2
        yi /= yi.max()
        ax.fill_between(xi, 0, yi, facecolor='b', alpha=0.25)
        # 2'
        yi = np.abs(arr[0, 2, :])**2
        yi /= yi.max()
        ax.fill_between(xi, 0, yi, facecolor='b', alpha=0.25)
        # 1
        yi = np.abs(arr[0, 0, :])**2
        yi /= yi.max()
        ax.fill_between(xi, 0, yi, facecolor='r', alpha=0.25)
        # plot density matrix terms
        # ag
        yi = np.abs(scan.sig[0, 1, :])**2
        yi /= yi.max()
        ax.plot(xi, yi, lw=2, c='g', label=r'$\mathsf{\rho_{ga}}$')
        # gg1 (actually aa)
        yi = np.abs(scan.sig[0, 3, :])**2
        yi /= yi.max()
        ax.plot(xi, yi, lw=2, c='m', label=r'$\mathsf{\rho_{aa}}$')
        # ag2
        yi = np.abs(scan.sig[0, 5, :])**2
        yi /= yi.max()
        ax.plot(xi, yi, lw=2, c='c', label=r'$\mathsf{\rho_{ag2}}$')
        ax.set_ylim(0, 1.15)
        ax.set_xlim(-7, 7)
        plt.grid()
        ax.set_xlabel(time_label(), fontsize=18)
        plt.setp(ax.get_yticklabels(), visible=False)
    # UL
    ax = plt.subplot(gs[0, 2])
    scan = plot(ax, -100, 100)
    # UR
    ax = plt.subplot(gs[0, 4])
    # LL
    ax = plt.subplot(gs[3, 2])
    # LR
    ax = plt.subplot(gs[3, 4])
    # 2D frequencies ----------------------------------------------------------
    filepath = os.path.join(directory, 'measured', 'dpr 1.0 TO all w1 w2 d1 d2.hdf5')
    d = wt.kit.read_h5(filepath)
    arr = np.sqrt(d['arr'])
    def plot(ax, zi):
        xi = normalize_frequency(d['w1'].copy())
        yi = normalize_frequency(d['w2'].copy())
        zi /= zi.max()
        xi, yi, zi = zoom_arrs(xi, yi, zi)
        levels = np.linspace(0, 1, 200)
        ax.contourf(xi, yi, zi, cmap=cmap, levels=levels)
        contour_levels = [0.25, 0.5, 0.75]
        ax.contour(xi, yi, zi, contour_levels, colors='k', alpha=0.5)
        ax.set_xlim(-4, 4)
        ax.set_ylim(-4, 4)
        ticks = [-2, 0, 2]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.grid()
        wt.artists.diagonal_line(xi, yi)
        ax.set_xlabel(frequency_label('1'))
        ax.set_ylabel(frequency_label('2'))
    # UL
    ax = plt.subplot(gs[1, 2])
    zi = arr[:, :, 15, 5].T
    plot(ax, zi)
    plt.setp(ax.get_xticklabels(), visible=False)
    # UR
    ax = plt.subplot(gs[1, 4])
    zi = arr[:, :, 10, 5].T
    plot(ax, zi)
    ax.set_ylabel(frequency_label('2'))
    plt.setp(ax.get_xticklabels(), visible=False)
    # LL
    ax = plt.subplot(gs[4, 2])
    zi = arr[:, :, 15, 10].T
    plot(ax, zi)
    ax.set_ylabel(frequency_label('2'))
    plt.setp(ax.get_xticklabels(), visible=False)
    # LR
    ax = plt.subplot(gs[4, 4])
    zi = arr[:, :, 10, 10].T
    plot(ax, zi)
    ax.set_ylabel(frequency_label('2'))
    plt.setp(ax.get_xticklabels(), visible=False)
    # left colorbar
    cax = plt.subplot(gs[:, 0])
    wt.artists.plot_colorbar(cax=cax, cmap=coolwarm, ticklocation='left')
    cax = plt.subplot(gs[:, -1])
    wt.artists.plot_colorbar(cax=cax, cmap='default')
    # finish
    wt.artists.savefig(output_path)
