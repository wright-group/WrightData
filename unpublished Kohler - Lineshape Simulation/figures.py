### import ####################################################################


import os
import itertools

import matplotlib
matplotlib.rcParams['font.size'] = 14
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
plt.close('all')

import scipy
from scipy import constants
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

cmap = wt.artists.colormaps['default']

props = dict(boxstyle='square', facecolor='white', alpha=0.8)

template_fname = r'npzs\dpr {0} TO {1} w1 w2 d1 d2 arr.npz'
dprs = [0.5, 1., 2.0]


def zoom_arrs(*args):
    return [scipy.ndimage.interpolation.zoom(arr, 5) for arr in args]


def normalize_frequency(arr):
    arr -= arr.mean()
    arr /= 2 * np.log(2) / (50*np.pi * 3e10) * 1e15
    return arr


def normalize_delay(arr):
    arr /= 50.
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
    return r'$\mathsf{t/w_t}$'


def frequency_label(sub=''):
    return r'$\mathsf{\nu_{' + sub + r'}/w_{\nu}}$'


freq_ticks = [-1, 0, 1]
delay_ticks = [-2, 0, 2]

### WMELs #####################################################################


output_path = os.path.join(directory, 'WMELs.png')

force_plotting = True

if not os.path.isfile(output_path) or force_plotting:
    # prepare for plot
    fig = plt.figure(figsize=[6.5, 6])
    gs = grd.GridSpec(3, 8, width_ratios=[1.5, 1, 1, 1, 1, 1, 1, 1])
    # large levels figure
    ax = plt.subplot(gs[:, 0])
    # plot energies
    energies = list(np.linspace(0, 0.85, 4)) + [1]
    print energies
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
    add_arrow([3, 2], 0.75, 8)
    add_arrow([2, 1], 0.75, 4)
    add_arrow([1, 0], 0.75, 2)
    state_text_buffer = 0.25
    state_names = ['n=0', 'n=1', 'n=2', 'n=3', 'n=N']
    for i in range(len(energies)):
        ax.text(-state_text_buffer, energies[i], state_names[i], fontsize=18, verticalalignment='center', horizontalalignment ='right')
    ax.text(-state_text_buffer, 0.925, '...', rotation=-90, fontsize=18, verticalalignment='center', horizontalalignment ='right')
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.025, 1.025)
    ax.axis('off')
    # highlights
    import matplotlib.patches as patches
    patch1 = patches.Rectangle((0.3525, 0.1775), 0.245, 0.175, transform=fig.transFigure, figure=fig, facecolor='y', alpha=0.25)
    patch2 = patches.Rectangle((0.689, 0.6476), 0.16, 0.175, transform=fig.transFigure, figure=fig, facecolor='y', alpha=0.25)
    fig.patches.extend([patch1, patch2])
    # pathway 1 ---------------------------------------------------------------
    energies = [0., 0.5, 1.]
    state_text_buffer = 0.25
    state_names = ['g', 'a', '2a']
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
    ax.text(text_buffer, 0.5, r'$\mathrm{\alpha}$', fontsize=25, verticalalignment='center', horizontalalignment='left')
    # pathway 6 beta
    ax = plt.subplot(gs[1, 7])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(1, [0, 1], 'bra', '2', color='b')
    wmel.add_arrow(2, [1, 2], 'ket', '1', color='r')
    wmel.add_arrow(3, [2, 1], 'out', '', color='r')
    ax.text(text_buffer, 0.5, r'$\mathrm{\beta}$', fontsize=25, verticalalignment='center', horizontalalignment='left')
    # pathway 6 gamma
    ax = plt.subplot(gs[2, 7])
    wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
    wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
    wmel.add_arrow(1, [0, 1], 'bra', '2', color='b')
    wmel.add_arrow(2, [1, 0], 'bra', '1', color='r')
    wmel.add_arrow(3, [1, 0], 'out', '', color='r')
    ax.text(text_buffer, 0.5, r'$\mathrm{\gamma}$', fontsize=25, verticalalignment='center', horizontalalignment='left')
    # labels
    dist = 0.05
    fig.text(0.075, 1-dist, 'a', bbox=props, horizontalalignment='left', verticalalignment='top', fontsize=18, alpha=1)
    fig.text(0.33, 1-dist, 'b', bbox=props, horizontalalignment='left', verticalalignment='top', fontsize=18, alpha=1)    
    # finish
    wt.artists.subplots_adjust(fig)
    plt.savefig(output_path, dpi=300, transparent=True,  pad_inches=1.)
    plt.close(fig)


### simulation overview #######################################################


output_path = os.path.join(directory, 'simulation_overview.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # prepare for plot
    fig, gs = wt.artists.create_figure(width='double', nrows=2, cols=[1, 1.5, 1, 'cbar'],
                                       hspace=0.75)
    # set of coherence times vs lab time --------------------------------------
    # calculating transients directly right here is easiest
    def rho(t,tau,mu,E, last_rho):
        f = np.zeros(t.shape, dtype=np.complex64)
        for i in range(t.size-1):
            dt = t[i+1]-t[i]
            df1 = 0.5 * E[i] * last_rho[i] - 1/tau * f[i]
            ftemp = f[i] + df1 * dt
            df2 = 0.5 * E[i+1] * last_rho[i+1] - 1/tau * ftemp
            f[i+1] = f[i] + 0.5 * dt * ( df1 + df2 )
        return f
    ax = plt.subplot(gs[0, 0])
    start_fine_int = -1000
    end_fine_int = 1000
    t1 = np.linspace(-10000, start_fine_int, num=10000)
    t2 = np.linspace(start_fine_int+1, end_fine_int-1, num=20000)
    t3 = np.linspace(end_fine_int,20000, num=10000)
    t = np.hstack([t1,t2,t3])
    # pulse
    function = wt.fit.Gaussian()
    yi = normalized_gauss(t, 50.)
    yi /= yi.max()
    xi = t / 50
    ax.fill_between(xi, 0, yi, facecolor='grey', alpha=0.25)
    # rhos
    taus = [25., 50., 100.]  # fs
    cs = ['b', 'g', 'r']
    labels = [r'$\mathsf{w_t/2}$',
              r'$\mathsf{w_t}$',
              r'$\mathsf{2w_t}$']
    for i, tau in enumerate(taus):
        mu = 1.
        FWHM = 50.  # fs
        xi = t / 50
        yi = rho(t,tau,mu,normalized_gauss(t, FWHM), np.ones(len(t)))
        yi = np.abs(yi)
        yi /= yi.max()
        ax.plot(xi, yi, c=cs[i], lw=2, label=labels[i])
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
    trive.exp.set_coord(trive.d2, 200.)
    trive.exp.set_coord(trive.ws, 7000.)
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
    H.TOs = [5]  # other pathways complicate interpretation...
    # scan
    d1 = trive.d1  # I seem to need an axis argument...
    d1.points = np.array([100])  # this is where you set D1, really
    scan = trive.exp.scan(trive.d1, H=H)
    scan.run(mp=False, autosave=False)
    # plot pulses
    arr = scan.efields(windowed=True)
    # TODO: decide on time convention (positive, negative)
    xi = np.linspace(-scan.early_buffer, scan.late_buffer, arr.shape[-1])
    xi /= 50.
    # -2
    yi = np.abs(arr[0, 1, :])**2
    yi /= yi.max()
    ax.fill_between(xi, 0, yi, facecolor='b', alpha=0.25)
    ax.text(-4, 0.1, '2', fontsize=18, ha='center', va='center')
    ax.text(-4, 0.3, 'x', fontsize=18, ha='center', va='center')
    # 2'
    yi = np.abs(arr[0, 2, :])**2
    yi /= yi.max()
    ax.fill_between(xi, 0, yi, facecolor='b', alpha=0.25)
    ax.text(-2, 0.1, '2\'', fontsize=18, ha='center', va='center')
    ax.text(-2, 0.3, 'y', fontsize=18, ha='center', va='center')
    # 1
    yi = np.abs(arr[0, 0, :])**2
    yi /= yi.max()
    ax.fill_between(xi, 0, yi, facecolor='r', alpha=0.25)
    ax.text(0, 0.1, '1', fontsize=18, ha='center', va='center')
    ax.text(0   , 0.3, 'z', fontsize=18, ha='center', va='center')
    # plot density matrix terms
    # ga
    yi = np.abs(scan.sig[0, 2, :])**2
    yi /= yi.max()
    ax.plot(xi, yi, lw=2, c='g', label=r'$\mathsf{\rho_{ga}}$')
    # aa
    yi = np.abs(scan.sig[0, 3, :])**2
    yi /= yi.max()
    ax.plot(xi, yi, lw=2, c='m', label=r'$\mathsf{\rho_{aa}}$')
    # ag2
    yi = np.abs(scan.sig[0, 5, :])**2
    yi /= yi.max()
    ax.plot(xi, yi, lw=2, c='c', label=r'$\mathsf{\rho_{ag2}}$')
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
    # TODO: make this correspond to actual calculations
    # FID
    FID = scan.sig[0, 5, :]
    t = np.linspace(-scan.early_buffer, scan.late_buffer, FID.size)
    # amplitudes
    wn_to_omega = NISE.lib.misc.wn_to_omega
    amps = FID * np.exp(1j*scan.get_color()[0]*wn_to_omega*t)
    amps = np.fft.fftshift(np.fft.fft(amps))
    amps = np.abs(amps)
    yi = amps
    yi /= yi.max()
    # hz
    d = (t.max()-t.min())/(t.size-1)
    d *= 1e-15  # seconds
    hz = np.fft.fftfreq(FID.size, d=d)  # Hz
    hz = np.fft.fftshift(hz)  # centered around zero Hz
    # wn
    xi = hz/(constants.c*1e2)
    xi /= 2 * np.log(2) / (50*np.pi * 3e10) * 1e15
    xi = scipy.ndimage.interpolation.zoom(xi, 10)
    yi = scipy.ndimage.interpolation.zoom(yi, 10)
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
    data = wt.data.from_pickle(os.path.join(directory, 'simulation_overview.p'))
    xi = normalize_frequency(data.w0.points)
    yi = normalize_frequency(data.w1.points)
    zi = data.simulation.values
    cmap = wt.artists.colormaps['default']
    levels= np.linspace(0, 1, 200)
    mappable = ax.contourf(xi, yi, zi, cmap=cmap, levels=levels)
    ax.contour(xi, yi, zi, 5, colors='k')
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
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.colors import colorConverter
    ax = plt.subplot(gs[1, 1], projection='3d')
    ax.view_init(35, 75 )
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    plt.setp(ax.get_zticklabels(), visible=False)
    xi = np.linspace(-2, 2, 1000)
    yi = xi
    function = wt.fit.Gaussian()
    zi = function.evaluate([0, 1/3., 1, 0], xi)
    verts = [zip(xi, yi, zi)]
    poly = Poly3DCollection(verts, facecolors=[colorConverter.to_rgba('k', alpha=0.1)])
    ax.add_collection3d(poly)
    ax.set_xlabel(r'$\mathsf{\nu_1}$', fontsize=18)
    ax.set_xlim3d(-2, 2)
    ax.yaxis.set_rotate_label(False)  # disable automatic rotation
    ax.set_ylabel(r'$\mathsf{\nu_2}$', fontsize=18, rotation=115)
    ax.set_ylim3d(-2, 2)
    ax.set_zlim3d(0, 1)
    wt.artists.corner_text('e', ax=ax, fontsize=18, background_alpha=1)
    # inhomogenious -----------------------------------------------------------
    ax = plt.subplot(gs[1, 2])
    data = wt.data.from_pickle(os.path.join(directory, 'simulation_overview_smeared.p'))
    xi = normalize_frequency(data.w0.points)
    yi = normalize_frequency(data.w1.points)
    zi = data.simulation.values
    cmap = wt.artists.colormaps['default']
    levels= np.linspace(0, 1, 200)
    mappable = ax.contourf(xi, yi, zi, cmap=cmap, levels=levels)
    ax.contour(xi, yi, zi, 5, colors='k')
    ax.grid()
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_xticks(freq_ticks)
    ax.set_yticks(freq_ticks)
    wt.artists.diagonal_line(xi, yi)
    wt.artists.corner_text('f', ax=ax, fontsize=18, background_alpha=1)
    # show inhomo linewidth applied
    center = 7000.
    width = 500.
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
    if False:  # messing around  
        plt.suptitle('NISE overview', fontsize=20)
        wt.artists.plot_margins(fig)
    # finish
    plt.savefig(output_path, dpi=300, transparent=True, pad_inches=1)
    plt.close('fig')


### hamiltonian ###############################################################


output_path = os.path.join(directory, 'hamiltonian.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    pass
    # TODO:
        

### FID v driven ##############################################################


output_path = os.path.join(directory, 'fid_v_driven.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # prepare figure
    aspects = [[[0, 0], 0.3],
               [[1, 0], 0.3],
               [[2, 0], 0.3],
               [[3, 0], 0.3],
               [[4, 0], 0.3]]
    fig, gs = wt.artists.create_figure(nrows=5, cols=[1, 'cbar'], aspects=aspects)
    # method to add data to axes (calculates on the fly)
    def add_data(ax, dpr, label, show_x=False):
        pulse_width = 50.
        slitwidth = 120.
        w_laser = 7050.
        # get axes from experiment
        w1 = trive.w1
        w2 = trive.w2
        ws = trive.ws
        d2 = trive.d2
        d1 = trive.d1
        # define points
        d2_points = np.zeros((1))
        d1_points = np.zeros((1))
        d1.points = d1_points
        d2.points = d2_points
        # set coords
        trive.exp.set_coord(w1, w_laser)
        trive.exp.set_coord(w2, w_laser)
        trive.exp.set_coord(trive.ss, pulse_width)
        trive.exp.timestep = 0.5
        # definitions external to the loop
        inhom_object = inhom.Inhom()
        m.Mono.slitwidth = slitwidth
        H1 = H0.Omega(tau_ag  = pulse_width * dpr,
                      tau_2aa = pulse_width * dpr,
                      tau_2ag = pulse_width * dpr,
                      TOs=[6])
        H1.out_group=[[1],[5],[6]]
        H1.wa_central = 7000.
        # exciton-exciton coupling
        H1.a_coupling = 75. # cm-1
        H1.mu_ag =  1.0
        H1.mu_2aa = H1.mu_ag # HO approx (1.414) vs. uncorr. electron approx. (1.)
        # change late buffer
        trive.exp.early_buffer = 250
        trive.exp.late_buffer = 500 #H0.Omega.tau_ag * 3 * max(dpr,1.)
        trive.exp.get_coords()
        # run
        out = trive.exp.scan(d1, d2, H=H1, inhom_object=inhom_object)
        out.run(autosave=False, mp=False)
        # process
        tprime = np.arange(-out.early_buffer, out.late_buffer, trive.exp.timestep)
        # convert the coherence to a rotating wave for easy phase analysis
        y = out.sig[0,0,0]
        y1 = np.exp(1j*H1.wa_central*NISE.lib.misc.wn_to_omega*tprime)*y
        y1 /= np.abs(y1).max()
        # try to convert signal to phase
        y2 = np.arctan(np.imag(y1)/np.real(y1))
        y3 = np.diff(y2)
        # filter for continuity
        for i in range(y3.size):
            if y3[i] > 3:
                print i, '>'
                y2[i+1:] -= np.pi
            elif y3[i] < -3:
                print i, '<'
                y2[i+1:] += np.pi
        # redefine the derivitive on y2 with disconitnuities removed
        y3 = np.diff(y2) 
        y2 /= NISE.lib.misc.wn_to_omega
        y3 /= out.timestep * NISE.lib.misc.wn_to_omega
        #y3[y3>1] = 0.
        y3 -= H1.wa_central
        y3 *= -1
        # plot ----------------------------------------------------------------
        ax.axvline(0, color='k', lw=2, zorder=10)
        # line
        ax.plot(tprime/pulse_width, np.abs(y1)/np.abs(y1).max(),
                linewidth=2., color='k', zorder=10)
        # colored inside
        color = (y3 - H1.wa_central) / (w_laser - H1.wa_central)
        import matplotlib.colors as colors
        import matplotlib.cm as cmx
        coolwarm = cm = plt.get_cmap('coolwarm') 
        cNorm  = colors.Normalize(vmin=0, vmax=1.)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
        for i in range(len(tprime[1:])):
            # calculate the color to use, given the instantaneous frequency
            colorVal = scalarMap.to_rgba(color[i])
            ax.plot([tprime[i+1]/pulse_width,tprime[i+1]/pulse_width], [0,np.abs(y1)[i]],
                     color=colorVal)
        # excitation pulse
        t = tprime/pulse_width
        yi = normalized_gauss(t, 1)
        yi /= yi.max()
        ax.plot(t, yi, color='grey', lw=2, zorder=15)
        # finish
        if show_x:
            ax.set_xlabel(r'$\mathsf{t / w_t}$', fontsize=18)
        else:
            plt.setp(ax.get_xticklabels(), visible=False)
        ax.set_ylabel(r'$\mathsf{|\rho_{10}(t)|}$', fontsize=18)
        ax.set_xlim(-1.5, 2.25)
        ax.set_ylim(0., 1.1)
        ax.set_xticks([-1, 0, 1, 2])
        ax.set_yticks([0,0.5,1.])
        plt.setp(ax.get_yticklabels(), visible=False)
        ax.grid()
        wt.artists.corner_text(label, fontsize=18, background_alpha=1.)
    # fill out axes
    add_data(plt.subplot(gs[0, 0]), 1e-2, r'$\mathsf{0}$')
    add_data(plt.subplot(gs[1, 0]), 0.5, r'$\mathsf{w_t/2}$')
    add_data(plt.subplot(gs[2, 0]), 1.0, r'$\mathsf{w_t}$')
    add_data(plt.subplot(gs[3, 0]), 2.0, r'$\mathsf{2w_t}$')
    add_data(plt.subplot(gs[4, 0]), 1e6, r'$\mathsf{\infty}$', show_x=True)
    # colorbar    
    from matplotlib import colorbar
    cax = plt.subplot(gs[:, -1])
    cb1 = colorbar.ColorbarBase(cax, cmap='coolwarm', ticks=[0,1.])
    plt.ylabel('$\mathsf{d\phi / dt (cm^{-1})}$', fontsize=18)
    plt.yticks(np.array([0,1.]), 
               (r'$\mathsf{\omega_{10}}$',
                r'$\mathsf{\omega_{laser}}$'), fontsize=18)
    # ????
    if False:
        plt.subplot(212)    
        plt.plot(tprime[1:]/pulse_width, y3,
                 linewidth=2)
        plt.plot(tprime/pulse_width,np.ones(tprime.shape)*H1.wa_central,
                 linewidth=2, linestyle='--', color='k')
        plt.plot(tprime/pulse_width,np.ones(tprime.shape)*w_laser,
                 linewidth=2, linestyle='--', color='k')
        plt.ylabel(r'$\mathsf{\frac{d\phi}{dt} (cm^{-1})}$', fontsize=14)
        plt.ylim([min(H1.wa_central,w_laser)-10.,max(H1.wa_central,w_laser)+10.])
        plt.xlabel(r'$\mathsf{t / \tau_{pulse}}$', fontsize=14)
        plt.grid()
    # finish
    plt.savefig(output_path, dpi=300, transparent=True, pad_inches=1)
    plt.close('fig')


### pathway I freq-freqs ######################################################


output_path = os.path.join(directory, 'pw1_2D_freqs.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    # preapre figure
    fig, gs = wt.artists.create_figure(width='double', nrows=2, 
                                       cols=[1, 1, 0.25, 1, 1, 'cbar'])
    # data
    title = None #r'$\mathsf{\tau=w_t/2}$'          
    filepath = template_fname.format('0.5', '1')
    npz = np.load(filepath)
    arr = np.sqrt(npz['arr'])
    # 2D delay
    ax = plt.subplot(gs[0:2, 0:2])
    xi = -npz['d1']/50.
    yi = -npz['d2']/50.
    zi = arr[20, 20].copy().T
    xi, yi, zi = zoom_arrs(xi, yi, zi)
    zi /= zi.max()
    levels = np.linspace(0, 1, 200)
    ax.contourf(xi, yi, zi, cmap=cmap, levels=levels)
    contour_levels = np.linspace(0.05, 1, 6)
    ax.contour(xi, yi, zi, contour_levels, colors='k')
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
    # preapre for frequency plots
    xi = normalize_frequency(npz['w1'])
    yi = normalize_frequency(npz['w2'])
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
    ax.contour(xi, yi, zi, contour_levels, colors='k')
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
    ax.contour(xi, yi, zi, contour_levels, colors='k')
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
    ax.contour(xi, yi, zi, contour_levels, colors='k')
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
    ax.contour(xi, yi, zi, contour_levels, colors='k')
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


output_path = os.path.join(directory, 'delay_space_ratios.png')

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
        filepath = template_fname.format(dpr, 'all')
        npz = np.load(filepath)
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
        ax.axvline(0, c='k', ls='-', lw=2)
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
        def plot_slice(TO):
            filepath = template_fname.format(dpr, TO)
            npz = np.load(filepath)
            arr = np.sqrt(npz['arr'][20, 20].T)
            arr /= amps_max
            ax.plot(xi, arr[10], lw=2)
        plot_slice(5)
        plot_slice(6)
        plot_slice(1)
        plot_slice(3)
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
    plot(0, dprs[0], manuals=manuals[0], yticks=True, title=r'$\mathsf{\tau = w_t/2}$')
    plot(1, dprs[1], manuals=manuals[1], yticks=False, title=r'$\mathsf{\tau = w_t}$')
    plot(2, dprs[2], manuals=manuals[2], yticks=False, title=r'$\mathsf{\tau = 2w_t}$')
    # colorbar
    ticks = np.linspace(0, 1, 11)
    cax = plt.subplot(gs[:, -1])
    matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, ticks=ticks)
    cax.set_ylabel('amplitude', fontsize=18)
    # finish
    plt.savefig(output_path, dpi=300, transparent=True, pad_inches=1)
    plt.close('fig')
