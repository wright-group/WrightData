### import ####################################################################


import os
import itertools

import matplotlib
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

props = dict(boxstyle='square', facecolor='white', alpha=0.8)

template_fname = r'npzs\dpr {0} TO {1} w1 w2 d1 d2 arr.npz'
dprs = [0.5, 1., 2.0]


def normalize_frequency(arr):
    arr -= arr.mean()
    arr /= 2 * np.log(2) / (50*np.pi * 3e10) * 1e15
    return arr


def normalize_delay(arr):
    arr /= 50.
    return arr


def time_label():
    return r'$\mathsf{t/w_t}$'


def frequency_label(sub=''):
    return r'$\mathsf{\nu_{' + sub + r'}/w_{\nu}}$'


freq_ticks = [-1, 0, 1]

### WMELs #####################################################################


output_path = os.path.join(directory, 'WMELs.png')

force_plotting = False

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
    # TODO: add dot dot dot, line for n=N
    state_names = ['n=0', '1', '2', '3']
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
    fig.text(0.075, 1-dist, 'a', bbox=props, horizontalalignment='left', verticalalignment='top', fontsize=16)
    fig.text(0.35, 1-dist, 'b', bbox=props, horizontalalignment='left', verticalalignment='top', fontsize=16)    
    # finish
    plt.savefig(output_path, dpi=300, transparent=True, bbox_inches=0.5)
    plt.close(fig)


### simulation overview #######################################################


output_path = os.path.join(directory, 'simulation_overview.png')

force_plotting = True

if not os.path.isfile(output_path) or force_plotting:
    # prepare for plot
    fig, gs = wt.artists.create_figure(width='double', nrows=2, cols=[1, 1.5, 1, 'cbar'],
                                       hspace=0.75)
    # set of coherence times vs lab time --------------------------------------
    # calculating transients directly right here is easiest
    def normalized_gauss(t, FWHM):  
        sigma = FWHM / (2.*np.sqrt(np.log(2.)))    
        out = np.exp(-0.5*(t/sigma)**2)
        out /= sigma * np.sqrt(2*np.pi)
        return out        
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
    ax.legend(fontsize=16, loc='upper right')
    ax.set_ylim(0, 1.15)
    ax.set_xlim(-3, 11)
    ax.set_xlabel(time_label(), fontsize=16)
    ax.grid()
    wt.artists.corner_text('a', ax=ax, fontsize=16, background_alpha=1)
    ax.set_ylabel('amplitude', fontsize=16)
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
    ax.text(-4, 0.1, '2', fontsize=16, ha='center', va='center')
    ax.text(-4, 0.3, 'x', fontsize=16, ha='center', va='center')
    # 2'
    yi = np.abs(arr[0, 2, :])**2
    yi /= yi.max()
    ax.fill_between(xi, 0, yi, facecolor='b', alpha=0.25)
    ax.text(-2, 0.1, '2\'', fontsize=16, ha='center', va='center')
    ax.text(-2, 0.3, 'y', fontsize=16, ha='center', va='center')
    # 1
    yi = np.abs(arr[0, 0, :])**2
    yi /= yi.max()
    ax.fill_between(xi, 0, yi, facecolor='r', alpha=0.25)
    ax.text(0, 0.1, '1', fontsize=16, ha='center', va='center')
    ax.text(0   , 0.3, 'z', fontsize=16, ha='center', va='center')
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
    ax.text(-4.1, 1.05, r'$\mathsf{\tau_{22^{\prime}}}$', ha='right', va='center', fontsize=16)
    ax.text(0.1, 1.1, r'$\mathsf{\tau_{21}}$', ha='left', va='center', fontsize=16)
    # finish
    ax.legend(fontsize=16, loc='right')
    ax.set_ylim(0, 1.15)
    ax.set_xlim(-7, 7)
    plt.grid()    
    wt.artists.corner_text('b', ax=ax, fontsize=16, background_alpha=1)
    ax.set_xlabel(time_label(), fontsize=16)
    plt.setp(ax.get_yticklabels(), visible=False)
    # FT of FID above and mono pass function ----------------------------------
    ax = plt.subplot(gs[0, 2])
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
    wt.artists.corner_text('c', ax=ax, fontsize=16, background_alpha=1)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.set_xlabel(frequency_label('out'), fontsize=16)
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
    wt.artists.corner_text('d', ax=ax, fontsize=16, background_alpha=1)
    ax.set_xlabel(frequency_label('1'), fontsize=16)
    ax.set_ylabel(frequency_label('2'), fontsize=16)
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
    ax.set_xlabel(r'$\mathsf{\nu_1}$', fontsize=16)
    ax.set_xlim3d(-2, 2)
    ax.yaxis.set_rotate_label(False)  # disable automatic rotation
    ax.set_ylabel(r'$\mathsf{\nu_2}$', fontsize=16, rotation=115)
    ax.set_ylim3d(-2, 2)
    ax.set_zlim3d(0, 1)
    wt.artists.corner_text('e', ax=ax, fontsize=16, background_alpha=1)
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
    wt.artists.corner_text('f', ax=ax, fontsize=16, background_alpha=1)
    # show inhomo linewidth applied
    center = 7000.
    width = 500.
    l = np.linspace(center-(width/2), center+(width/2), 100)
    l = normalize_frequency(l)
    ax.plot(l, l, c='k', lw=5)
    ax.set_xlabel(frequency_label('1'), fontsize=16)
    ax.set_ylabel(frequency_label('2'), fontsize=16)
    # colorbar ----------------------------------------------------------------
    ticks = np.linspace(0, 1, 11)
    cax = plt.subplot(gs[:, -1])
    matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, ticks=ticks)
    cax.set_ylabel('amplitude', fontsize=16)
    if False:    
        plt.suptitle('NISE overview', fontsize=20)
        wt.artists.plot_margins(fig)
    # finish
    plt.savefig(output_path, dpi=300, transparent=True, pad_inches=1)
    plt.close('fig')


### hamiltonian ###############################################################


# TODO:

