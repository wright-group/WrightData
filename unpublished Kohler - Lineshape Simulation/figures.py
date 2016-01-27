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
    plt.savefig(output_path, dpi=300, transparent=True)
    plt.close(fig)


### simulation overview #######################################################


output_path = os.path.join(directory, 'simulation_overview.png')

force_plotting = True

if not os.path.isfile(output_path) or force_plotting:
    # prepare for plot
    fig = plt.figure(figsize=[15, 8])
    gs = grd.GridSpec(2, 4, width_ratios=[20, 30, 20, 1])
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
    # TODO: use actually relevant parameters
    tau = 200.
    mu = 1.
    FWHM = np.linspace(-1.6, 1.6, num=10)
    FWHM = 10**FWHM * tau
    start_fine_int = -250
    end_fine_int = 500
    t1 = np.linspace(-10000, start_fine_int, num=10000)
    t2 = np.linspace(start_fine_int+1, end_fine_int-1, num=20000)
    t3 = np.linspace(end_fine_int,20000, num=10000)
    t = np.hstack([t1,t2,t3])
    rho1 = [rho(t,tau,mu,normalized_gauss(t,FWHMi), np.ones(len(t))) for FWHMi in FWHM]
    ax = plt.subplot(gs[0, 0])    
    # plot pulse
    xi = t
    function = wt.fit.Gaussian()
    yi = function.evaluate([0, 50, 1, 0], xi)
    plt.plot(xi, yi)
    # plot rhos
    for i, rhoi in enumerate(rho1):
        ax.plot(t,np.abs(rhoi)*2, 'k', linewidth=2, alpha=0.8)
    ax.set_xlim(-500, 500)
    ax.grid()
    wt.artists.corner_text('a', ax=ax, fontsize=16)
    # evolution of density matrix terms in pw5 --------------------------------
    ax = plt.subplot(gs[0, 1])
    trive.exp.set_coord(trive.d2, 200.)
    trive.exp.set_coord(trive.ws, 7000.)
    trive.exp.set_coord(trive.ss, 50.)
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
    # -2
    yi = np.abs(arr[0, 1, :])**2
    yi /= yi.max()
    plt.plot(yi, c='k')
    # 2'
    yi = np.abs(arr[0, 2, :])**2
    yi /= yi.max()
    plt.plot(yi, c='k')
    # 1
    yi = np.abs(arr[0, 0, :])**2
    yi /= yi.max()
    plt.plot(yi, c='k')
    # plot density matrix terms
    # ga
    yi = np.abs(scan.sig[0, 2, :])**2
    yi /= yi.max()
    plt.plot(yi)
    # aa
    yi = np.abs(scan.sig[0, 3, :])**2
    yi /= yi.max()
    plt.plot(yi)
    # ag2
    yi = np.abs(scan.sig[0, 5, :])**2
    yi /= yi.max()
    plt.plot(yi)
    # finish    
    plt.ylim(0, 1.1)
    plt.grid()    
    wt.artists.corner_text('b', ax=ax, fontsize=16)
    # FT of FID above and mono pass function ----------------------------------
    ax = plt.subplot(gs[0, 2])
    # TODO: make this correspond to actual calculations
    function = wt.fit.Gaussian()
    xi = np.linspace(-100, 100)
    yi = function.evaluate([0, 40, 1, 0], xi)
    ax.plot(xi, yi, c='k', lw=2)
    ax.axvline(-10, c='k', ls='--')
    ax.axvline(10, c='k', ls='--')
    ax.set_ylim(0, 1.1)
    ax.grid()
    wt.artists.corner_text('c', ax=ax, fontsize=16)
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
    wt.artists.diagonal_line(xi, yi)
    wt.artists.corner_text('d', ax=ax, fontsize=16)
    # representation of kernal ------------------------------------------------    
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.collections import PolyCollection
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.colors import colorConverter        
    ax = plt.subplot(gs[1, 1], projection='3d')
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    plt.setp(ax.get_zticklabels(), visible=False)
    def cc(arg):
        return colorConverter.to_rgba(arg, alpha=0.6)
    xi = np.linspace(-2, 2, 1000)
    yi = xi
    function = wt.fit.Gaussian()
    zi = function.evaluate([0, 1/3., 1, 0], xi)
    verts = [zip(xi, yi, zi)]
    poly = Poly3DCollection(verts, facecolors=[cc('grey')])
    ax.add_collection3d(poly)
    ax.set_xlabel(r'$\mathsf{\bar{\nu_1}}$', fontsize=16)
    ax.set_xlim3d(-2, 2)
    ax.set_ylabel(r'$\mathsf{\bar{\nu_2}}$', fontsize=16)
    ax.set_ylim3d(-2, 2)
    ax.set_zlim3d(0, 1)
    wt.artists.corner_text('e', ax=ax, fontsize=16)
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
    wt.artists.diagonal_line(xi, yi)
    wt.artists.corner_text('f', ax=ax, fontsize=16)
    # show inhomo linewidth applied
    center = 7000.
    width = 500.
    l = np.linspace(center-(width/2), center+(width/2), 100)
    l = normalize_frequency(l)
    ax.plot(l, l, c='k', lw=5)
    # colorbar ----------------------------------------------------------------
    plt.colorbar(mappable=mappable, cax=plt.subplot(gs[:, 3]))

    # finish
    plt.savefig(output_path, dpi=300, transparent=True)
    plt.close('fig')


### hamiltonian ###############################################################


# TODO:

