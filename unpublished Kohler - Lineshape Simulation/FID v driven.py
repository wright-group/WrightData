# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 11:48:52 2015

illustrate the regime of FID and driven

@author: Dan
"""

from time import strftime
import NISE as n
import numpy as np
import os
import matplotlib.pyplot as plt

t = n.experiments.trive
m = n.lib.measure
H = n.hamiltonians.H0
inhom = n.hamiltonians.params.inhom
folder_path = os.path.expanduser("~") + r'\Documents\NISE\1S params search'
n.lib.scan.Scan.output_folder = folder_path

pulse_width = 50.
slitwidth = 120.
dpr = 0.75
w_laser = 7050.

if __name__ ==  '__main__':

    # dephasing : pulse ratio
    #dpr = np.array([0.5, 1]).astype(float)
    
    w1 = t.w1
    w2 = t.w2
    ws = t.ws
    d2 = t.d2
    d1 = t.d1
    
    #"""
    d2_points = np.zeros((1)) #np.linspace(-100,150,num=21)
    d1_points = np.zeros((1)) #np.linspace(-100,150,num=21)
    d1.points = d1_points
    d2.points = d2_points
    
    t.exp.set_coord(w1, w_laser)
    t.exp.set_coord(w2, w_laser)
    #"""
    """
    w1_points = np.linspace(5500,8500,num=41)
    w2_points = np.linspace(5500,8500,num=41)
    d2_points = np.array([-100,-50, 0, 50, 150])
    w1.points = w1_points
    w2.points = w2_points
    d2.points = d2_points
    
    t.exp.set_coord(d1, 0.)
    #"""

    t.exp.set_coord(t.ss, pulse_width)
    t.exp.timestep = 0.5

    # definitions external to the loop
    
    inhom_object = inhom.Inhom()
    
    m.Mono.slitwidth = slitwidth
    
    H1 = H.Omega(tau_ag  = pulse_width * dpr,
                 tau_2aa = pulse_width * dpr,
                 tau_2ag = pulse_width * dpr,
                 TOs=[6])
    H1.out_group=[[1],[5],[6]]
    #H1.D = 8.
    H1.wa_central = 7000.
    # exciton-exciton coupling
    H1.a_coupling = 75. # cm-1

    H1.mu_ag =  1.0
    H1.mu_2aa = H1.mu_ag # HO approx (1.414) vs. uncorr. electron approx. (1.)

    # change late buffer
    t.exp.late_buffer = H.Omega.tau_ag * 3 * max(dpr,1.)
    t.exp.get_coords()

    out = t.exp.scan(d1, d2, H=H1, inhom_object=inhom_object)
    out.run(autosave=False, mp=False)#, chunk=True)
    tprime = np.arange(-out.early_buffer, out.late_buffer, t.exp.timestep)
    
    plt.close('all')
    f1 = plt.figure()
    ax1 = f1.add_axes([0.1, 0.1, 0.6, 0.8])
    #plt.subplot(111)
    plt.grid()
    # convert the coherence to a rotating wave for easy phase analysis
    y = out.sig[0,0,0]
    y1 = np.exp(1j*H1.wa_central*n.lib.misc.wn_to_omega*tprime)*y
    y1 = y1 / np.abs(y1).max()
    #plt.plot(tprime/pulse_width, np.real(y1),
    #         linewidth=2)
    #plt.plot(tprime/pulse_width, np.imag(y1),
    #         linewidth=2)
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
    y2 /= n.lib.misc.wn_to_omega
    y3 /= out.timestep * n.lib.misc.wn_to_omega

    #y3[y3>1] = 0.
    y3 -= H1.wa_central
    y3 *= -1
    plt.ylabel(r'$\mathsf{|\rho_{10}(t) | (a.u.)}$', fontsize=14)
    plt.yticks([0,0.5,1.])
    plt.xlabel(r'$\mathsf{t / \tau_{pulse}}$', fontsize=14)
    color = (y3 - H1.wa_central) / (w_laser - H1.wa_central)
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    coolwarm = cm = plt.get_cmap('coolwarm') 
    cNorm  = colors.Normalize(vmin=0, vmax=1.)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    for i in range(len(tprime[1:])):
        # calculate the color to use, given the instantaneous frequency
        colorVal = scalarMap.to_rgba(color[i])
        plt.plot([tprime[i+1]/pulse_width,tprime[i+1]/pulse_width], [0,np.abs(y1)[i]],
                 color=colorVal)
    plt.plot(tprime/pulse_width, np.abs(y1)/np.abs(y1).max(),
             linewidth=4., color='k')
    plt.ylim(0.,1.1)

    from matplotlib import colorbar
    ax2 = f1.add_axes([0.75, 0.1, 0.025, 0.8])
    cb1 = colorbar.ColorbarBase(ax2, cmap='coolwarm', ticks=[0,1.])
    plt.ylabel('$\mathsf{d\phi / dt (cm^{-1})}$', fontsize=14)
    plt.yticks(np.array([0,1.]), 
               (r'$\mathsf{\omega_{10}}$',
                r'$\mathsf{\omega_{laser}}$'), fontsize=16)
               #[min(H1.wa_central,w_laser),
               # max(H1.wa_central,w_laser)])

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
    
