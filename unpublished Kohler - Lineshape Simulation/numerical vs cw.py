# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 17:40:30 2016

for comparing numerical solutions at special coordinates with 
calculations under simple conditions

@author: Dan
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import fscolors as f
from scipy.ndimage.interpolation import zoom
#from scipy.signal import fftconvolve

data_folder = os.getcwd()

# wavenumbers to relative units 
# supposedly this is the fwhm bandwidth (intensity scale) a 50fs pulse
wn_to_rel = (2 * np.log(2) / (50e-15*np.pi * 3e10))**-1
# omega / (2pi * c) = wavenumber
G_factor = (2 * np.pi * 3e10)**-1
# use either 1/50, 1/100, or 1/25 
Gamma_10 = (50e-15)**-1 * G_factor
Gamma_11 = (50e-14)**-1 * G_factor # keep fixed at the pulse width

Gamma_10 *= wn_to_rel
Gamma_11 *= wn_to_rel

def import_npzs():
    fname = r'\smeared fwhm_maj 0 dpr TO w1 w2 d1 d2 arr.npz'
    x = np.load(data_folder+r'\data'+fname)
    global w1
    global w2
    global d1
    global d2
    w1 = x['w1']
    w2 = x['w2']
    w1 = w1 - w1.mean()
    w2 = w2 - w2.mean()
    #w1 /= 2 * np.log(2) / (50*np.pi * 3e10) * 1e15
    #w2 /= 2 * np.log(2) / (50*np.pi * 3e10) * 1e15
    w1 *= wn_to_rel
    w2 *= wn_to_rel
    d1 = x['d1'] / 50.
    d2 = x['d2'] / 50.
    # keys:  dpr,TO,w1,w2,d1,d2,arr
    print 'array imported'
    return x['arr']

data = import_npzs()

# redefine w1 and w2 to have twice as many elements
w = np.linspace(w1.min(), w1.max(), num = w1.size * 3)
# keys:  dpr,TO,w1,w2,d1,d2,arr
num_sim = data[1,-1] # dpr index 0 is CW-ish, dpr index 2 is impulsive-ish

### define the cw equations
def gen_f_cw(l,  # l = lambda1 * lambda2 * lambda3 determines the sign
         kappas, # kappas = [k1,k2,k3]
         gammas, # gammas = [G1,G2,G3],
         ws      # ws = [w1,w2,w3] are FID frequencies
    ):
    def f(omegas # omegas = [omegax, omegay, omegaz] are efield frequencies
        ):
        out = (gammas[2] + 1j*(kappas[0]*omegas[0]+kappas[1]*omegas[1]+kappas[2]*omegas[2]-ws[0]-ws[1]-ws[2]) )
        out *= (gammas[1] + 1j*(kappas[0]*omegas[0]+kappas[1]*omegas[1]                    -ws[0]-ws[1]) )
        out *= (gammas[0] + 1j*(kappas[0]*omegas[0]                                        -ws[0]) )
        out **= -1
        return out
    return f

def f_cw(omegas, # omegas = [omegax, omegay, omegaz] are efield frequencies
         l, # l = lambda1 * lambda2 * lambda3 determines the sign
         kappas, # kappas = [k1,k2,k3]
         gammas, # gammas = [G1,G2,G3],
         ws, # ws = [w1,w2,w3] are FID frequencies
    ):
    out = np.zeros((w.size, w.size, w.size), dtype=np.complex)
    out += (gammas[2] + 1j*(kappas[0]*omegas[0]+kappas[1]*omegas[1]+kappas[2]*omegas[2]-ws[0]-ws[1]-ws[2]) )
    out *= (gammas[1] + 1j*(kappas[0]*omegas[0]+kappas[1]*omegas[1]                    -ws[0]-ws[1]) )
    out *= (gammas[0] + 1j*(kappas[0]*omegas[0]                                        -ws[0]) )
    out **= -1
    return out

### apply pulse widths 
# to avoid error, to explicit integrals for now--don't need to be fast
def C(w):
    """
    spectral envelope
    """
    return np.exp(-2*np.sqrt(2)*np.log(2)*w**2)

def pixel_average(w1,w2, Gamma_11, V_VI_only=False):
    """
    generate the signal for this w1,w2 coordinate
    integrate over spectral degrees of freedom
    """
    b = np.linspace(-2,2,num=21)
    c = b.copy()
    # out[a,b]
    alpha_nr1 = gen_f_cw(1,
                   [-1,  1, -1],
                   [Gamma_10, Gamma_11, Gamma_10],
                   [0,0,0])
    alpha_nr2 = gen_f_cw(1,
                   [-1,  1, -1],
                   [Gamma_10, Gamma_11, Gamma_10],
                   [0,0,0])
    alpha_r1  = gen_f_cw(1,
                   [1, -1, -1],
                   [Gamma_10, Gamma_11, Gamma_10],
                   [0,0,0])
    alpha_r2  = gen_f_cw(1,
                   [1, -1, -1],
                   [Gamma_10, Gamma_11, Gamma_10],
                   [0,0,0])
    # in order of b,c
    alpha_5 =  alpha_r2([w2+b[:,None],w2+c[None,:],w1+b[:,None]-c[None,:]])
    alpha_6 = alpha_nr2([w2+c[None,:],w2+b[:,None],w1+b[:,None]-c[None,:]])
    alpha_1 = alpha_nr1([w1+b[:,None]-c[None,:],w2+b[:,None],w2+c[None,:]])
    alpha_3 =  alpha_r1([w2+b[:,None],w1+b[:,None]-c[None,:],w2+c[None,:]]) 
    if V_VI_only:
        out = alpha_5 + alpha_6
    else:
        out = alpha_1 + alpha_3 + alpha_5 + alpha_6
    out *= C(b[:,None]-c[None,:]) * C(b[:,None]) * C(c[None,:])
    return out.sum().sum()

convolved1 = np.zeros((w.size,w.size),dtype=np.complex)
convolved2 = convolved1.copy()
convolved3 = convolved1.copy()
convolved4 = convolved1.copy()
for ind in np.ndindex(convolved1.shape):
    convolved1[ind] = pixel_average(w[ind[0]], w[ind[1]], 
                                    (50e-12)**-1 * G_factor * wn_to_rel)
    convolved2[ind] = pixel_average(w[ind[0]], w[ind[1]], 
                                    (50e-15)**-1 * G_factor * wn_to_rel)
    convolved3[ind] = pixel_average(w[ind[0]], w[ind[1]], 
                                    (50e-12)**-1 * G_factor * wn_to_rel,
                                    V_VI_only=True)
    convolved4[ind] = pixel_average(w[ind[0]], w[ind[1]], 
                                    (50e-15)**-1 * G_factor * wn_to_rel,
                                    V_VI_only=True)
                                    
convolved1 /= np.abs(convolved1).max()
convolved2 /= np.abs(convolved2).max()
convolved3 /= np.abs(convolved3).max()
convolved4 /= np.abs(convolved4).max()

# wx,wy,wz
alpha_nr1 = f_cw([w[:,None,None],w[None,:,None],w[None,None,:]],
               1,
               [-1,  1, -1],
               [Gamma_10, Gamma_11, Gamma_10],
               [0,0,0])
alpha_nr2 = f_cw([w[:,None,None],w[None,:,None],w[None,None,:]],
               1,
               [-1,  1, -1],
               [Gamma_10, Gamma_11, Gamma_10],
               [0,0,0])
alpha_r1 = f_cw([w[:,None,None],w[None,:,None],w[None,None,:]],
               1,
               [1, -1, -1],
               [Gamma_10, Gamma_11, Gamma_10],
               [0,0,0])
alpha_r2 = f_cw([w[:,None,None],w[None,:,None],w[None,None,:]],
               1,
               [1, -1, -1],
               [Gamma_10, Gamma_11, Gamma_10],
               [0,0,0])

# order this as w1,w2
alpha_1 = alpha_nr1.transpose(0,1,2).diagonal(axis1=1,axis2=2)
alpha_3 =  alpha_r1.transpose(1,0,2).diagonal(axis1=1,axis2=2) 
alpha_5 =  alpha_r2.transpose(2,0,1).diagonal(axis1=1,axis2=2)
alpha_6 = alpha_nr2.transpose(2,1,0).diagonal(axis1=1,axis2=2)

### Case 1:  compare pulse overlap 2D frequency, no inhomogeneity
plt.close('all')
ticks=np.linspace(-2,2,num=5)
f1 = plt.figure(figsize=(14,8))

gs1 = GridSpec(4,4)

ax241 = f1.add_subplot(gs1[0:2,0], aspect='equal')
cw_sim = alpha_1 + alpha_3 + alpha_5 + alpha_6
cw_sim /= np.abs(cw_sim).max()
w_zoom_cw = zoom(w,1)
convolved1_zoom = zoom(np.abs(convolved1),1)
plt.contourf(w_zoom_cw, w_zoom_cw, np.abs(convolved1_zoom.T), 
             levels=np.linspace(0,1,num=256), cmap=f.cm.chw2)
plt.contour(w_zoom_cw, w_zoom_cw, np.abs(convolved1_zoom.T), 
            levels=np.linspace(0,1,num=6)[1:-1], colors='k', linewidths=2)
plt.plot(ticks, ticks, 'k', linestyle='-', linewidth=1)
plt.plot(ticks, -ticks, 'k', linestyle='--', linewidth=1)
plt.title(r'Steady State, $\Gamma_{11} \rightarrow \infty$')
plt.xlim(-2,2)
plt.ylim(-2,2)
plt.xticks(ticks)
plt.yticks(ticks)
plt.grid(b=True)


ax242 = f1.add_subplot(gs1[0:2,1], aspect='equal')
plt.contourf(w, w, np.abs(convolved2.T),
    levels=np.linspace(0,1,num=256), cmap=f.cm.chw2)
plt.contour(w, w, np.abs(convolved2.T),
    levels=np.linspace(0,1,num=6)[1:-1], colors='k', linewidths=2)
plt.plot(ticks, ticks, 'k', linestyle='-', linewidth=1)
plt.plot(ticks, -ticks, 'k', linestyle='--', linewidth=1)
plt.title('Steady State, $\Gamma_{11} \Delta_t = 1$')
plt.xlim(-2,2)
plt.ylim(-2,2)
plt.xticks(ticks)
plt.yticks(ticks)
plt.grid(b=True)


ax243 = f1.add_subplot(gs1[0:2,2], aspect='equal')
# keys:  w1,w2,d1,d2,arr
num_sim1 = zoom(num_sim[:,:,10,10], 3)**0.5
num_sim1 /= np.abs(num_sim1).max()
w_zoom = np.linspace(w1.min(), w1.max(), num = w1.size * 3)
plt.contourf(w_zoom, w_zoom, np.abs(num_sim1.T), 
             levels=np.linspace(0,1,num=256), cmap=f.cm.chw2)
plt.contour(w_zoom, w_zoom, np.abs(num_sim1.T), 
            levels=np.linspace(0,1,num=6)[1:-1], colors='k', linewidths=2)
plt.plot(ticks, ticks, 'k', linestyle='-', linewidth=1)
plt.plot(ticks, -ticks, 'k', linestyle='--', linewidth=1)
plt.title('Actual')
plt.xlim(-2,2)
plt.ylim(-2,2)
plt.xticks(ticks)
plt.yticks(ticks)
plt.grid(b=True)

f1.add_subplot(gs1[0,3])
#plt.title('diagonal traces')
# keys:  dpr,TO,w1,w2,d1,d2,arr
plt.plot(w_zoom, np.abs(num_sim1).diagonal(),
         color='b', linewidth=3, linestyle='-', alpha=0.5)
plt.plot(w, np.abs(convolved2).diagonal(),
         color='g', linewidth=3, linestyle='-', alpha=0.5)
plt.plot(w, np.abs(convolved1).diagonal(),
         color='r', linewidth=3, linestyle='-', alpha=0.5)
plt.xlim(-2,2)
plt.ylim(-.1,1.1)
plt.xticks(ticks)
plt.grid()

f1.add_subplot(gs1[1,3])
#plt.title('antidiagonal traces')
plt.plot(w_zoom, np.abs(num_sim1[::-1]).diagonal(),
         color='b', linewidth=3, linestyle='--', alpha=0.5)
plt.plot(w, np.abs(convolved1[::-1]).diagonal(),
         color='r', linewidth=3, linestyle='--', alpha=0.5)
plt.plot(w, np.abs(convolved2[::-1]).diagonal(),
         color='g', linewidth=3, linestyle='--', alpha=0.5)
plt.xlim(-2,2)
plt.ylim(-.1,1.1)
plt.xticks(ticks)
#plt.yticks(ticks)
plt.grid(b=True)


### Case 2:  pathway 5/6 (2 and 2' overlapped), no inhomogeneity
ax245 = f1.add_subplot(gs1[2:,0], aspect='equal')
plt.contourf(w, w, np.abs(convolved3.T), 
             levels=np.linspace(0,1,num=256), cmap=f.cm.chw2)
plt.contour(w, w, np.abs(convolved3.T), 
            levels=np.linspace(0,1,num=6)[1:-1], colors='k', linewidths=2)
plt.plot(ticks, ticks, 'k', linestyle='-', linewidth=1)
plt.plot(ticks, -ticks, 'k', linestyle='--', linewidth=1)
plt.xlim(-2,2)
plt.ylim(-2,2)
plt.xticks(ticks)
plt.yticks(ticks)
plt.grid(b=True)

ax246 = f1.add_subplot(gs1[2:,1], aspect='equal')
plt.contourf(w, w, np.abs(convolved4.T), 
             levels=np.linspace(0,1,num=256), cmap=f.cm.chw2)
plt.contour(w, w, np.abs(convolved4.T), 
            levels=np.linspace(0,1,num=6)[1:-1], colors='k', linewidths=2)
plt.plot(ticks, ticks, 'k', linestyle='-', linewidth=1)
plt.plot(ticks, -ticks, 'k', linestyle='--', linewidth=1)
plt.xlim(-2,2)
plt.ylim(-2,2)
plt.xticks(ticks)
plt.yticks(ticks)
plt.grid(b=True)


ax247 = f1.add_subplot(gs1[2:,2], aspect='equal')
# keys:  w1,w2,d1,d2,arr
num_sim2 = zoom(num_sim[:,:,10,-1]**0.5,3)
num_sim2 /= np.abs(num_sim2).max()
plt.contourf(w_zoom, w_zoom, np.abs(num_sim2.T), 
             levels=np.linspace(0,1,num=256), cmap=f.cm.chw2)
plt.contour(w_zoom, w_zoom, np.abs(num_sim2.T), 
            levels=np.linspace(0,1,num=6)[1:-1], colors='k', linewidths=2)
plt.plot(ticks, ticks, 'k', linestyle='-', linewidth=1)
plt.plot(ticks, -ticks, 'k', linestyle='--', linewidth=1)
plt.xlim(-2,2)
plt.ylim(-2,2)
plt.xticks(ticks)
plt.yticks(ticks)
plt.grid(b=True)

# keys:  dpr,TO,w1,w2,d1,d2,arr
f1.add_subplot(gs1[2,3])
plt.plot(w_zoom, np.abs(num_sim2).diagonal(),
         color='b', linewidth=3, linestyle='-', alpha=0.5)
plt.plot(w, np.abs(convolved3).diagonal(),
         color='r', linewidth=3, linestyle='-', alpha=0.5)
plt.plot(w, np.abs(convolved4).diagonal(),
         color='g', linewidth=3, linestyle='-', alpha=0.5)
plt.xlim(-2,2)
plt.ylim(-.1,1.1)
plt.xticks(ticks)
#plt.yticks(ticks)
plt.grid(b=True)
f1.add_subplot(gs1[3,3])
plt.plot(w_zoom, np.abs(num_sim2[::-1]).diagonal(),
         color='b', linewidth=3, linestyle='--', alpha=0.5)
plt.plot(w, np.abs(convolved3[::-1]).diagonal(),
         color='r', linewidth=3, linestyle='--', alpha=0.5)
plt.plot(w, np.abs(convolved4[::-1]).diagonal(),
         color='g', linewidth=3, linestyle='--', alpha=0.5)
plt.xlim(-2,2)
plt.ylim(-.1,1.1)
plt.xticks(ticks)
#plt.yticks(ticks)
plt.grid(b=True)

for ax, color in zip([ax241, ax245, ax242, ax246, ax243, ax247], 
                     ['red', 'red', 'green', 'green', 'blue', 'blue']):
    plt.setp(ax.spines.values(), color=color, linewidth=3, alpha=0.7)

f1.savefig(r'figures\num vs cw pulse overlap.png', transparent=True, dpi=200)

### Case 3:  pathway 1, no overlap, no inhomogeneity (compare with impulsive)



