# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 12:41:07 2015

make figure for M-factor corrections:
show for fsb 19 and fsb 25 to show how both can be made to have relative 
agreement

@author: Dan
"""

import fscolors as f
import numpy as np
import os
from scipy.interpolate import interp1d

xticks = np.array([6500, 7500, 8500])
yticks = xticks

path = os.path.expanduser(r'~')
path1 = '\\'.join([path, 
           'Desktop', 
           'incoming from fs table',
           'fsb 25',
           '2014.09.11'])
path2 = '\\'.join([path, 
           'Desktop', 
           'incoming from fs table',
           'fsb 19',
           '2014.09.14'])

#-------------- fsb 25 raw
zi00_file = '\\'.join([path1, 
           'tot rank 81 raw NRMSD 1 0.00%.npz'])
npz00 = np.load(zi00_file)
d2 = npz00['d2']
w1 = npz00['w1'][1:]
w2 = npz00['w2']
zi00 = npz00['zis']

#-------------- fsb 25 without m-factors
zi01_file = '\\'.join([path1, 
           'tot rank 81 without m-factors NRMSD 0.00%.npz'])
npz01 = np.load(zi01_file)
zi01 = npz01['arr_3']

#-------------- fsb 25 with M-factors
zi02_file = os.path.expanduser(r'~')
zi02_file = '\\'.join([path1, 
           'tot rank 81 with m-factors NRMSD 0.00%.npz'])
npz02 = np.load(zi02_file)
zi02 = npz02['arr_3']

#-------------- fsb 19 raw
zi10_file = '\\'.join([path2,
                       'tot raw rank 54 NRMSD 0.00%.npz'])
npz10 = np.load(zi10_file)
w11 = npz10['w1'][1:]
w21 = npz10['w2']
d21 = npz10['d2']
zi10 = npz10['zi']
#-------------- fsb 19 with M-factors
zi11_file = '\\'.join([path2,
                       'tot without m-factor rank 54 NRMSD 0.00%.npz'])
npz11 = np.load(zi11_file)
zi11 = npz11['zi']
#-------------- fsb 19 without M-factors
zi12_file = '\\'.join([path2,
                       'tot with m-factor rank 54 NRMSD 0.00%.npz'])
npz12 = np.load(zi12_file)
zi12 = npz12['zi']


abs_file = os.path.expanduser(r'~')
abs_file = '\\'.join([abs_file,
           'Documents', 
           'Colloidal Synthesis Characterization',
           'DK-2014.05.21-PbSe 025'])
abs_file += r'\DK-2014.09.17-PbSe 025.4.a.txt'

l,a = np.loadtxt(abs_file, unpack = True, skiprows=18)
A = interp1d(1e7/l, a)

def M(w1,w2):
    a1 = A(w1)
    a2 = A(w2)
    out = 10**(-a1[None,:] / 2.) * (1 - 10**(-a2[:,None])) #* w1[None,:]
    out /= a2[:,None] * np.log(10)
    return np.abs(out)

M12 = M(w1,w2)
M12max = M12.max()
M12 /= M12max

def process_zi(arr, i):
    zi = np.ma.masked_less(arr[i],0)
    zimax = zi.max()
    zi /= zimax
    zi = zi.filled(0)
    return zi

z00 = process_zi(zi00,-45)[1:]
z01 = process_zi(zi02,-45)[1:]
z02 = process_zi(zi01,-45)[1:]
z10 = process_zi(zi10,-1)[1:]
z11 = process_zi(zi11,-1)[1:]
z12 = process_zi(zi12,-1)[1:]

import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size':16
})
alpha=0.3
linewidth=4
plt.figure(figsize=(10,7))
#-------------- 131
plt.subplot(231, aspect='equal')
plt.contourf(w1, w2, np.sqrt(z00.T), 
             levels=np.linspace(0,1,num=256),
             cmap=f.Dat.mycm)
plt.contour(w1, w2, np.sqrt(z00.T), 
             levels=np.linspace(0-1e-6,1,num=11),
             colors='k')
plt.plot([min(w1), max(w1)],[6620,6620], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([min(w1), max(w1)],[8340,8340], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([6620,6620],[min(w2), max(w2)], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([8340,8340],[min(w2), max(w2)], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([max(min(w1),min(w2)), min(max(w1),max(w2))],
          [max(min(w1),min(w2)), min(max(w1),max(w2))], 'k:')
plt.title(r'$\sqrt{S_{\mathsf{TrEE}}^{(0)}} $', fontsize=16)
plt.yticks(yticks)
plt.ylabel(r'$\mathsf{\bar\nu_2 (cm^{-1})}$', fontsize=16)
plt.xticks(xticks, visible=False)
plt.grid(b=True)

#-------------- 232
plt.subplot(232, aspect='equal')
plt.contourf(w1,w2,np.sqrt(z01.T), 
             levels=np.linspace(0,1,num=256),
             cmap=f.Dat.mycm)
plt.contour(w1, w2, np.sqrt(z01.T), 
             levels=np.linspace(0-1e-6,1,num=11),
             colors='k')
plt.plot([min(w1), max(w1)],[6620,6620], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([min(w1), max(w1)],[8340,8340], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([6620,6620],[min(w2), max(w2)], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([8340,8340],[min(w2), max(w2)], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([max(min(w1),min(w2)), min(max(w1),max(w2))],
          [max(min(w1),min(w2)), min(max(w1),max(w2))], 'k:')
plt.title(r'$\sqrt{S_{\mathsf{TrEE}}^{(0)} / I_1 I_2^2} $', fontsize=16)
plt.yticks(yticks, visible=False)
plt.xticks(xticks, visible=False)
plt.grid(b=True)

#-------------- 233
plt.subplot(233, aspect='equal')
plt.contourf(w1,w2,np.sqrt(z02.T),
             levels=np.linspace(0,1,num=256),
             cmap=f.Dat.mycm)
plt.contour(w1, w2, np.sqrt(z02.T), 
             levels=np.linspace(0-1e-6,1,num=11),
             colors='k')
plt.plot([min(w1), max(w1)],[6620,6620], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([min(w1), max(w1)],[8340,8340], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([6620,6620],[min(w2), max(w2)], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([8340,8340],[min(w2), max(w2)], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([max(min(w1),min(w2)), min(max(w1),max(w2))],
          [max(min(w1),min(w2)), min(max(w1),max(w2))], 'k:')
plt.title(r'$S_{\mathsf{TrEE}}$', fontsize=16)
plt.yticks(yticks, visible=False)
plt.xticks(xticks, visible=False)
plt.grid(b=True)
#-------------- 234
plt.subplot(234, aspect='equal')
plt.contourf(w11, w21, np.sqrt(z10.T), 
             levels=np.linspace(0,1,num=256),
             cmap=f.Dat.mycm)
plt.contour(w11, w21, np.sqrt(z10.T), 
             levels=np.linspace(0-1e-6,1,num=11),
             colors='k')
plt.plot([min(w11), max(w11)],[7570,7570], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([7570,7570],[min(w21), max(w21)], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([max(min(w11),min(w21)), min(max(w11),max(w21))],
          [max(min(w11),min(w21)), min(max(w11),max(w21))], 'k:')
#plt.title(r'$\sqrt{S_{\mathsf{TrEE}}^{(0)}} $', fontsize=16)
plt.yticks(yticks)
plt.ylabel(r'$\mathsf{\bar\nu_2 (cm^{-1})}$', fontsize=16)
plt.xticks(xticks)
plt.grid(b=True)

#-------------- 235
plt.subplot(235, aspect='equal')
plt.contourf(w11,w21,np.sqrt(z11.T), 
             levels=np.linspace(0,1,num=256),
             cmap=f.Dat.mycm)
plt.contour(w11, w21, np.sqrt(z11.T), 
             levels=np.linspace(0-1e-6,1,num=11),
             colors='k')
plt.plot([min(w11), max(w11)],[7570,7570], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([7570,7570],[min(w21), max(w21)], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([max(min(w11),min(w21)), min(max(w11),max(w21))],
          [max(min(w11),min(w21)), min(max(w11),max(w21))], 'k:')
plt.yticks(yticks, visible=False)
plt.xticks(xticks)
plt.grid(b=True)
plt.xlabel(r'$\mathsf{\bar\nu_1 (cm^{-1})}$', fontsize=16)
#-------------- 236
plt.subplot(236, aspect='equal')
plt.contourf(w11,w21,np.sqrt(z12.T),
             levels=np.linspace(0,1,num=256),
             cmap=f.Dat.mycm)
plt.contour(w11, w21, np.sqrt(z12.T), 
             levels=np.linspace(0-1e-6,1,num=11),
             colors='k')
plt.plot([min(w11), max(w11)],[7570,7570], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([7570,7570],[min(w21), max(w21)], linewidth=linewidth, color='k', alpha=alpha)
plt.plot([max(min(w11),min(w21)), min(max(w11),max(w21))],
          [max(min(w11),min(w21)), min(max(w11),max(w21))], 'k:')
plt.yticks(yticks, visible=False)
plt.xticks(xticks)
plt.grid(b=True)

plt.tight_layout()

#plt.savefig('M-factors.png', transparent=True, dpi=300)


