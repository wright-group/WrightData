# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 19:47:52 2016

look at measurable properties of inhomogeneity in the samples

@author: Dan
"""

import numpy as np
import matplotlib.pyplot as plt
import os

inhom = np.array([0,0.5,1,2])
dprs = np.array([0.5,1,2])
data_folder = os.getcwd() + r'\data'
fname = r'\smeared fwhm_maj {0} dpr TO w1 w2 d1 d2 arr.npz'

def load_npz():
    arr = []
    for inhom_i in inhom:
        fname_i = fname.format(inhom_i)
        dic = np.load(data_folder + fname_i)
        arr.append(dic['arr'])
    return np.array(arr)

d1 = np.linspace(-200,200, num=21)
d2 = d1.copy()
w1 = np.linspace(5500,8500,num=41)
w2 = w1.copy()
# arr[inhom,dpr,TO,w1,w2,d1,d2]
try: arr.size # a fairly good check that we have the right arr
except NameError:
    arr = load_npz()

### analyze 3peps 
# select the array on resonance, at the latest 2d value
# as of right now, analyzing on intensity scale, because that is
# how it is typically analyzed
arr1 = arr[:,:,:,20,20,:,-1]

# go by fitting?
def gauss(p, x):
    A, mu, w = p
    z = (mu - x) / (np.sqrt(2) * w)
    out = A * np.exp(-z**2)
    return out

fits = np.zeros((4,3,3))

def erf1(p,x,y):
    return y - gauss(p,x)

from scipy.optimize import leastsq

for ind in np.ndindex(arr1.shape[:-2]):
    A0 = arr1[ind][-1].max()
    mu0 = 0
    w0 = 40.
    fits[ind] = leastsq(erf1, [A0,mu0,w0], args=(d1, arr1[ind][-1]), 
                   full_output=False)[0]

# second param (fits[...,1]) is the peak shift we want
"""
# go by moment?
mom1 = np.zeros(arr1.shape[:-2])
for ind in np.ndindex(arr1.shape[:-2]):
    mom1[ind] = np.dot(d1,arr1[ind][-1]) / arr1.sum()
"""
colors=['k','r','g','b']
plt.close('all')
plt.figure()
plt.subplot(121)
for i in range(inhom.size):
    plt.scatter(dprs,#+x[i] , 
                fits[i,:,1], s=100, alpha=0.7, color=colors[i])
    plt.plot(dprs,#+x[i] , 
                fits[i,:,1], linewidth=2, alpha=0.7, color=colors[i])
plt.xlabel(r'$\mathsf{\Gamma_{10} / w}$')
plt.ylabel(r'Peak Shift (fs)')
plt.xlim([0,2.25])
plt.grid()
plt.subplot(122)
for i in range(dprs.size):
    plt.scatter(inhom, fits[:,i,1], s=100, alpha=0.7, color=colors[i])
    plt.plot(inhom, fits[:,i,1], linewidth=2, alpha=0.7, color=colors[i])
plt.xlabel(r'$\mathsf{\sigma / w}$')
#plt.ylabel(r'Peak Shift (fs)')
plt.xlim([0,2.25])
plt.grid()

### calculate ellipticity
# fit gaussian along antidiagonal and diagonal of 2D freq plots
# only interested at tau_22=0?
# analyzing on amplitude scale, because this is how we usually do this
arr2 = arr[:,:,-1,:,:,10,-1]**0.5
diags = arr2.diagonal(axis1=-1, axis2=-2)
antidiags = arr2[...,::-1].diagonal(axis1=-1,axis2=-2)

ad_fits = np.zeros((inhom.size,dprs.size,3))
d_fits = np.zeros((inhom.size,dprs.size,3))
for ind in np.ndindex(diags.shape[:-1]):
    print ind
    A0 = diags[ind].max()
    mu0 = w1.mean()
    w0 = 100.
    d_fits[ind] = leastsq(erf1, [A0,mu0,w0], args=(w1, diags[ind]), 
                   full_output=False)[0]
    ad_fits[ind] = leastsq(erf1, [A0,mu0,w0], args=(w1, antidiags[ind]), 
                   full_output=False)[0]

ellipticity = d_fits[...,-1]**2 - ad_fits[...,-1]**2
ellipticity /= d_fits[...,-1]**2 + ad_fits[...,-1]**2

plt.figure()
plt.subplot(111, aspect='equal')
st_ratio = np.dot(inhom[:,None], dprs[None,::-1]**-1) # tau_10 / sigma (w/W)
for i in range(3):
    plt.scatter(ellipticity[:,i], fits[:,i,1]/50, 
                s=100*st_ratio[:,i], alpha=0.7, color='k')
    # lines track constant dpr
    plt.plot(ellipticity[:,i],
             fits[:,i,1] / 50, 
             linewidth=2, alpha=0.7, color=colors[i+1],
             label=r'$\tau_{10} = $' + r'$ {0} w$'.format(dprs[i]))
for i in [0.5, 1, 2]: # places where our st ratio has multiple nontrival locations
    index = np.where(st_ratio==i)
    xi = ellipticity[index]
    yi = fits[...,1] 
    yi = yi[index]/50.
    print xi.shape
    print yi.shape
    plt.plot(xi,yi,
             color='k', linewidth=2, alpha=0.1)

# zero delay inhomogeneity characteristics--from "inhomogeneity vs. delay.py"
x_zd = [0.487957698403, 0.393573654103 , 0.294256140624]
y_zd = [0.391359942215, 0.542088186528, 0.800523068354]
plt.scatter(x_zd,y_zd, marker='*', s=300, color=colors[1:])

plt.xlabel(r'Ellipticity (a.u.)')
plt.ylabel(r'Peak Shift (w)')
plt.ylim(0,0.9)
plt.xlim(0,0.6)
plt.grid()
plt.legend(loc=2, fontsize=12)

### try to plot against inhomogeneity:dephasing ratio
"""
plt.figure()

plt.subplot(211)
plt.scatter(np.ravel(ts_ratio), np.ravel(ellipticity),  
            s=100, alpha=0.7)
plt.xlabel(r'$\mathsf{\tau_{10} / \sigma_{inhom} (w/W)}$')
plt.ylabel(r'Ellipticity')
plt.grid()
plt.subplot(212)
plt.scatter(np.ravel(ts_ratio), np.ravel(fits[...,1]/50.),  
            s=100, alpha=0.7)
plt.xlabel(r'$\mathsf{\tau_{10} / \sigma_{inhom} (w/W)}$')
plt.ylabel(r'Peak shift (w)')
plt.grid()
"""