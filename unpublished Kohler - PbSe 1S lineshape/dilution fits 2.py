# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 17:16:13 2015

simulate the heterodyne equation for demonstrative purposes with 
also:
combine the important part of workup.py (the scaling of integrated intensity)  
    into this figure

@author: Dan
"""
import numpy as np
import matplotlib.pyplot as plt

gamma = 200.
w0 =7625.
# table III of Levenson and Bloembergen
chi_sol = 1.85e-14 # cm^3/erg 

def Delta(nu, gamma, ket_side=True):
    if ket_side:
        out = complex(-1,0)*nu + complex(0,-1)*gamma
        out **= -1
        out_max = np.abs(out).max()
        out /= out_max # so max of this function is 1 (easy for scaling later)
    return -out

def fitfunc(p,x):
    a, gamma_qd = p
    nu_1, c = x
    res = Delta(nu_1,gamma)
    out = a * np.abs(res[None,:] * gamma_qd + c[:,None]**-1 * chi_sol)
    return out
    
linfit = lambda p,x: p[0]*x + p[1]

def errfunc_maker(fit):
    def errfunc(p,y,x):
        """
        robust to 1d or 2d datasets
        """
        fitter = fit(p,x)
        #return (y[1:]-fit[1:]).reshape(1,-1)[0]
        return (y-fitter).reshape(1,-1)[0]
    return errfunc

from scipy.optimize import leastsq

# import the data
data = np.load('nonres_interference.npz')
nu_1 = data['w1']
nu_1 -= w0#data['w1'].mean()
y = data['zis']

x0=[1e30, 1e-30]
od = np.array([0.79, 0.43, 0.18, 0.1, 0.06])
c = np.array(od * 780 / 150e-15 * 10)
p0=leastsq(errfunc_maker(fitfunc), x0, args=(y, [nu_1, c]), full_output=True)
err_tot = errfunc_maker(fitfunc)(p0[0],y,[nu_1,c])
print 'total error', (err_tot**2).sum()
print p0[0][1]
plt.close('all')
plt.figure(figsize=(12,6))
plt.subplot(121)
plt.title('simulation fit method')
sig_guess = fitfunc(p0[0],[nu_1,c])
plt.ylabel(r'$\left|\sigma_{\mathrm{tot}}^{(3)}\right|$ '
            +r'$\mathsf{(a.u.)}$' , 
            fontsize=16)
plt.ylim([0,6])
plt.xlabel(r'$\mathsf{\bar\nu_1 (cm^{-1})}$', fontsize=16)
plt.xticks([6500, 7500, 8500])
colors = ['b', 'g', 'r', 'teal', 'purple']
for i, color  in enumerate(colors):
    if True:#i !=0:
        plt.plot(nu_1+w0, y[i], color, 
                 linewidth=3)
        plt.plot(nu_1+w0, sig_guess[i], color,
                 alpha=0.4, linewidth=10)
plt.grid(b=True)
#"""
plt.subplot(122)
plt.plot(1e32/c**2, (sig_guess**2).sum(axis=1), 'blue',
         linewidth=5, alpha=0.6)
ysum = np.ma.masked_less(y,0)
ysum **= 2
ysum[ysum.mask==True] **= 2
ysum[ysum.mask==True] *= -1
ysum.mask = False
ysum = ysum.sum(axis=-1)
plt.scatter(1e32/c**2, ysum, c='r', s=400, alpha=0.5)	
pline = leastsq(errfunc_maker(linfit), [100,60], args=(ysum, 1e32/c**2))
plt.plot(1e32/c**2, pline[0][0]*1e32/c**2+pline[0][1], 'k--',
         linewidth=5, alpha=0.6)
plt.grid(b=True)
plt.yticks(visible=False)
plt.ylabel(r'$\int{\left|\sigma_{\mathrm{tot}}^{(3)}\right|^2} d\bar\nu_1 \mathsf{(a.u.)}$', fontsize=16)
plt.xlabel(r'$\mathsf{N_{\mathrm{QD}}^{-2} \times 10^{32} (cm^{-3})}$', fontsize=16)
plt.tight_layout()
old_yticks = plt.yticks()[0]
plt.yticks(np.linspace(old_yticks[0], old_yticks[-1], num=7))
print 'gamma_qd = ', p0[0][1], 'cm^6/erg'
print 'n_critical = ', chi_sol / p0[0][1], 'cm^-3'
#"""
"""
res = Delta(nu_1, gamma) #* gamma

sig_tot = np.abs(res[None,:] * gamma_qd  + c[:,None]**-1 * chi_sol)/gamma_qd
plt.close('all')
plt.figure(figsize=(12,6))
plt.subplot(121)
plt.plot(nu_1,sig_tot.T, linewidth=5, alpha=0.5)
plt.grid(b=True)
plt.subplot(122)
plt.plot(1e32*c**-2, (sig_tot**2).sum(axis=1), linewidth=5)
plt.grid(b=True)
#plt.figure()
#plt.plot(nu_1, np.real(res))
"""
plt.savefig('dilution.png', transparent=True, dpi=300)
