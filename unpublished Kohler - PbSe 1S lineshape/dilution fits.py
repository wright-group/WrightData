# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 19:30:53 2014

shortcut clips for processing data from the dilution study

@author: Dan
"""

import fscolors as f
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
figsize=(12,6)
#--------------------Abs spectra--------------------
#--------------------WIGNERS--------------------
xsmooth=0

abs_file = os.path.expanduser('~')
abs_file = '\\'.join([abs_file,
           'Documents', 
           'Colloidal Synthesis Characterization',
           'DK-2014.02.05-PbSe 019'])
abs_file += r'\fs19 2014.09.01.trimmed.txt'
l,a = np.loadtxt(abs_file, unpack = True, skiprows=18)
a -= a.min()
A = interp1d(1e7/l, a)

# import the opa power smoothness curves
opa_folder = os.path.expanduser('~')
opa_folder = '\\'.join([opa_folder,
                          'Desktop',
                          'incoming from fs table',
                          'fsb 19',
                          '2014.09.14',
                          'calibrations'])
opa1_file = opa_folder + r'\opa1 power smoothness (1) 4000 ms wait [1x51] Freq.dat'
#opa2_file = opa_folder + r'\opa2 power smoothness (1) 4000 ms wait [1x51] Freq.dat'
opa1_dat = np.loadtxt(opa1_file)
#opa2_dat = np.loadtxt(opa2_file)

f_opa1 = interp1d(1e7/opa1_dat[:,1], opa1_dat[:,16])
#f_opa2 = interp1d(opa1_dat[:,3], opa1_dat[:,16])

I_1 = f_opa1(np.linspace(6256,8699,num=51))

files_folder = os.path.expanduser('~')
files_folder = '\\'.join([files_folder,
                          'Desktop',
                          'incoming from fs table',
                          'fsb 19',
                          '2014.09.19 dilution study'])
filename = r'\fsb 19.1 wigner (2) 60 ms wait [56x51] Wigner.dat'
filename = files_folder + filename
darkest_sample_path = os.path.expanduser('~')
darkest_sample_path = '\\'.join([darkest_sample_path,
                                 'Desktop',
                                 'incoming from fs table',
                                 'fsb 19',
                                 '2014.09.14'])
darkest_sample_path += r'\w2_d2 high res 2014.09.09 (12) 500 ms wait [29x51] Wigner.dat'

fs = []
od = np.array([0.79, 0.43, 0.18, 0.1, 0.06])
c = od * 780 / 150e-15 * 10
# get number density from this:
# concentration ratios as measured from gaussian peak fits (select amplitude)
n = np.array([0.625, 0.347, 0.147, 0.0764, 0.0475])
n = n / n[0]

names = [r'\fsb 19.1 wigner (2) 60 ms wait [56x51] Wigner.dat',
         r'\fsb 19.2 wigner (1) 60 ms wait [56x51] Wigner.dat',
         r'\fsb 19.3 wigner (1) 60 ms wait [56x51] Wigner.healed.dat',
         r'\fsb 19.4 wigner (2) 60 ms wait [56x49] Wigner.healed.dat'
         #r'\.dat',
         ]

def M(w, ax):
    """
     1D M-factor
     ax controls how to scale A
    """
    a2 = A(7500)*ax
    a1 = A(w)*ax
    m = 10**(-a1/2.) * (1- 10**(-a2))
    m /= a2 * np.log(10)
    return m

if True:
    fi = f.Dat(filepath=darkest_sample_path, 
               xvar='w1', yvar='d2', gridfactor=1)
    fi.level(-4)
    fi.smooth(x=xsmooth)
    fi.zi /= I_1
    zii = fi.zi.copy()
    zii = np.ma.masked_less(zii,0)
    zii **= 0.5
    zii.mask = False
    zii = np.ma.masked_greater(zii,0)
    zii *= -1
    zii **= 0.5
    zii *= -1
    zii.mask = False    
    fi.zi_norm = zii
    fi.zfit()
    #plt.plot2d(alt_zi='amp')
    #plt.close()
    fs.append(fi)

    for name in names:
        fi = f.Dat(filepath=files_folder + name, 
                   xvar='w1', yvar='d2', gridfactor=1)
        fi.level(-3)
        fi.smooth(x=xsmooth)
        if fi.zi.shape[-1]==51:
            fi.zi /= I_1
        else:
            fi.zi /= I_1[1:-1]
        zii = fi.zi.copy()
        zii = np.ma.masked_less(zii,0)
        zii **= 0.5
        zii.mask = False
        zii = np.ma.masked_greater(zii,0)
        zii *= -1
        zii **= 0.5
        zii *= -1
        zii.mask = False    
        fi.zi_norm = zii
        fi.zfit()
        #fi.plot2d(alt_zi='amp')
        plt.close()
        fs.append(fi)

num=2
plt.figure(figsize=figsize)
plt.subplot(121)
z0m = n[0]/1.6
plt.plot(fs[0].xi, fs[0].zi_norm[2]/z0m, 
         linewidth=5, 
         alpha=0.6, 
         label=r'$\mathrm{OD_{1S}} = $' + '{0}'.format(round(od[0],2)))

for i in range(len(fs))[1:]:
    fsi = fs[i]
    zim = n[i]#fsi.zi_norm[num].max()
    zii = fsi.zi_norm[num]/zim
    plt.plot(fsi.xi, zii, 
             linewidth=5, 
             alpha=0.6, 
             label=r'{0}'.format(round(od[i],2)))
plt.ylabel(r'$\sigma^{(3)}$ '+r'$\mathsf{(a.u.)}$' , 
            fontsize=16)
plt.xlabel(r'$\mathsf{\bar\nu_1 (cm^{-1})}$', fontsize=16)
plt.xticks([6500, 7500, 8500])
plt.title(r'raw')
# overlay absorbance
plt.plot(fs[0].xi, A(fs[0].xi)*5, 'k', linewidth=10, alpha=0.3)
plt.grid(b=True)
plt.ylim(0,4.5)
l1 = plt.legend(loc=0, fontsize=12)
plt.title('raw')


plt.subplot(122)

y_corrected = fs[0].zi_norm[2]/z0m * 1/M(fs[0].xi, n[0])
plt.plot(fs[0].xi, y_corrected, linewidth=5, alpha=0.6)

for i in range(len(fs))[1:]:
    y_corrected = fs[i].zi_norm[2] * 1/(M(fs[i].xi, n[i])*n[i])
    plt.plot(fs[i].xi, y_corrected, linewidth=5, alpha=0.6)
plt.grid()
plt.xlabel(r'$\mathsf{\bar\nu_1 (cm^{-1})}$', fontsize=16)
plt.xticks([6500, 7500, 8500])
plt.xlim([fs[0].xi.min(), fs[0].xi.max()])
plt.ylim(0,4.5)
plt.yticks(visible=False)
plt.title('corrected')

#plt.savefig('dilution1.png', transparent=True, dpi=300)
if True:    
    num=43
    plt.figure(figsize=figsize)
    z0m = n[0]/1.6 #fs[0].zi_norm[2].max()
    """
    plt.subplot(121)
    #plt.title(r'$\tau_{21} = 30$ ' + 'fs')
    plt.plot(fs[0].xi,0.5*(fs[0].zi_norm[14]+fs[0].zi_norm[15])/z0m, 
             linewidth=5, 
             alpha=0.6, 
             label=r'$\mathrm{OD_{1S}} = $' + '{0}'.format(round(od[0],2)))
    for i in range(len(fs))[1:]:
        fsi = fs[i]
        zim = n[i]#fsi.zi_norm[num].max()
        plt.plot(fsi.xi,fsi.zi_norm[num]/zim, 
                 linewidth=5, 
                 alpha=0.6, 
                 label=r'{0}'.format(round(od[i],2)))
    plt.ylabel(r'$\sigma^{(3)}$ '+r'$\mathsf{(a.u.)}$' , 
                fontsize=16)
    #plt.ylabel(r'$\sqrt{S_{\mathsf{TrEE}}^{(0)}} / N$ '+r'$\mathsf{(a.u.)}$' , 
    #            fontsize=16)
    plt.xlabel(r'$\mathsf{\bar\nu_1 (cm^{-1})}$', fontsize=16)
    plt.xticks([6500, 7500, 8500])
    plt.ylim([0,6])
    plt.grid()
    plt.legend(loc=0, fontsize=12)
    """
    
    plt.subplot(121)
    plt.ylabel(r'$\left|\sigma_{\mathrm{tot}}^{(3)}\right|$ '
                +r'$\mathsf{(a.u.)}$' , 
                fontsize=16)
    # save this data to a dictionary of arrays
    non_res_data = np.zeros((len(fs),fs[0].xi.size-2))
    corrected0 = 0.5*(fs[0].zi_norm[14]+fs[0].zi_norm[15])/z0m* 1/M(fs[0].xi, n[0]).copy()
    non_res_data[0] = corrected0[1:-1] 
    plt.plot(fs[0].xi,corrected0, 
             linewidth=5, 
             alpha=0.6, 
             label=r'$\mathrm{OD_{1S}} = $' + '{0}'.format(round(od[0],2)))
    for i in range(len(fs))[1:]:
        fsi = fs[i]
        zim = n[i]#fsi.zi_norm[num].max()
        correctedi = fsi.zi_norm[num]/zim * 1/M(fs[i].xi, n[i]).copy()
        stri = 'z{0}'.format(i)
        if correctedi.size == 51:
            non_res_data[i] = correctedi[1:-1]
        else:
            non_res_data[i] = correctedi
        plt.plot(fsi.xi, correctedi, 
                 linewidth=5, 
                 alpha=0.6, 
                 label=r'{0}'.format(round(od[i],2)))
    plt.yticks(visible=False)
    plt.xlabel(r'$\mathsf{\bar\nu_1 (cm^{-1})}$', fontsize=16)
    plt.legend(loc=0, fontsize=12)
    plt.xticks([6500, 7500, 8500])
    plt.ylim([0,6])
    plt.grid()

    plt.subplot(122)
    # for non_res_data, I want the intensity-level data!
    # if it's negative, the square should also be negative
    #y = (non_res_data**2.00001).sum(axis=-1)
    cinv2 = 1/(c**2) * 1e32
    y = np.ma.masked_less(non_res_data,0)
    y **= 2
    y[y.mask==True] *= -1
    y[y.mask==True] **= 2
    y[y.mask==True] *= -1
    y.mask = False
    y = y.sum(axis=-1)
    y *= fs[0].xi[1]-fs[0].xi[0]
    #y /= 49
    plt.scatter(cinv2, y, linewidth=5, alpha=0.6)
    plt.grid()
    from scipy.optimize import leastsq
    fitfunc = lambda p,x: p[0]*x + p[1]
    errfunc = lambda p,x,y: y - fitfunc(p,x)
    fit1 = leastsq(errfunc, 
                   [50,50], 
                   args=(cinv2,y))
    print fit1
    plt.plot(np.linspace(cinv2.max(), cinv2.min()), 
             fitfunc(fit1[0],
                     np.linspace(cinv2.max(), cinv2.min())), 
             'r', linewidth=5, alpha=0.6)    
    plt.yticks(visible=False)
    target_val = fit1[0][1]/800 * 2400 + fit1[0][1]
    plt.plot([cinv2.max(), cinv2.min()],[target_val,target_val],'k:', linewidth=3)    
    plt.ylabel(r'$\int{\left|\sigma_{\mathrm{tot}}^{(3)}\right|^2} d\bar\nu_1 \mathsf{(a.u.)}$', fontsize=16)
    plt.xlabel(r'$\mathsf{N_{\mathrm{QD}}^{-2} \times 10^{32} (cm^{-3})}$', fontsize=16)
    #plt.ylim([0,6])
    #plt.tight_layout()

# take those params from the fit and calculate an approximate non-linear 
# susceptibility
chi_sol = (fit1[0][0] / 2347) # intensity scale
FWHM = 400 # parameter to adjust
chi_qd_peak = (fit1[0][1] / FWHM) # intensity scale
print 'chi_sol (average):  ', chi_sol
print 'chi_qd (peak):  ', chi_qd_peak
intd_chi_qd = chi_qd_peak*2347 
n_critical = intd_chi_qd / fit1[0][0]
n_critical **= -0.5
n_critical *= 1e16
print 'n* = ', n_critical
sigma_qd = n_critical**-1 * 1.85e-14
print 'sigma_qd (peak):  ', sigma_qd
# save the datasets for some modelling
np.savez('nonres_interference.npz',
         n = n, w1=fs[0].xi[1:-1], zis=non_res_data)  

#if False:
#    # to plot signals across different concentrations
#    import matplotlib.pyplot as plt
#    plt.plot(f1.yi, f1.zi.sum(axis=1)/.43, label='1')
#    plt.plot(f2.yi, f2.zi.sum(axis=1)/.18, label='2')
#    plt.plot(f3.yi, f3.zi.sum(axis=1)/.1, label='3')
#    plt.plot(f4.yi, f4.zi.sum(axis=1)/.06, label='4')
#    plt.xlabel(r'$\mathsf{\tau_{21}}$ (fs)')
#    plt.grid()
#    plt.legend()

