# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 21:29:19 2015

calculate the integrated coherence for a resonant excitation
find the trend as we go from driven limit to impulsive limit and between

NOTE:  does not work for only one coherence--needs to be applied to several 
coherences

@author: Dan
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

tau = 200.
mu = 1.
FWHM = np.linspace(-1.6,1.6, num=10)
FWHM = 10**FWHM * tau

start_fine_int = -250
end_fine_int = 500
t1 = np.linspace(-10000, start_fine_int, num=10000)
t2 = np.linspace(start_fine_int+1, end_fine_int-1, num=20000)
t3 = np.linspace(end_fine_int,20000, num=10000)
t = np.hstack([t1,t2,t3])
#t = np.linspace(-5000,10000, num=160000)

def int_rho(t, rho):
    """
    reasonably integrate rho using trapezoid method
    """
    dt = t[1:] - t[:-1]
    out = 0.5 * (rho[1:] + rho[:-1]) * dt
    return np.abs(out).sum()

def CW_limit(mu, tau, pulse_FWHM):
    sigma = pulse_FWHM / (2*np.sqrt(2*np.log(2.)))
    gamma = 1 / tau
    out = mu**3 / gamma**3 * (8*np.sqrt(3)*2*sigma**2 * np.pi)**-1
    return out

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
        #f[i+1] = f[i] + df * dt
    return f

rho1 = [rho(t,tau,mu,normalized_gauss(t,FWHMi), np.ones(len(t))) for FWHMi in FWHM]
rho2 = [rho(t,tau,mu,normalized_gauss(t,FWHMi), rho1[i]) for i,FWHMi in enumerate(FWHM)]
rho3 = [rho(t,tau,mu,normalized_gauss(t,FWHMi), rho2[i]) for i,FWHMi in enumerate(FWHM)]

plt.close('all')
plt.figure()
gs = gridspec.GridSpec(2, 3, wspace=0.05, hspace=0.10)
plt.subplot(gs[0,0])
plt.title('first coherence')
for i,rhoi in enumerate(rho1):
    plt.plot(t,np.abs(rhoi), 
             'r', linewidth=2, alpha=1 - 0.9*i/len(FWHM))
#plt.yticks(visible=False)
plot_max = max([np.abs(arr).max() for arr in rho1])
plt.ylim(-0.1*plot_max, 1.1*plot_max)
plt.grid()
plt.xlim(-100,1200)

plt.subplot(gs[0,1])
plt.title('second coherence')
for i,rhoi in enumerate(rho2):
    plt.plot(t,np.abs(rhoi), 
             'r', linewidth=2, alpha=1 - 0.9*i/len(FWHM))
plt.grid()
#plt.yticks(visible=False)
plt.xlim(-100,1200)
plot_max = max([np.abs(arr).max() for arr in rho2])
plt.ylim(-0.1*plot_max, 1.1*plot_max)

plt.subplot(gs[0,2])
plt.title('third coherence')
for i,rhoi in enumerate(rho3):
    plt.plot(t,np.abs(rhoi), 
             'r', linewidth=2, alpha=1 - 0.9*i/len(FWHM))
plt.grid()
#plt.yticks(visible=False)
plt.xlim(-100,1200)
plot_max = max([np.abs(arr).max() for arr in rho3])
plt.ylim(-0.1*plot_max, 1.1*plot_max)


y1 = [int_rho(t,rhoi) for rhoi in rho1]
y2 = [int_rho(t,rhoi) for rhoi in rho2]
y3 = [int_rho(t,rhoi) for rhoi in rho3]

plt.subplot(gs[1,0])
plt.scatter(FWHM/tau,y1/max(y1),
            s=50, alpha=0.8)
plt.xlim(0.9*FWHM.min() / tau, 1.1*FWHM.max()/tau)
plt.ylim(0.0005, 2)
plt.yscale('log')
plt.xscale('log')
#plt.xlabel(r'$\mathsf{\frac{FWHM}{\tau}}$', fontsize = 18)
#plt.tight_layout()
plt.grid()

plt.subplot(gs[1,1])
plt.scatter(FWHM/tau,y2/max(y2),
            s=50, alpha=0.8)
plt.xlim(0.9*FWHM.min() / tau, 1.1*FWHM.max()/tau)
plt.ylim(0.0005, 2)
plt.yticks(visible=False)
plt.yscale('log')
plt.xscale('log')
#plt.xlabel(r'$\mathsf{\frac{FWHM}{\tau}}$', fontsize = 18)
#plt.tight_layout()
plt.grid()

plt.subplot(gs[1,2])
plt.scatter(FWHM/tau,y3/max(y3),
            s=50, alpha=0.8)
ytemp = CW_limit(mu, tau, FWHM) / max(y3)
plt.plot(FWHM/tau, ytemp, 
         color='g', linewidth=3, alpha=0.6)
plt.plot(FWHM/tau, np.ones(len(FWHM)) * tau * 0.167 / (8* max(y3))  , 
         color='r', linewidth=3, alpha=0.6)
plt.xlim(0.9*FWHM.min() / tau, 1.1*FWHM.max()/tau)
plt.ylim(0.0005, 2.)
#plt.ylim(-0.1, 1.1)
plt.yticks(visible=False)
plt.yscale('log')
plt.xscale('log')
#plt.xlabel(r'$\mathsf{\frac{FWHM}{\tau}}$', fontsize = 18)
#plt.tight_layout()
plt.grid()

plt.savefig('single coherence.png')

### an important detail in predicting the impulsive limite yield:
#       compute the integral of the impulsive limit for peak height:
#       non-trivial!
"""
tprime = np.linspace(-600, 600, num=5000)
dt = tprime[1]-tprime[0]
gauss1 = normalized_gauss(tprime, 200.)
gauss2 = np.cumsum(gauss1) * dt * gauss1
gauss3 = np.cumsum(gauss2) * dt * gauss1
plt.figure()
plt.plot(tprime, np.cumsum(gauss1)*dt)
plt.plot(tprime, np.cumsum(gauss2)*dt)
plt.plot(tprime, np.cumsum(gauss3)*dt)

print gauss1.sum() * dt
print gauss2.sum() * dt
print gauss3.sum() * dt
#"""