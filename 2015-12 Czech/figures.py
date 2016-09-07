import os
import copy

import numpy as np
from scipy import signal

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import matplotlib.colors as mplcolors
import matplotlib.patheffects as PathEffects
from mpl_toolkits.axes_grid1 import make_axes_locatable
matplotlib.rcParams['font.size'] = 18
plt.ioff()

import fscolors as f  

import WrightTools as wt
box_path = wt.kit.get_box_path()
tuning_path = os.path.join(box_path, 'MX2', '2015.02.08', 'Post-optimization data', '2015.01.29 OPA tuning')
wright_cm = wt.artists.colormaps['wright']

fluencies = np.array([653, 359, 285, 227, 90])

def gauss(x, p):
    # p = amp, mu, sigma
    a, mu, sigma = p
    return a*np.exp(-(x-mu)**2 / (2*np.abs(sigma)**2))

a_params = [1.81054482E-02, 1.44717631E04, 3.03607786E02]
b_params = [1.89676390E-02, 1.57011616E04, 4.38342016E02]

A_energy = a_params[1] #wn
B_energy = b_params[1] #wn
energies = [14000, 14500, 15100, 15700, 16200]
A_eV = wt.units.converter(A_energy, 'wn', 'eV')
B_eV = wt.units.converter(B_energy, 'wn', 'eV')

data_pickle_path = os.path.join(box_path, 'MX2', '2015.02.08', 'Post-optimization data', 'Fluence Study', 'data.p')
data = wt.data.from_pickle(data_pickle_path)
data.convert('eV')
data_at_fluence = data.chop('w2', 'w1', 'd2', {'fluence': [227, 'uJ per sq. cm']})[0]
data_at_fluence.level(0, 'd2', -3)
data_at_fluence.smooth([2, 2, 0])
data_at_fluence.normalize()
data_at_fluence.scale()
data_min = data_at_fluence.channels[0].znull
data_max = data_at_fluence.channels[0].zmax

abs_data_path = os.path.join(box_path, 'MX2', 'outside spectra', 'thin film abs.p')
abs_data = wt.data.from_pickle(abs_data_path)
abs_data.convert('eV')

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    #subax.xaxis.set_tick_params(labelsize=x_labelsize)
    #subax.yaxis.set_tick_params(labelsize=y_labelsize)
    plt.setp(subax.get_xticklabels(), visible=False)
    plt.setp(subax.get_yticklabels(), visible=False)
    for tic in subax.xaxis.get_major_ticks(): tic.tick1On = tic.tick2On = False
    for tic in subax.yaxis.get_major_ticks(): tic.tick1On = tic.tick2On = False
    return subax

### intro / overview figures ##################################################

if 0:
    #AFM playground
    plt.close('all')
    afm_path = os.path.join(box_path, 'MX2', 'papers', 'CMDS of MoS2 thin films', 'data', 'AFM.txt')
    AFM_data = np.genfromtxt(afm_path, skiprows = 163)*1000
    AFM_data.shape = (256, 256)
    xi = np.linspace(0, 256*.215, 256)
    yi = np.linspace(0, 256*.215, 256)
    lower_cutoff = 80
    upper_cutoff = 220
    AFM_data = AFM_data[lower_cutoff:upper_cutoff, lower_cutoff:upper_cutoff]
    xi = xi[lower_cutoff:upper_cutoff]
    yi = yi[lower_cutoff:upper_cutoff]
    
    i = 50
    a = np.ones((i , i ))    
    a /= a.sum() 
    AFM_slow = signal.convolve2d(AFM_data, a, mode = 'valid')
    AFM_data = AFM_data[i/2:(-i/2)+1, i/2:(-i/2)+1]
    xi = xi[i/2:(-i/2)+1]    
    yi = yi[i/2:(-i/2)+1]
    AFM_data = AFM_data - AFM_slow
    
    '''
    i = 3
    a = np.ones((i , i ))    
    a /= a.sum() 
    AFM_data = signal.convolve2d(AFM_data, a, mode = 'valid')
    xi = xi[i/2:(-i/2)+1]    
    yi = yi[i/2:(-i/2)+1]    
    '''
    
    levels = np.linspace(AFM_data.min(), AFM_data.max()/3, 200)
    cmap = plt.get_cmap('afmhot')
    cmap.set_over(color = 'r')
    fig = plt.figure()
    ax = plt.subplot()
    plt.contourf(xi, yi, AFM_data, levels, cmap = cmap, clip = False)
    plt.colorbar()
    plt.xlim(25, 40)
    plt.ylim(25, 40)
    plt.grid()
    plt.axhline(yi[82])

    plt.figure()
    plt.plot(xi, AFM_data[82])
    #plt.show()
    plt.savefig('AFM.tif', transparent = True, dpi = 300)


if False:
    #material introduction

    plt.figure(figsize = (10, 15))

    ax1 = plt.subplot2grid((3, 2), (0, 0))#, colspan = 2)
    ax2 = plt.subplot2grid((3, 2), (0, 1))
    ax3 = plt.subplot2grid((3, 2), (1, 0))
    ax4 = plt.subplot2grid((3, 2), (1, 1))
    ax5 = plt.subplot2grid((3, 2), (2, 0), colspan = 2)
    
    #images
    '''
    image1_path = os.path.join(box_path, 'MX2', 'papers', 'CMDS of MoS2 thin films', 'data', 'images of thin film', 'Picture1.png')
    image2_path = os.path.join(box_path, 'MX2', 'papers', 'CMDS of MoS2 thin films', 'data', 'images of thin film', 'Picture2.png')
    image1 = matplotlib.image.imread(image1_path)
    image2 = matplotlib.image.imread(image2_path)
    image1 = image1[0:200, 0:290]
    spacer = np.zeros((6, 290, 4))
    image2 = image2[0:200, 0:290]
    image = np.append(image1, spacer, axis = 0)
    image = np.append(image, image2, axis = 0)
    npad = ((12,12),(70,70),(0,0))
    image = np.pad(image, pad_width=npad, mode='constant', constant_values=0)
    ax1.imshow(image)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), visible=False)
    '''
    

    image1_path = os.path.join(box_path, 'MX2', 'papers', 'CMDS of MoS2 thin films', 'data', 'images of thin film', 'Picture1.png')
    image2_path = os.path.join(box_path, 'MX2', 'papers', 'CMDS of MoS2 thin films', 'data', 'images of thin film', 'Picture2.png')
    image1 = matplotlib.image.imread(image1_path)
    image2 = matplotlib.image.imread(image2_path)
    image1 = image1[0:200, 0:290]
    spacer = np.zeros((6, 290, 4))
    image2 = image2[0:200, 0:290]
    image = np.append(image1, spacer, axis = 0)
    image = np.append(image, image2, axis = 0)
    npad = ((12,12),(70,70),(0,0))
    image = np.pad(image, pad_width=npad, mode='constant', constant_values=0)
    
    npad = ((45+12,45+12),(12,12),(0,0))
    image1 = np.pad(image1, pad_width=npad, mode='constant', constant_values=0)   
    ax1.imshow(image1)
    ax1.text(150, 260, '1.5 cm', ha = 'center', va = 'top', size = 20)
    ax1.annotate("", xy=(35, 250), xycoords='data', xytext=(270, 235), textcoords='data',
                 arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", lw = 1.5))
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), visible=False)
    for tic in ax1.xaxis.get_major_ticks(): tic.tick1On = tic.tick2On = False
    for tic in ax1.yaxis.get_major_ticks(): tic.tick1On = tic.tick2On = False

    npad = ((45+12,45+12),(12,12),(0,0))
    image2 = np.pad(image2, pad_width=npad, mode='constant', constant_values=0)   
    ax2.imshow(image2)
    ax2.text(150, 260, '1.5 cm', ha = 'center', va = 'top', size = 20)
    ax2.annotate("", xy=(30, 250), xycoords='data', xytext=(295, 250), textcoords='data',
                 arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", lw = 1.5))
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    for tic in ax2.xaxis.get_major_ticks(): tic.tick1On = tic.tick2On = False
    for tic in ax2.yaxis.get_major_ticks(): tic.tick1On = tic.tick2On = False

    #raman
    raman_data_path = os.path.join(box_path, 'MX2', 'papers', 'CMDS of MoS2 thin films', 'data', 'Raman_Spectroscopy_Sample_1.3.3.3_1.csv')
    raman_data = np.genfromtxt(raman_data_path, delimiter = ',').T
    xi = raman_data[0]
    zi = raman_data[1]
    ax3.plot(xi, zi, lw = 1.5)
    ax3.axvline(382.4, color = 'k', alpha = 1, lw = 1, ls = '--')
    ax3.axvline(407.5, color = 'k', alpha = 1, lw = 1, ls = '--')
    ax3.set_xlim(365, 425)
    ax3.set_ylim(0, 5500)
    '''
    ax3.text(395, 4600, '$\mathsf{\Delta = 25.1}$', ha = 'center', va = 'bottom', size = 20)
    ax3.annotate("", xy=(382.4, 4500), xycoords='data', xytext=(407.5, 4500), textcoords='data',
                 arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", lw = 1.5))
    '''
    ax3.set_ylabel('Raman (counts)')
    ax3.set_xlabel('Raman Shift $\mathsf{(cm^{-1})}$')

    #TEM
    TEM_image_path = os.path.join(box_path, 'MX2', 'papers', 'CMDS of MoS2 thin films', 'data', 'TEM', 'TEM.png')
    diffraction_image_path = os.path.join(box_path, 'MX2', 'papers', 'CMDS of MoS2 thin films', 'data', 'TEM', 'Diffraction.png')
    image1 = matplotlib.image.imread(TEM_image_path)
    image2 = matplotlib.image.imread(diffraction_image_path)
    ax4.imshow(image1)
    rect = [0.55,0.55,0.4,0.4]
    ax4_sub = add_subplot_axes(ax4, rect)
    ax4_sub.imshow(image2)
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), visible=False)
    plt.setp(ax4_sub.get_xticklabels(), visible=False)
    plt.setp(ax4_sub.get_yticklabels(), visible=False)  
   
    #absorbance
    abs_xi = abs_data.axes[0].points
    abs_zi = abs_data.channels[0].values
    ax5.plot(abs_xi, abs_zi, lw = 1.5)    
    ax5.axvline(A_eV, color = 'k', alpha = 1, lw = 1, ls = '--')
    ax5.axvline(B_eV, color = 'k', alpha = 1, lw = 1, ls = '--')
    ax5.set_ylim(0, 0.065)
    ax5.set_ylabel('Absorbance (OD)', color = 'b')
    ax5.set_xlabel('$\mathsf{\hslash \omega \,\, (eV)}$')
    for tl in ax5.get_yticklabels(): tl.set_color('b')
    
    #simulated
    abs_xi_wn = wt.units.converter(abs_xi, 'eV', 'wn')
    a_zi = gauss(abs_xi_wn, a_params)
    b_zi = gauss(abs_xi_wn, b_params)
    remainder = abs_zi-(a_zi+b_zi)
    ax5.plot(abs_xi, a_zi, 'k--')
    ax5.plot(abs_xi, b_zi, 'k--')
    ax5.plot(abs_xi, remainder, 'k--')
    #ax5.plot(abs_xi, a_zi+b_zi, 'k:')
    
    #PL
    PL_file = r'C:\Users\Blaise\Box Sync\Wright Shared\MX2\outside spectra\PL\MoS2_ThinFilm_III__532nm__1.0mW_2.CSV'
    PL_path = os.path.join(box_path, 'MX2', 'outside spectra', 'PL', 'MoS2_ThinFilm_III__532nm__1.0mW_2.CSV')
    PL_data = np.genfromtxt(PL_file, delimiter = ',').T
    pl_xi = PL_data[0]
    pl_xi = wt.units.converter(pl_xi, 'wn', 'eV')
    pl_zi = PL_data[1]
    ax6 = ax5.twinx()
    ax6.plot(pl_xi, pl_zi, lw = 1.5, color ='g')
    ax6.set_ylabel('Photoluminescence (counts)', color = 'g')
    for tl in ax6.get_yticklabels(): tl.set_color('g')
    ax6.set_ylim(15, 150)
    
    #OPA
    OPA2_dat = os.path.join(tuning_path, 'OPA2', 'OPA2 Mixer2 NON-SH-NON-Sig OPA2_SHS_tunetest_2015.01.27 VAI0 0 ms wait [41x15] Pts Test.dat')
    dat = f.Dat(filepath = OPA2_dat, xvar = 'w2', yvar = 'wm', cols = 'v2', colortune = True)
    OPA2_xi = wt.units.converter(dat.xi, 'wn', 'eV')
    OPA2_yi = dat.yi
    OPA2_zi = dat.zi
    OPA2_zi -= OPA2_zi.min()
    OPA2_zi /= OPA2_zi.max()
    slice_index = 9
    xi = wt.units.converter(OPA2_yi, 'wn', 'eV') + OPA2_xi[slice_index]
    zi = OPA2_zi[:, slice_index] * 0.06
    ax5.plot(xi, zi, color = 'r', lw = 1.5)
    
    ax5.set_xlim(1.65, 2.15)
    
    #labels
    props = dict(boxstyle='square', facecolor='white', alpha=0.75)
    ax1.text(0.05, 0.95, 'a', transform=ax1.transAxes, fontsize=20, verticalalignment='top', horizontalalignment='left', bbox=props)
    ax2.text(0.05, 0.95, 'b', transform=ax2.transAxes, fontsize=20, verticalalignment='top', horizontalalignment='left', bbox=props)    
    ax3.text(0.05, 0.95, 'c', transform=ax3.transAxes, fontsize=20, verticalalignment='top', horizontalalignment='left', bbox=props)
    ax4.text(0.05, 0.95, 'd', transform=ax4.transAxes, fontsize=20, verticalalignment='top', horizontalalignment='left', bbox=props)
    ax5.text(0.025, 0.95, 'e', transform=ax5.transAxes, fontsize=20, verticalalignment='top', horizontalalignment='left', bbox=props)
        
    #plt.tight_layout()        
        
    plt.savefig('material overview.tif', transparent = True, dpi = 300)
    plt.close('all')
    
if False:
    #spectroscopy introduction
    si_alpha = 0.25

    plt.figure(figsize = (10, 10))

    #gs = grd.GridSpec(3, 3)

    ax1 = plt.subplot2grid((2, 2), (0, 0))
    ax2 = plt.subplot2grid((2, 2), (0, 1))
    ax3 = plt.subplot2grid((2, 2), (1, 0))
    ax4 = plt.subplot2grid((2, 2), (1, 1))    
    #ax1_subA = add_subplot_axes(ax1, [0, 0.55, 1, 0.45])
    #ax1_subB = add_subplot_axes(ax1, [0, 0, 1, 0.45])
    '''
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), visible=False)
    for tic in ax1.xaxis.get_major_ticks(): tic.tick1On = tic.tick2On = False
    for tic in ax1.yaxis.get_major_ticks(): tic.tick1On = tic.tick2On = False
    for child in ax1.get_children():
        if isinstance(child, matplotlib.spines.Spine):
            child.set_color('r')
            child.set_linewidth(0)
    '''
    #mask
    '''
    scalar = 9.
    mp = [[-0.0874*scalar, 0*scalar], #w1
          [0.0175*scalar, 0*scalar], #w2'
          [0*scalar, 0.0262*scalar]] #w2    
    label_offset = 0
    mec = matplotlib.colors.ColorConverter().to_rgba('k', alpha=1)
    #    
    ax1_subA.annotate('$\mathsf{1}$', [mp[0][0], mp[0][1]], color='k', size = 20, va = 'center', ha = 'center')
    ax1_subA.plot(mp[0][0], mp[0][1], color='r', marker='o', lw = 0, markeredgecolor = mec, markersize = 30, alpha = si_alpha)
    #
    ax1_subA.annotate('$\mathsf{2}$', [mp[2][0], mp[2][1]], color='k', size = 20, va = 'center', ha = 'center')
    ax1_subA.plot(mp[2][0], mp[2][1], color='b', marker='o', lw = 0, markeredgecolor = 'k', markersize = 30, alpha = si_alpha)
    #    
    ax1_subA.annotate('$\mathsf{2^{\prime}}$', [mp[1][0], mp[1][1]], color='k', size = 20, va = 'center', ha = 'center')
    ax1_subA.plot(mp[1][0], mp[1][1], color='b', marker='o', lw = 0, markeredgecolor = None, markersize = 30, alpha = si_alpha)
    #
    out_x = mp[0][0] + mp[1][0] - mp[2][0]
    out_y = mp[0][1] + mp[1][1] - mp[2][1]
    ax1_subA.plot(out_x, out_y, color='r', marker='*', lw = 0, markeredgecolor = None, markersize = 30, alpha = si_alpha)
    #
    ax1_subA.plot(0, 0, color='k', marker='o', lw = 0, markeredgecolor = None, markersize = 5)
    #    
    aspect = 0.45
    ax1_subA.set_xlim(-1.25, 0.75) 
    ax1_subA.set_ylim(-1*aspect, 1*aspect) 
    '''
   
    #pulse delay space
    N = 1000
    T = 1.0 / 50.0
    x = np.linspace(-N*T, N*T, N)
    tau = np.linspace(-0.5, 0.5, 1000)
    sig = 2
    #
    def pulse(t):
        gaus = np.exp(-((t-0.0)**2)/((sig)**2))
        return gaus
    #
    #plt.subplot(ax1_subA)
    ax1.set_ylim(-0.2, 1.3)
    ax1.set_xlim(-23, 23)
    ax1.plot(x, pulse(x-12), 'k', alpha = si_alpha)
    ax1.fill_between(x, 0, pulse(x-12), facecolor='b', alpha = si_alpha)
    ax1.plot(x, pulse(x), 'k', alpha = si_alpha)
    ax1.fill_between(x, 0, pulse(x), facecolor='b' , alpha = si_alpha)
    ax1.plot(x, pulse(x+12), 'k', alpha = si_alpha)
    ax1.fill_between(x, 0, pulse(x+12), facecolor='r', alpha = si_alpha)
    ax1.annotate("", xy=(-12.0, 1.06), xytext=(0.0, 1.06), arrowprops=dict(arrowstyle="->"))
    ax1.annotate("", xy=(0.0, 1.06), xytext=(12.0, 1.06), arrowprops=dict(arrowstyle="<-"))
    ax1.text(-6.0, 1.08, r'$\mathsf{\tau_{21}}$', fontsize=25, ha = 'center', va = 'bottom')
    ax1.text(6.0, 1.08, r'$\mathsf{\tau_{22^{\prime}}}$', fontsize=25, ha = 'center', va = 'bottom')
    ax1.annotate("", xy=(-20.0, 0.0), xytext=(21.0, 0.0), arrowprops=dict(arrowstyle="<|-", lw=2.0, facecolor='black'), size=40.0)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), visible=False)
    ax1.text(0, -0.05, 'direction of travel', fontsize=18, ha = 'center', va = 'top')
    for tic in ax1.xaxis.get_major_ticks(): tic.tick1On = tic.tick2On = False
    for tic in ax1.yaxis.get_major_ticks(): tic.tick1On = tic.tick2On = False
    #
    ax1.annotate('$\mathsf{1}$', [-12, 0.2], color='k', size = 30, va = 'center', ha = 'center')
    ax1.annotate('$\mathsf{2}$', [0, 0.2], color='k', size = 30, va = 'center', ha = 'center')
    ax1.annotate('$\mathsf{2^{\prime}}$', [12, 0.2], color='k', size = 30, va = 'center', ha = 'center')
    
    #2D delay
    fake_2D_path = os.path.join(box_path, 'MX2', 'papers', 'CMDS of MoS2 thin films', 'figures', 'fake 2D delay', 'data.p')
    data = wt.data.from_pickle(fake_2D_path)
    data.zoom(2)
    data.chop('d0', 'd2')
    data.normalize()
    xi = data.axes[0].points
    yi = data.axes[1].points
    zi = data.channels[0].values.T
    levels = np.linspace(0, 1, 200)
    ax2.contourf(xi, yi, zi, levels, cmap = wright_cm)
    ax2.set_xlim(-175, 175)
    ax2.set_ylim(-175, 175)
    ax2.grid()
    diag_min = max(min(xi), min(yi))
    diag_max = min(max(xi), max(yi))
    ax2.plot([diag_min, diag_max],[diag_min, diag_max], 'k-')
    ax2.axhline(0, color = 'k')
    ax2.axvline(0, color = 'k')
    plt.subplot(ax3)
    plt.xticks(rotation = 45)
    ax2.set_ylabel(r'$\mathsf{\tau_{21}}$', fontsize = 25)
    ax2.set_xlabel(r'$\mathsf{\tau_{22^{\prime}}}$', fontsize = 25)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    #ax2.text(0.95, 0.05, '$\mathsf{\hslash \omega_1 \, = \, \hslash \omega_2}$', transform=ax2.transAxes, fontsize=20, verticalalignment='bottom', horizontalalignment='right', bbox=props)    
    factor = 200  
    ax2.text(-0.5*factor, 0.5*factor, 'I', fontsize = 35, verticalalignment = 'center', horizontalalignment = 'center', path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")])
    ax2.text(0.25*factor, 0.6*factor, 'II', fontsize = 35, verticalalignment = 'center', horizontalalignment = 'center', path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")])
    ax2.text(-0.6*factor, -0.25*factor, 'III', fontsize = 35, verticalalignment = 'center', horizontalalignment = 'center', path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")])
    ax2.text(0.6*factor, 0.25*factor, 'IV', fontsize = 35, verticalalignment = 'center', horizontalalignment = 'center', path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")])
    ax2.text(-0.25*factor, -0.6*factor, 'V', fontsize = 35, verticalalignment = 'center', horizontalalignment = 'center', path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")])
    ax2.text(0.5*factor, -0.5*factor, 'VI', fontsize = 35, verticalalignment = 'center', horizontalalignment = 'center', path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")])


    #band structure
    xi = np.linspace(-1.25, 1.25, 100)
    yi1 = xi**2
    yi2 = -xi**2 - 0.9
    yi3 = -xi**2 - 1.2
    ax3.plot(xi, yi1, color = 'k', lw = 1.5)
    ax3.plot(xi, yi2, color = 'k', lw = 1.5)
    ax3.plot(xi, yi3, color = 'k', lw = 1.5)
    ax3.set_ylim(-2, 1)
    ax3.annotate("", xy=(-0.1, 0), xycoords='data', xytext=(-0.1, -0.9), textcoords='data',
                 arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", lw = 1.5, color = 'r'))
    ax3.text(-0.2, -0.5, 'A', ha = 'right', va = 'center', size = 25, color = 'r')
    ax3.annotate("", xy=(0.1, 0), xycoords='data', xytext=(0.1, -1.2), textcoords='data',
                 arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", lw = 1.5, color = 'b'))
    ax3.text(0.2, -0.5, 'B', ha = 'left', va = 'center', size = 25, color = 'b')
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    for tic in ax3.xaxis.get_major_ticks(): tic.tick1On = tic.tick2On = False
    for tic in ax3.yaxis.get_major_ticks(): tic.tick1On = tic.tick2On = False

    #2D frequency
    ax4.set_ylabel('$\mathsf{\hslash \omega_2}$', fontsize = 25)
    ax4.set_xlabel('$\mathsf{\hslash \omega_1 \, = \, \hslash \omega_m}$', fontsize = 25)
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), visible=False)
    fake_2D_path = os.path.join(box_path, 'MX2', 'papers', 'CMDS of MoS2 thin films', 'figures', 'fake 2D frequency', 'data.p')
    data = wt.data.from_pickle(fake_2D_path)
    data.zoom(2)
    data.chop('w0', 'w1')
    data.normalize()
    xi = data.axes[0].points
    yi = data.axes[1].points
    zi = data.channels[0].values.T
    levels = np.linspace(0, 1, 200)
    ax4.contourf(xi, yi, zi, levels, cmap = wright_cm)
    #ax4.set_xlim(6000, 10000)
    #ax4.set_ylim(6000, 10000)
    diag_min = max(min(xi), min(yi))
    diag_max = min(max(xi), max(yi))
    ax4.plot([diag_min, diag_max],[diag_min, diag_max], color = 'k', lw = 1, ls = '--')
    ax4.axhline(6500, color = 'k', lw = 1, ls = '--')
    ax4.axhline(7500, color = 'k', lw = 1, ls = '--')
    ax4.axvline(6500, color = 'k', lw = 1, ls = '--')
    ax4.axvline(7500, color = 'k', lw = 1, ls = '--')
    offset = 150
    ax4.text(6500 + offset, 6500 - offset, 'AA', fontsize=25, ha = 'left', va = 'top')
    ax4.text(7500 + offset, 6500 - offset, 'BA', fontsize=25, ha = 'left', va = 'top')
    ax4.text(6500 + offset, 7500 - offset, 'AB', fontsize=25, ha = 'left', va = 'top')
    ax4.text(7500 + offset, 7500 - offset, 'BB', fontsize=25, ha = 'left', va = 'top')
    
    
    '''
    #band structure
    ax4_sub = add_subplot_axes(ax4, [0.7, 0.7, 0.5, 0.5])
    xi = np.linspace(-1.25, 1.25, 100)
    yi1 = xi**2
    yi2 = -xi**2 - 0.9
    yi3 = -xi**2 - 1.2
    ax4_sub.plot(xi, yi1, color = 'k', lw = 1.5)
    ax4_sub.plot(xi, yi2, color = 'k', lw = 1.5)
    ax4_sub.plot(xi, yi3, color = 'k', lw = 1.5)
    ax4_sub.set_ylim(-2, 1)
    ax4_sub.annotate("", xy=(-0.1, 0), xycoords='data', xytext=(-0.1, -0.9), textcoords='data',
                 arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", lw = 1.5, color = 'r'))
    ax4_sub.text(-0.2, -0.5, 'A', ha = 'right', va = 'center', size = 20, color = 'r')
    ax4_sub.annotate("", xy=(0.1, 0), xycoords='data', xytext=(0.1, -1.2), textcoords='data',
                 arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", lw = 1.5, color = 'b'))
    ax4_sub.text(0.2, -0.5, 'B', ha = 'left', va = 'center', size = 20, color = 'b')
    #plt.setp(ax3.get_xticklabels(), visible=False)
    #plt.setp(ax3.get_yticklabels(), visible=False)
    #for tic in ax3.xaxis.get_major_ticks(): tic.tick1On = tic.tick2On = False
    #for tic in ax3.yaxis.get_major_ticks(): tic.tick1On = tic.tick2On = False
    '''

    #labels
    props = dict(boxstyle='square', facecolor='white', alpha=0.5)
    ax1.text(0.05, 0.95, 'a', transform=ax1.transAxes, fontsize=20, verticalalignment='top', horizontalalignment='left', bbox=props)
    #ax1_subB.text(0.05, 0.9, 'B', transform=ax1_subB.transAxes, fontsize=20, verticalalignment='top', horizontalalignment='left', bbox=props)
    ax2.text(0.05, 0.95, 'b', transform=ax2.transAxes, fontsize=20, verticalalignment='top', horizontalalignment='left', bbox=props)    
    ax3.text(0.05, 0.95, 'c', transform=ax3.transAxes, fontsize=20, verticalalignment='top', horizontalalignment='left', bbox=props)
    ax4.text(0.05, 0.95, 'd', transform=ax4.transAxes, fontsize=20, verticalalignment='top', horizontalalignment='left', bbox=props)

    #plt.tight_layout()

    plt.savefig('spectroscopy overview.tif', transparent = True, dpi = 300)
    plt.close('all')

### 227 uJ figures #############################################################

if False:
    data = data_at_fluence.chop('w2', 'w1', verbose = False)[0]
    artist = wt.artists.mpl_2D(data_at_fluence, xaxis = 'w2', yaxis = 'd2')
    artist.plot(normalize_slices = 'horizontal')
    plt.show()
    
if False:
    #227 uJ 9 2D frequencies

    plt.figure(figsize=(15.5, 15))
    gs0 = grd.GridSpec(1, 2, width_ratios=(30, 1))

    gs00 = grd.GridSpecFromSubplotSpec(3, 3, width_ratios=(1,1,1), height_ratios=(1,1,1), subplot_spec=gs0[0], hspace = 0.05, wspace = 0.05)
    ax1 = plt.subplot(gs00[0])
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2 = plt.subplot(gs00[1])
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax3 = plt.subplot(gs00[2])
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    ax4 = plt.subplot(gs00[3])
    plt.setp(ax4.get_xticklabels(), visible=False)
    ax5 = plt.subplot(gs00[4])
    plt.setp(ax5.get_xticklabels(), visible=False)
    plt.setp(ax5.get_yticklabels(), visible=False)
    ax6 = plt.subplot(gs00[5])
    plt.setp(ax6.get_xticklabels(), visible=False)
    plt.setp(ax6.get_yticklabels(), visible=False)
    ax7 = plt.subplot(gs00[6])
    ax8 = plt.subplot(gs00[7])
    plt.setp(ax8.get_yticklabels(), visible=False)
    ax9 = plt.subplot(gs00[8])
    plt.setp(ax9.get_yticklabels(), visible=False)
    subplots = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
    
    gs01 = grd.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1])
    ax_cb = plt.subplot(gs01[0])
    
    #get list of data objects
    delays = [80, 0, -80, -160, -240, -320, -400, -480, -560]
    data_objects = []
    for delay in delays:
        data_object = data_at_fluence.chop('w2', 'w1', {'d2': [delay, 'fs']}, verbose = False)[0]   
        data_objects.append(data_object)    
    
    for i in range(9):
        #get subplot
        subplot = subplots[i]
        #get data
        xi = data_objects[i].axes[1].points
        yi = data_objects[i].axes[0].points
        zi = data_objects[i].channels[0].values
        levels = np.linspace(data_min, data_max, 200)
        if True:
            X, Y, Z = wt.artists.pcolor_helper(xi, yi, zi)
            plt.subplot(subplot)
            cax = plt.pcolor(X, Y, Z, cmap = wright_cm, vmin = 0, vmax = 1)
            plt.xlim(xi.min(), xi.max())
            plt.ylim(yi.min(), yi.max())
        else:
            plt.subplot(subplot)
            cax = subplot.contourf(xi, yi, zi, levels, cmap = wright_cm)
        diag_min = max(min(xi), min(yi))
        diag_max = min(max(xi), max(yi))
        plt.plot([diag_min, diag_max],[diag_min, diag_max], color = 'k', lw = 1, ls = '--') #'k:')
        props = dict(boxstyle='square', facecolor='white', alpha=0.5)
        subplot.text(0.925, 0.075, str(delays[i]) + ' fs', transform=subplot.transAxes, fontsize=20, 
                     verticalalignment='bottom', horizontalalignment='right', bbox=props)
        #plt.grid()
        #vertical and horizontal lines
        line_alpha = 1
        plt.axhline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        plt.axhline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        plt.axvline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        plt.axvline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        #sideplots
        divider = make_axes_locatable(subplot)
        #x
        axCorrx = divider.append_axes('top', 0.75, pad=0.0, sharex=subplot)
        axCorrx.autoscale(False)
        axCorrx.set_adjustable('box-forced')
        plt.setp(axCorrx.get_xticklabels(), visible=False)
        axCorrx.get_yaxis().set_visible(False)
        axCorrx.axhline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorrx.axhline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorrx.axvline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorrx.axvline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        x_ax_int = zi.sum(axis = 0)
        x_ax_int -= x_ax_int.min()
        x_ax_int /= x_ax_int.max()
        axCorrx.set_ylim([-0.1,1.1])
        axCorrx.plot(xi, x_ax_int, color = 'b', lw = 1.5)
        #x_abs
        s_xi = abs_data.axes[0].points
        s_zi = abs_data.channels[0].values
        xlim =  axCorrx.get_xlim()
        min_index = np.argmin(abs(s_xi - min(xlim)))
        max_index = np.argmin(abs(s_xi - max(xlim)))
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi - min(s_zi_in_range)
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi / max(s_zi_in_range) 
        axCorrx.plot(s_xi, s_zi, lw = 1.5, color = 'k')
        #y
        axCorry = divider.append_axes('right', 0.75, pad=0.0, sharey=subplot)
        axCorry.autoscale(False)
        axCorry.set_adjustable('box-forced')
        plt.setp(axCorry.get_yticklabels(), visible=False)
        axCorry.get_xaxis().set_visible(False)
        axCorry.axhline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorry.axhline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorry.axvline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorry.axvline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        y_ax_int = zi.sum(axis = 1)
        y_ax_int = (y_ax_int - min(y_ax_int)) / (max(y_ax_int) - min(y_ax_int))
        axCorry.set_xlim([-0.1,1.1])
        axCorry.plot(y_ax_int, yi, color = 'b', lw = 1.5)
        #y_abs
        s_xi = abs_data.axes[0].points
        s_zi = abs_data.channels[0].values
        xlim =  axCorry.get_ylim()
        min_index = np.argmin(abs(s_xi - min(xlim)))
        max_index = np.argmin(abs(s_xi - max(xlim)))
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi - min(s_zi_in_range)
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi / max(s_zi_in_range) 
        axCorry.plot(s_zi, s_xi, lw = 1.5, color = 'k')
        
    #colorbar
    plt.subplot(ax_cb)
    plt.colorbar(cax, cax = ax_cb, ticks=np.linspace(0, 1, 11))
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.075, right=0.95, top=0.95, bottom=0.055)
    
    
    plt.figtext(0.01, 0.5, '$\mathsf{\hslash \omega_2 \,\, (eV)}$', rotation = 90, size = 30)
    plt.figtext(0.5, 0.01, '$\mathsf{\hslash \omega_1 \, = \, \hslash \omega_m \,\, (eV)}$', size = 30, horizontalalignment = 'center')
    
    plt.savefig('longtime 2D frequencies.tif', transparent=True, dpi=300)
    plt.close()

if False:
    #227 uJ 6 2D frequencies

    plt.figure(figsize=(15.5, 10))
    gs0 = grd.GridSpec(1, 2, width_ratios=(30, 1))

    gs00 = grd.GridSpecFromSubplotSpec(2, 3, width_ratios=(1,1,1), height_ratios=(1,1), subplot_spec=gs0[0], hspace = 0.05, wspace = 0.05)
    ax1 = plt.subplot(gs00[0])
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2 = plt.subplot(gs00[1])
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax3 = plt.subplot(gs00[2])
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    ax4 = plt.subplot(gs00[3])
    ax5 = plt.subplot(gs00[4])
    plt.setp(ax5.get_yticklabels(), visible=False)
    ax6 = plt.subplot(gs00[5])
    plt.setp(ax6.get_yticklabels(), visible=False)
    subplots = [ax1, ax2, ax3, ax4, ax5, ax6]
    
    gs01 = grd.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1])
    ax_cb = plt.subplot(gs01[0])
    
    #get list of data objects
    delays = [120, 80, 40, 0, -40, -80]
    data_objects = []
    for delay in delays:
        data_object = data_at_fluence.chop('w2', 'w1', {'d2': [delay, 'fs']}, verbose = False)[0]   
        data_object.normalize()
        data_objects.append(data_object)    
    
    for i in range(6):
        #get subplot
        subplot = subplots[i]
        #get data
        xi = data_objects[i].axes[1].points
        yi = data_objects[i].axes[0].points
        zi = data_objects[i].channels[0].values
        levels = np.linspace(0, 1, 200)
        if True:
            X, Y, Z = wt.artists.pcolor_helper(xi, yi, zi)
            plt.subplot(subplot)
            cax = plt.pcolor(X, Y, Z, cmap = wright_cm, vmin = 0, vmax = 1)
            plt.xlim(xi.min(), xi.max())
            plt.ylim(yi.min(), yi.max())
        else:
            plt.subplot(subplot)
            cax = subplot.contourf(xi, yi, zi, levels, cmap = wright_cm)
        diag_min = max(min(xi), min(yi))
        diag_max = min(max(xi), max(yi))
        plt.plot([diag_min, diag_max],[diag_min, diag_max], color = 'k', lw = 1, ls = '--') #'k:')
        props = dict(boxstyle='square', facecolor='white', alpha=0.5)
        subplot.text(0.925, 0.075, str(delays[i]) + ' fs', transform=subplot.transAxes, fontsize=20, 
                     verticalalignment='bottom', horizontalalignment='right', bbox=props)
        #plt.grid()
        #vertical and horizontal lines
        line_alpha = 1
        plt.axhline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        plt.axhline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        plt.axvline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        plt.axvline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        #sideplots
        divider = make_axes_locatable(subplot)
        #x
        axCorrx = divider.append_axes('top', 0.75, pad=0.0, sharex=subplot)
        axCorrx.autoscale(False)
        axCorrx.set_adjustable('box-forced')
        plt.setp(axCorrx.get_xticklabels(), visible=False)
        axCorrx.get_yaxis().set_visible(False)
        axCorrx.axhline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorrx.axhline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorrx.axvline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorrx.axvline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        x_ax_int = zi.sum(axis = 0)
        x_ax_int = (x_ax_int - min(x_ax_int)) / (max(x_ax_int) - min(x_ax_int))
        axCorrx.set_ylim([-0.1,1.1])
        axCorrx.plot(xi, x_ax_int, color = 'b', lw = 1.5)
        #x_abs
        s_xi = abs_data.axes[0].points
        s_zi = abs_data.channels[0].values
        xlim =  axCorrx.get_xlim()
        min_index = np.argmin(abs(s_xi - min(xlim)))
        max_index = np.argmin(abs(s_xi - max(xlim)))
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi - min(s_zi_in_range)
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi / max(s_zi_in_range) 
        axCorrx.plot(s_xi, s_zi, lw = 1.5, color = 'k')
        #y
        axCorry = divider.append_axes('right', 0.75, pad=0.0, sharey=subplot)
        axCorry.autoscale(False)
        axCorry.set_adjustable('box-forced')
        plt.setp(axCorry.get_yticklabels(), visible=False)
        axCorry.get_xaxis().set_visible(False)
        axCorry.axhline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorry.axhline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorry.axvline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorry.axvline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        y_ax_int = zi.sum(axis = 1)
        y_ax_int = (y_ax_int - min(y_ax_int)) / (max(y_ax_int) - min(y_ax_int))
        axCorry.set_xlim([-0.1,1.1])
        axCorry.plot(y_ax_int, yi, color = 'b', lw = 1.5)
        #y_abs
        s_xi = abs_data.axes[0].points
        s_zi = abs_data.channels[0].values
        xlim =  axCorry.get_ylim()
        min_index = np.argmin(abs(s_xi - min(xlim)))
        max_index = np.argmin(abs(s_xi - max(xlim)))
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi - min(s_zi_in_range)
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi / max(s_zi_in_range) 
        axCorry.plot(s_zi, s_xi, lw = 1.5, color = 'k')
        
    #colorbar
    plt.subplot(ax_cb)
    plt.colorbar(cax, cax = ax_cb, ticks=np.linspace(0, 1, 11))
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.075, right=0.95, top=0.95, bottom=0.085)
    
    
    plt.figtext(0.01, 0.5, '$\mathsf{\hslash \omega_2 \,\, (eV)}$', rotation = 90, size = 30)
    plt.figtext(0.5, 0.01, '$\mathsf{\hslash \omega_1 \, = \, \hslash \omega_m \,\, (eV)}$', size = 30, horizontalalignment = 'center')
    
    plt.savefig('earlytime 2D frequencies.tif', transparent=True, dpi=300)
    plt.close()

if False:
    #227 uJ w1 wigners

    plt.figure(figsize=(10, 15))

    gs0 = grd.GridSpec(1, 2, width_ratios=(30, 1), hspace = 0.05, wspace = 0.05)

    gs00 = grd.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs0[0], hspace = 0, wspace = 0.05)
    
    ax1 = plt.subplot(gs00[0])
    ax2 = plt.subplot(gs00[1])
    ax3 = plt.subplot(gs00[2])
    ax4 = plt.subplot(gs00[3])
    ax5 = plt.subplot(gs00[4])
    subplots = [ax1, ax2, ax3, ax4, ax5]
    
    gs01 = grd.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1])
    ax_cb = plt.subplot(gs01[0])
    
    #get list of data objects
    data_objects = []
    for energy in energies:
        data_object = data_at_fluence.chop('d2', 'w1', {'w2': [energy, 'wn']}, verbose = False)[0]
        data_object.normalize()
        data_objects.append(data_object)
        
    for i in range(5):
        subplot = subplots[i]
        #get data
        xi = data_objects[i].axes[1].points
        yi = data_objects[i].axes[0].points
        zi = data_objects[i].channels[0].values
        levels = np.linspace(0, 1, 200)
        if True:
            X, Y, Z = wt.artists.pcolor_helper(xi, yi, zi)
            plt.subplot(subplot)
            cax = plt.pcolor(X, Y, Z, cmap = wright_cm, vmin = levels.min(), vmax = levels.max())
            plt.xlim(xi.min(), xi.max())
            plt.ylim(yi.min(), yi.max())
        else:
            cax = subplot.contourf(xi, yi, zi, levels, cmap = wright_cm)
        yticks = [200, 100, 0, -100, -200, -300, -400, -500]
        plt.yticks(yticks)
        for val in yticks:
            plt.axhline(val, color = 'k', alpha = 1, lw = 1, ls = ':')
        #vertical lines
        plt.axvline(A_eV, color = 'k', alpha = 1, lw = 1, ls = '--')
        plt.axvline(B_eV, color = 'k', alpha = 1, lw = 1, ls = '--')
        w_eV = wt.units.converter(energies[i], 'wn', 'eV')
        plt.axvline(w_eV, color = 'k', alpha = 1, lw = 3)
            
    #colorbar
    plt.subplot(ax_cb)
    plt.colorbar(cax, cax = ax_cb)
    
    plt.subplots_adjust(bottom = 0.06, left = 0.15)    
    
    plt.figtext(0.01, 0.5, r'$\mathsf{\tau_{21} \,\, (fs)}$', rotation = 90, size = 30)
    plt.figtext(0.5, 0.01, '$\mathsf{\hslash \omega_1 \, = \, \hslash \omega_m \,\, (eV)}$', size = 30, horizontalalignment = 'center')
    
    plt.savefig('w1 wigners.tif', transparent=True, dpi=300)
    plt.close()
    
if False:
    #227 uJ w2 wigners

    plt.figure(figsize=(10, 15))

    gs0 = grd.GridSpec(1, 2, width_ratios=(30, 1), hspace = 0.05, wspace = 0.05)

    gs00 = grd.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs0[0], hspace = 0, wspace = 0.05)
    
    ax1 = plt.subplot(gs00[0])
    ax2 = plt.subplot(gs00[1])
    ax3 = plt.subplot(gs00[2])
    ax4 = plt.subplot(gs00[3])
    ax5 = plt.subplot(gs00[4])
    subplots = [ax1, ax2, ax3, ax4, ax5]
    
    gs01 = grd.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1])
    ax_cb = plt.subplot(gs01[0])
    
    #get list of data objects
    data_objects = []
    for energy in energies:
        data_object = data_at_fluence.chop('d2', 'w2', {'w1': [energy, 'wn']}, verbose = False)[0]
        data_object.normalize()        
        data_objects.append(data_object)
        
    for i in range(5):
        subplot = subplots[i]
        #get data
        xi = data_objects[i].axes[1].points
        yi = data_objects[i].axes[0].points
        zi = data_objects[i].channels[0].values
        levels = np.linspace(0, 1, 200)
        if True:
            X, Y, Z = wt.artists.pcolor_helper(xi, yi, zi)
            plt.subplot(subplot)
            cax = plt.pcolor(X, Y, Z, cmap = wright_cm, vmin = levels.min(), vmax = levels.max())
            plt.xlim(xi.min(), xi.max())
            plt.ylim(yi.min(), yi.max())
        else:
            cax = subplot.contourf(xi, yi, zi, levels, cmap = wright_cm)
        yticks = [200, 100, 0, -100, -200, -300, -400, -500]
        plt.yticks(yticks)
        for val in yticks:
            plt.axhline(val, color = 'k', alpha = 1, lw = 1, ls = ':')
        #vertical lines
        plt.axvline(A_eV, color = 'k', alpha = 1, lw = 1, ls = '--')
        plt.axvline(B_eV, color = 'k', alpha = 1, lw = 1, ls = '--')
        w_eV = wt.units.converter(energies[i], 'wn', 'eV')
        plt.axvline(w_eV, color = 'k', alpha = 1, lw = 3)
            
    #colorbar
    plt.subplot(ax_cb)
    plt.colorbar(cax, cax = ax_cb)
    
    plt.subplots_adjust(bottom = 0.06, left = 0.15)    
    
    plt.figtext(0.01, 0.5, r'$\mathsf{\tau_{21} \,\, (fs)}$', rotation = 90, size = 30)
    plt.figtext(0.5, 0.01, '$\mathsf{\hslash \omega_2 \,\, (eV)}$', size = 30, horizontalalignment = 'center')
    
    plt.savefig('w2 wigners.tif', transparent=True, dpi=300)
    plt.close()
    
   
### power study figures ########################################################


if False:
    #power study wigners

    plt.figure(figsize=(20, 15))

    gs0 = grd.GridSpec(1, 2, width_ratios=(30, 1), hspace = 0.05, wspace = 0.05)

    gs00 = grd.GridSpecFromSubplotSpec(5, 2, subplot_spec=gs0[0], hspace = 0, wspace = 0)
    
    axa1 = plt.subplot(gs00[0, 0])
    axa2 = plt.subplot(gs00[1, 0])
    axa3 = plt.subplot(gs00[2, 0])
    axa4 = plt.subplot(gs00[3, 0])
    axa5 = plt.subplot(gs00[4, 0])
    subplots_a = [axa1, axa2, axa3, axa4, axa5, 14480]
    
    axb1 = plt.subplot(gs00[0, 1])
    axb2 = plt.subplot(gs00[1, 1])
    axb3 = plt.subplot(gs00[2, 1])
    axb4 = plt.subplot(gs00[3, 1])
    axb5 = plt.subplot(gs00[4, 1])
    subplots_b = [axb1, axb2, axb3, axb4, axb5, 15690]
    for subplot in subplots_b[:5]:
        plt.setp(subplot.get_yticklabels(), visible=False)
    
    gs01 = grd.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1])
    ax_cb = plt.subplot(gs01[0])
    
    #get list of data objects

    for subplots in [subplots_a, subplots_b]:

        data_objects = data.chop('d2', 'w1', {'w2': [subplots[5], 'wn']}, verbose = True)
        for i in range(5):
            subplot = subplots[i]
            plt.subplot(subplot)
            #get data
            data_objects[i].level(0, 'd2', -3)
            xi = data_objects[i].axes[1].points
            yi = data_objects[i].axes[0].points
            zi = data_objects[i].channels[0].values
            zi -= zi.min()
            zi /= zi.max()
            levels = np.linspace(zi.min(), zi.max(), 200)
            if True:
                X, Y, Z = wt.artists.pcolor_helper(xi, yi, zi)
                plt.subplot(subplot)
                cax = plt.pcolor(X, Y, Z, cmap = wright_cm, vmin = levels.min(), vmax = levels.max())
                plt.xlim(xi.min(), xi.max())
                plt.ylim(yi.min(), yi.max())
            else:
                cax = subplot.contourf(xi, yi, zi, levels, cmap = wright_cm)
            yticks = [200, 100, 0, -100, -200, -300, -400, -500]
            plt.yticks(yticks)
            for val in yticks:
                plt.axhline(val, color = 'k', alpha = 1, lw = 1, ls = ':')
            #lines
            plt.axvline(A_eV, color = 'k', alpha = 1, lw = 1, ls = '--')
            plt.axvline(B_eV, color = 'k', alpha = 1, lw = 1, ls = '--')
            w_eV = wt.kit.unit_converter(subplots[5], 'wn', 'eV')
            plt.axvline(w_eV, color = 'k', alpha = 1, lw = 3)
            props = dict(boxstyle='square', facecolor='white', alpha=0.5)
            subplot.text(0.025, 0.925, str(fluencies[i]), transform=subplot.transAxes, fontsize=20, 
                         verticalalignment='top', horizontalalignment='left', bbox=props)        
            
    #colorbar$\mathsf{\hslash \omega_2 \,\, (eV)}$
    plt.subplot(ax_cb)
    plt.colorbar(cax, cax = ax_cb)
    
    plt.subplots_adjust(bottom = 0.05, left = 0.075)    
    
    plt.figtext(0.01, 0.5, r'$\mathsf{\tau_{21} \,\, (fs)}$', rotation = 90, size = 30)
    plt.figtext(0.5, 0.01, '$\mathsf{\hslash \omega_1 \, = \, \hslash \omega_m \,\, (eV)}$', size = 30, horizontalalignment = 'center')
    
    plt.savefig('power study wigners.png', transparent = True, dpi = 300)
    plt.close('all')
    


### supporting information figures #############################################

def cubehelix_r(gamma=1.0, s=0.5, r=-1.5, h=1.0):
    def get_color_function(p0, p1):
        def color(x):
            # Apply gamma factor to emphasise low or high intensity values
            xg = x ** gamma

            # Calculate amplitude and angle of deviation from the black
            # to white diagonal in the plane of constant
            # perceived intensity.
            a = h * xg * (1 - xg) / 2

            phi = 2 * np.pi * (s / 3 + r * x)
            out = xg + a * (p0 * np.cos(phi) + p1 * np.sin(phi))

            return out[::-1]
        return color
    return {
            'red': get_color_function(-0.14861, 1.78277),
            'green': get_color_function(-0.29227, -0.90649),
            'blue': get_color_function(1.97294, 0.0),
    }

from matplotlib.colors import LinearSegmentedColormap
cm1_r = cubehelix_r(s=0.5, r=0.9, h=1.4, gamma=1)
a_r = LinearSegmentedColormap('chw1_r', cm1_r)

if False:
    #227 uJ movie
    plt.close('all')
    datas = data_at_fluence.chop('w2', 'w1')
    datas = datas[::-1]
    
    for i, data in zip(range(len(datas)), datas):
        d2 = np.round(data.constants[2].points)
        
        fig = plt.figure(figsize=(8, 7))
    
        gs = grd.GridSpec(1, 2, width_ratios=[20, 1], wspace=0.1)            
            
        subplot_main = plt.subplot(gs[0])
        subplot_main.patch.set_facecolor('grey')
        
        #main figure
        line_alpha = 1
        xi = data.w1.points
        yi = data.w2.points
        zi = data.channels[0].values
        X, Y, Z = wt.artists.pcolor_helper(xi, yi, zi)
        
        cax = plt.pcolor(X, Y, Z, vmin = 0, vmax = 1, cmap = 'gnuplot2_r')
        plt.xlim(min(xi), max(xi))
        plt.ylim(min(yi), max(yi))
        diag_min = max(min(xi), min(yi))
        diag_max = min(max(xi), max(yi))
        plt.plot([diag_min, diag_max],[diag_min, diag_max], color = 'k', lw = 1, ls = '--') #'k:')
        plt.axhline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        plt.axhline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        plt.axvline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        plt.axvline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        plt.ylabel('$\mathsf{\hslash \omega_2 \,\, (eV)}$', size = 18)
        plt.xlabel('$\mathsf{\hslash \omega_1 \, = \, \hslash \omega_m \,\, (eV)}$', size = 18)

        #sideplots
        subplot = subplot_main
        divider = make_axes_locatable(subplot)
        #x
        axCorrx = divider.append_axes('top', 0.75, pad=0.0, sharex=subplot)
        axCorrx.autoscale(False)
        axCorrx.set_adjustable('box-forced')
        plt.setp(axCorrx.get_xticklabels(), visible=False)
        axCorrx.get_yaxis().set_visible(False)
        axCorrx.axhline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorrx.axhline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorrx.axvline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorrx.axvline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        x_ax_int = zi.sum(axis = 0)
        x_ax_int = (x_ax_int - min(x_ax_int)) / (max(x_ax_int) - min(x_ax_int))
        axCorrx.set_ylim([-0.1,1.1])
        axCorrx.plot(xi, x_ax_int, color = 'b', lw = 1.5)
        #x_abs
        s_xi = abs_data.axes[0].points
        s_zi = abs_data.channels[0].values
        xlim =  axCorrx.get_xlim()
        min_index = np.argmin(abs(s_xi - min(xlim)))
        max_index = np.argmin(abs(s_xi - max(xlim)))
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi - min(s_zi_in_range)
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi / max(s_zi_in_range) 
        axCorrx.plot(s_xi, s_zi, lw = 1.5, color = 'k')
        #y
        axCorry = divider.append_axes('right', 0.75, pad=0.0, sharey=subplot)
        axCorry.autoscale(False)
        axCorry.set_adjustable('box-forced')
        plt.setp(axCorry.get_yticklabels(), visible=False)
        axCorry.get_xaxis().set_visible(False)
        axCorry.axhline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorry.axhline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorry.axvline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorry.axvline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        y_ax_int = zi.sum(axis = 1)
        y_ax_int = (y_ax_int - min(y_ax_int)) / (max(y_ax_int) - min(y_ax_int))
        axCorry.set_xlim([-0.1,1.1])
        axCorry.plot(y_ax_int, yi, color = 'b', lw = 1.5)
        #y_abs
        s_xi = abs_data.axes[0].points
        s_zi = abs_data.channels[0].values
        xlim =  axCorry.get_ylim()
        min_index = np.argmin(abs(s_xi - min(xlim)))
        max_index = np.argmin(abs(s_xi - max(xlim)))
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi - min(s_zi_in_range)
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi / max(s_zi_in_range) 
        axCorry.plot(s_zi, s_xi, lw = 1.5, color = 'k')
        
        plt.suptitle(r'$\mathsf{\tau_{21} \, = \, %d \, fs}$'%d2)
        
        
        #colorbar
        subplot_cb = plt.subplot(gs[1])
        plt.colorbar(cax, cax = subplot_cb)
        
        
        
        plt.savefig(r'movie\{}.png'.format(str(i).zfill(3)), transparent = False, dpi = 300)
        plt.close('all')
        print i

if False:
    #data workup example
    data_raw = wt.data.from_pickle(data_pickle_path)

    plt.figure(figsize=(15.5, 10))
    gs0 = grd.GridSpec(1, 2, width_ratios=(30, 1))

    gs00 = grd.GridSpecFromSubplotSpec(2, 3, width_ratios=(1,1,1), height_ratios=(1,1), subplot_spec=gs0[0], hspace = 0.05, wspace = 0.05)
    ax1 = plt.subplot(gs00[0])
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2 = plt.subplot(gs00[1])
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax3 = plt.subplot(gs00[2])
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    ax4 = plt.subplot(gs00[3])
    ax5 = plt.subplot(gs00[4])
    plt.setp(ax5.get_yticklabels(), visible=False)
    ax6 = plt.subplot(gs00[5])
    plt.setp(ax6.get_yticklabels(), visible=False)
    subplots = [ax1, ax2, ax3, ax4, ax5, ax6]
    
    gs01 = grd.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1])
    ax_cb = plt.subplot(gs01[0])
    
    #get list of data objects  
    data_at_delay = data_raw.chop('w2', 'w1', {'d2': [0, 'fs'], 'fluence': [227, 'uJ per sq. cm']})[0]
    data_at_delay.convert('eV')
    
    zis = []
    #detector voltages
    zi = data_at_delay.channels[4].values.copy()
    zi -= zi.min()
    zi /= zi.max()
    zis.append(zi)
    #chopped portion
    zi = data_at_delay.channels[0].values.copy()
    zi -= zi.min()
    zi /= zi.max()
    zis.append(zi)
    #invarient portion
    points = data_raw.d2.points[-3:]    
    zis_to_average = []
    for point in points:
        zis_to_average.append(data_raw.chop('w2', 'w1', {'d2': [point, 'fs'], 'fluence': [227, 'uJ per sq. cm']})[0].channels[0].values)
    zis_to_average = np.array(zis_to_average)
    print zis_to_average.shape
    a = np.average(zis_to_average, axis = 0)
    a -= a.min()
    a /= a.max()
    zis.append(a)
    #subtracted portion
    data_raw.level(0, 'd2', -3)
    data_at_delay = data_raw.chop('w2', 'w1', {'d2': [0, 'fs'], 'fluence': [227, 'uJ per sq. cm']})[0]  
    data_at_delay.normalize()
    zis.append(data_at_delay.channels[0].values.copy())
    #smoothed
    data_at_delay.smooth(2)
    data_at_delay.normalize()
    zis.append(data_at_delay.channels[0].values.copy())
    #amplitude
    data_at_delay.scale()
    data_at_delay.normalize()
    zis.append(data_at_delay.channels[0].values.copy())
    
    for i in range(6):
        #get subplot
        subplot = subplots[i]
        #get data
        xi = data_at_delay.w1.points
        yi = data_at_delay.w2.points
        zi = zis[i]
        levels = np.linspace(0, 1, 200)
        if True:
            X, Y, Z = wt.artists.pcolor_helper(xi, yi, zi)
            plt.subplot(subplot)
            cax = plt.pcolor(X, Y, Z, cmap = wright_cm, vmin = levels.min(), vmax = levels.max())
            plt.xlim(xi.min(), xi.max())
            plt.ylim(yi.min(), yi.max())
        else:
            plt.subplot(subplot)
            cax = subplot.contourf(xi, yi, zi, levels, cmap = wright_cm)
        diag_min = max(min(xi), min(yi))
        diag_max = min(max(xi), max(yi))
        plt.plot([diag_min, diag_max],[diag_min, diag_max], color = 'k', lw = 1, ls = '--') #'k:')
        props = dict(boxstyle='square', facecolor='white', alpha=0.5)
        labels = ['a', 'b', 'c', 'd', 'e', 'f']
        subplot.text(0.1, 0.875, str(labels[i]), transform=subplot.transAxes, fontsize=20, 
                     verticalalignment='bottom', horizontalalignment='right', bbox=props)
        #plt.grid()
        #vertical and horizontal lines
        line_alpha = 1
        plt.axhline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        plt.axhline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        plt.axvline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        plt.axvline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        #sideplots
        divider = make_axes_locatable(subplot)
        #x
        axCorrx = divider.append_axes('top', 0.75, pad=0.0, sharex=subplot)
        axCorrx.autoscale(False)
        axCorrx.set_adjustable('box-forced')
        plt.setp(axCorrx.get_xticklabels(), visible=False)
        axCorrx.get_yaxis().set_visible(False)
        axCorrx.axhline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorrx.axhline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorrx.axvline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorrx.axvline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        x_ax_int = zi.sum(axis = 0)
        x_ax_int = (x_ax_int - min(x_ax_int)) / (max(x_ax_int) - min(x_ax_int))
        axCorrx.set_ylim([-0.1,1.1])
        axCorrx.plot(xi, x_ax_int, color = 'b', lw = 1.5)
        #x_abs
        s_xi = abs_data.axes[0].points
        s_zi = abs_data.channels[0].values
        xlim =  axCorrx.get_xlim()
        min_index = np.argmin(abs(s_xi - min(xlim)))
        max_index = np.argmin(abs(s_xi - max(xlim)))
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi - min(s_zi_in_range)
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi / max(s_zi_in_range) 
        axCorrx.plot(s_xi, s_zi, lw = 1.5, color = 'k')
        #y
        axCorry = divider.append_axes('right', 0.75, pad=0.0, sharey=subplot)
        axCorry.autoscale(False)
        axCorry.set_adjustable('box-forced')
        plt.setp(axCorry.get_yticklabels(), visible=False)
        axCorry.get_xaxis().set_visible(False)
        axCorry.axhline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorry.axhline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorry.axvline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorry.axvline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        y_ax_int = zi.sum(axis = 1)
        y_ax_int = (y_ax_int - min(y_ax_int)) / (max(y_ax_int) - min(y_ax_int))
        axCorry.set_xlim([-0.1,1.1])
        axCorry.plot(y_ax_int, yi, color = 'b', lw = 1.5)
        #y_abs
        s_xi = abs_data.axes[0].points
        s_zi = abs_data.channels[0].values
        xlim =  axCorry.get_ylim()
        min_index = np.argmin(abs(s_xi - min(xlim)))
        max_index = np.argmin(abs(s_xi - max(xlim)))
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi - min(s_zi_in_range)
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi / max(s_zi_in_range) 
        axCorry.plot(s_zi, s_xi, lw = 1.5, color = 'k')
    
    #colorbar
    plt.subplot(ax_cb)
    plt.colorbar(cax, cax = ax_cb, ticks=np.linspace(0, 1, 11))
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.075, right=0.95, top=0.95, bottom=0.075)
    
    plt.figtext(0.01, 0.5, '$\mathsf{\hslash \omega_2 \,\, (eV)}$', rotation = 90, size = 30)
    plt.figtext(0.4, 0.01, '$\mathsf{\hslash \omega_1 \, = \, \hslash \omega_m \,\, (eV)}$', size = 30)
    
    plt.savefig('data workup.tif', transparent=True, dpi=300)
    plt.close()

if False:
    #OPAs  
    
    plt.figure(figsize = (15, 10))
    
    gs1 = grd.GridSpec(2, 2, width_ratios = [30, 1])
    gs1.update(wspace=0.05, hspace=0)
    ax1 = plt.subplot(gs1[0, 0])
    ax2 = plt.subplot(gs1[1, 0], sharex = ax1)
    ax3 = plt.subplot(gs1[0:2, 1])
    
    #OPA1
    OPA1_dat = os.path.join(tuning_path, 'OPA1', 'OPA1 3505437107 Test Mixer2 NON-SH-NON-Sig [41x15] Topas.dat')
    dat = f.Dat(filepath = OPA1_dat, xvar = 'w1', yvar = 'wm', cols = 'v2', colortune = True)
    OPA1_xi = wt.units.converter(dat.xi, 'wn', 'eV')
    OPA1_yi = dat.yi
    OPA1_zi = dat.zi
    OPA1_zi -= OPA1_zi.min()
    OPA1_zi /= OPA1_zi.max()
    levels = np.linspace(0, 1, 200)
    ax1.contourf(OPA1_xi, OPA1_yi, OPA1_zi, levels, cmap = wright_cm)
    ax1.grid()
    ax1.set_xlim(1.61, 2.13)
    ax1.set_ylim(-1250, 1250)
    ax1.axhline(0, color ='k', lw = 1)
    plt.setp(ax1.get_xticklabels(), visible=False)
    props = dict(boxstyle='square', facecolor='white', alpha=0.5)
    ax1.text(0.02, 0.95, 'OPA1', transform=ax1.transAxes, fontsize=20, 
                 verticalalignment='top', horizontalalignment='left', bbox=props)    
    '''
    divider = make_axes_locatable(ax1)
    #x
    axCorrx = divider.append_axes('top', 0.75, pad=0.0, sharex=ax1)
    axCorrx.grid()
    axCorrx.autoscale(False)
    axCorrx.set_adjustable('box-forced')
    plt.setp(axCorrx.get_xticklabels(), visible=False)
    axCorrx.get_yaxis().set_visible(False)
    x_ax_int = OPA1_zi.sum(axis = 0)
    x_ax_int = (x_ax_int - min(x_ax_int)) / (max(x_ax_int) - min(x_ax_int))
    axCorrx.set_ylim([-0.1,1.1])
    axCorrx.plot(OPA1_xi, x_ax_int, color = 'b', lw = 1.5)
    '''
    
    #OPA2
    OPA2_dat = os.path.join(tuning_path, 'OPA2', 'OPA2 Mixer2 NON-SH-NON-Sig OPA2_SHS_tunetest_2015.01.27 VAI0 0 ms wait [41x15] Pts Test.dat')
    dat = f.Dat(filepath = OPA2_dat, xvar = 'w2', yvar = 'wm', cols = 'v2', colortune = True)
    OPA2_xi = wt.units.converter(dat.xi, 'wn', 'eV')
    OPA2_yi = dat.yi
    OPA2_zi = dat.zi
    OPA2_zi -= OPA2_zi.min()
    OPA2_zi /= OPA2_zi.max()
    cax = ax2.contourf(OPA2_xi, OPA2_yi, OPA2_zi, levels, cmap = wright_cm)
    ax2.grid()
    #ax2.set_xlim(1.61, 2.13)
    ax2.set_ylim(-1250, 1250)
    ax2.axhline(0, color ='k', lw = 1)
    ax2.text(0.02, 0.95, 'OPA2', transform=ax2.transAxes, fontsize=20, 
                 verticalalignment='top', horizontalalignment='left', bbox=props) 
    '''
    divider = make_axes_locatable(ax2)
    #x
    axCorrx = divider.append_axes('top', 0.75, pad=0.0, sharex=ax2)
    axCorrx.grid()
    axCorrx.autoscale(False)
    axCorrx.set_adjustable('box-forced')
    plt.setp(axCorrx.get_xticklabels(), visible=False)
    axCorrx.get_yaxis().set_visible(False)
    x_ax_int = OPA2_zi.sum(axis = 0)
    x_ax_int = (x_ax_int - min(x_ax_int)) / (max(x_ax_int) - min(x_ax_int))
    axCorrx.set_ylim([-0.1,1.1])
    axCorrx.plot(OPA2_xi, x_ax_int, color = 'b', lw = 1.5)
    '''
    
    plt.figtext(0.01, 0.65, '$\mathsf{detuning \,\, (cm^{-1})}$', rotation = 90, size = 30)
    plt.figtext(0.35, 0.01, '$\mathsf{OPA \, setpoint \,\, (eV)}$', size = 30)
    
    #colorbar
    plt.subplot(ax3)
    plt.colorbar(cax, cax = ax3, ticks=np.linspace(0, 1, 11))
    
    #plt.show()
    plt.savefig('OPA tune tests.tif', transparent=True, dpi=300)
    plt.close('all')


if False:
    # epi vs transmissive

    fig = plt.figure(figsize = (16, 6))
    
    gs0 = grd.GridSpec(1, 2, width_ratios=(1, 2), wspace = 0.2)

    gs00 = grd.GridSpecFromSubplotSpec(1, 2, width_ratios = [30, 1], subplot_spec=gs0[0], wspace = 0.05)
    ax0 = plt.subplot(gs00[0])
    plt.setp(ax0.get_xticklabels(), visible=False)
    plt.setp(ax0.get_yticklabels(), visible=False)

    gs01 = grd.GridSpecFromSubplotSpec(1, 2, width_ratios = [30, 1], subplot_spec=gs0[1], wspace = 0.05)
    ax1 = plt.subplot(gs01[0])
    plt.setp(ax1.get_yticklabels(), visible=False)
    ax2 = plt.subplot(gs01[1])
    
    #mask
    scalar = 9.
    mp = [[-0.0874*scalar, 0*scalar], #w1
          [0.0175*scalar, 0*scalar], #w2'
          [0*scalar, 0.0262*scalar]] #w2    
    label_offset = 0
    mec = matplotlib.colors.ColorConverter().to_rgba('k', alpha=1)
    #    
    si_alpha = 0.25
    ax0.annotate('$\mathsf{1}$', [mp[0][0], mp[0][1]], color='k', size = 20, va = 'center', ha = 'center')
    ax0.plot(mp[0][0], mp[0][1], color='r', marker='o', lw = 0, markeredgecolor = mec, markersize = 30, alpha = si_alpha)
    #
    ax0.annotate('$\mathsf{2}$', [mp[2][0], mp[2][1]], color='k', size = 20, va = 'center', ha = 'center')
    ax0.plot(mp[2][0], mp[2][1], color='b', marker='o', lw = 0, markeredgecolor = 'k', markersize = 30, alpha = si_alpha)
    #    
    ax0.annotate('$\mathsf{2^{\prime}}$', [mp[1][0], mp[1][1]], color='k', size = 20, va = 'center', ha = 'center')
    ax0.plot(mp[1][0], mp[1][1], color='b', marker='o', lw = 0, markeredgecolor = None, markersize = 30, alpha = si_alpha)
    #
    out_x = mp[0][0] + mp[1][0] - mp[2][0]
    out_y = mp[0][1] + mp[1][1] - mp[2][1]
    ax0.plot(out_x, out_y, color='r', marker='*', lw = 0, markeredgecolor = None, markersize = 30, alpha = si_alpha)
    #
    ax0.plot(0, 0, color='k', marker='o', lw = 0, markeredgecolor = None, markersize = 5)
    #    
    aspect = 4.35/4.
    ax0.set_xlim(-1.25, 0.75) 
    ax0.set_ylim(-1*aspect, 1*aspect)
    ax0.grid()
    
    #transmissive
    tra_dat_filepath = os.path.join(box_path, 'MX2', '2014.09.29 group work', '2014.10.09 zero delay workup', 'w1=w2=15820cm-1_absolute_zerodelay_MoS2onQuartz_transmissiveTRIEE  (2) 200 ms wait [21x21] Delay.dat')
    tra_data = wt.data.from_COLORS(tra_dat_filepath)
    tra_data.zoom(2)
    tra_data.smooth(3)
    tra_data.level(0, 'd2', -4)
    tra_data.scale()
    tra_data.normalize()
    xi = tra_data.d1.points - 25.
    yi = tra_data.d2.points - 10.
    zi = tra_data.channels[0].values.T
    zi -= 0.05
    plt.subplot(ax1)
    levels = np.linspace(0, 1, 200)
    cax = ax1.contourf(xi, yi, zi, levels, cmap = wright_cm)
    ax1.set_xlim(-175, 175)
    ax1.set_ylim(-175, 175)
    ax1.grid()
    diag_min = max(min(xi), min(yi))
    diag_max = min(max(xi), max(yi))
    ax1.plot([diag_min, diag_max],[diag_min, diag_max], 'k:')
    plt.xticks(rotation = 45)
    
    #epi
    epi_dat_filepath = os.path.join(box_path, 'MX2', '2014.09.29 group work', '2014.10.04 TrEE workup', 'absolute zero_quartz_mos2 (1) 0 ms wait [56x21] Delay.dat')
    epi_data = wt.data.from_COLORS(epi_dat_filepath)
    epi_data.zoom(2)
    epi_data.smooth(3)
    epi_data.level(0, 'd2', -4)
    epi_data.scale()
    epi_data.normalize()
    xi = epi_data.d1.points + 10.
    yi = epi_data.d2.points + 20.
    zi = epi_data.channels[0].values.T
    divider = make_axes_locatable(ax1)
    bbox = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    axCorry = divider.append_axes('right', 4., pad=0.0, sharey=ax1)
    #plt.subplot(axCorry)
    plt.setp(axCorry.get_yticklabels(), visible=False)
    axCorry.contourf(xi, yi, zi, levels, cmap = wright_cm)
    axCorry.set_xlim(-175, 175)
    axCorry.set_ylim(-175, 175)
    axCorry.grid()
    diag_min = max(min(xi), min(yi))
    diag_max = min(max(xi), max(yi))
    axCorry.plot([diag_min, diag_max],[diag_min, diag_max], 'k:')
    plt.xticks(rotation = 45)
    
    
    #colorbar
    plt.subplot(ax2)
    plt.colorbar(cax, cax = ax2, ticks=np.linspace(0, 1, 11))
    
    #plt.tight_layout()
    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.225)
    
    plt.figtext(0.36, 0.65, r'$\mathsf{\tau_{21} \,\, (fs)}$', rotation = 90, size = 30)
    plt.figtext(0.63, 0.025, r'$\mathsf{\tau_{22^{\prime}} \,\, (fs)}$', size = 30)
    
    props = dict(boxstyle='square', facecolor='white', alpha=0.5)
    ax1.text(0.05, 0.95, 'transmissive', transform=ax1.transAxes, fontsize=20, 
                 verticalalignment='top', horizontalalignment='left', bbox=props)
    axCorry.text(0.05, 0.95, 'reflective', transform=axCorry.transAxes, fontsize=20, 
                 verticalalignment='top', horizontalalignment='left', bbox=props)
    
    #labels
    props = dict(boxstyle='square', facecolor='white', alpha=0.5)
    ax0.text(-0.05, 0.975, 'a', transform=ax0.transAxes, fontsize=20, verticalalignment='top', horizontalalignment='right', bbox=props)
    ax1.text(-0.2, 0.975, 'b', transform=ax1.transAxes, fontsize=20, verticalalignment='top', horizontalalignment='right', bbox=props)

                 
    plt.savefig('epi vs transmissive.tif', transparent=False, dpi=300)
    plt.close('all')
    
    for axis in [ax0, ax1, axCorry]:
        bbox = axis.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        print bbox.width, bbox.height

if False:
    # dynamics

    fig = plt.figure(figsize = (33, 10))
    gs = grd.GridSpec(1, 2, width_ratios=(1, 2), wspace = 0.2)
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    # 2D frequency
    data_2D_freq = data_at_fluence.chop('w2', 'w1', {'d2': [-600, 'fs']})[0]
    data_2D_freq.normalize()
    xi = data_2D_freq.axes[1].points
    yi = data_2D_freq.axes[0].points
    zi = data_2D_freq.channels[0].values
    X, Y, Z = wt.artists.pcolor_helper(xi, yi, zi)
    plt.subplot(ax0)
    cax = plt.pcolor(X, Y, Z, cmap = wright_cm, vmin = 0, vmax = 1)
    plt.xlim(xi.min(), xi.max())
    plt.ylim(yi.min(), yi.max())
    #ax0.set_aspect('equal')
    diag_min = max(min(xi), min(yi))
    diag_max = min(max(xi), max(yi))
    plt.plot([diag_min, diag_max],[diag_min, diag_max], color = 'k', lw = 1, ls = '--')
    line_alpha = 1
    plt.axhline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
    plt.axhline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
    plt.axvline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
    plt.axvline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')

    # delay traces
    colors_list = [[A_eV, A_eV],
                   [A_eV, B_eV],
                   [B_eV, B_eV],
                   [B_eV, A_eV],
                   [A_eV, 2.04],
                   [B_eV, 2.04]]
    from matplotlib.pyplot import cm
    cs = cm.rainbow(np.linspace(0,1,len(colors_list)))
    for c, colors in zip(cs, colors_list):
        data_delay = data_at_fluence.chop('d2', {'w1': [colors[0], 'eV'], 
                                                 'w2': [colors[1], 'eV']})[0]
        data_delay.normalize()
        xi = data_delay.d2.points
        zi = data_delay.channels[0].values
        plt.subplot(ax1)
        plt.plot(xi, zi, lw=4, alpha=0.75, c=c)
        # also in ax0
        ax0.scatter(colors[0], colors[1], s=1000, edgecolor='w', linewidth=3, c=c)
    plt.grid()
    plt.xlim(250, -600)
    
    plt.savefig('dynamics.png', transparent = True, dpi = 300)
    plt.close('all')


if False:
    # zerotune
    zerotune_folderpath = os.path.join(box_path, 'MX2', '2015.02.08', 'Post-optimization data', '2015.02.04_optimized ZDC', 'Zerotune')
    w1d2_filepath = os.path.join(zerotune_folderpath, 'OPA1_D2_ZDC__ Mixer2 [29x24] Zero.fit')
    w2d1_filepath = os.path.join(zerotune_folderpath, 'OPA2_D1_ZDC__Mixer2 [29x24] Zero.fit')
    w2d2_filepath = os.path.join(zerotune_folderpath, 'OPA2_D2_ZDC__Mixer2 [29x24] Zero.fit')
    filepaths = [w1d2_filepath, w2d1_filepath, w2d2_filepath]   
    
    plt.figure(figsize=(10, 10))
    gs = grd.GridSpec(3, 1, hspace=0)
    
    ax0 = plt.subplot(gs[0])
    plt.setp(ax0.get_xticklabels(), visible=False)
    ax1 = plt.subplot(gs[1])
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2 = plt.subplot(gs[2])
    axes = [ax0, ax1, ax2]    
    
    yticks = [None, None, None]
    ylims = [None, None, None]
    yticks[0] = [-150, -100, -50, 0, 50, 100]
    ylims[0] = [-200, 150]
    yticks[1] = [-30, -20, -10, 0, 10]
    ylims[1] = [-40, 20]
    yticks[2] = [-100, -50, 0, 50, 100, 150, 200, 250]
    ylims[2] = [-150, 300]
    
    ylabels = [None, None, None]
    ylabels[0] = r'$\mathsf{\tau_{21} \/ (fs)}$'
    ylabels[1] = r'$\mathsf{\tau_{22^{\prime}} \/ (fs)}$'
    ylabels[2] = r'$\mathsf{\tau_{21} \/ (fs)}$'
    
    texts = ['OPA1', 'OPA2', 'OPA2']

    for i in range(3):
        
        fit = np.genfromtxt(filepaths[i]).T
        xi = wt.units.converter(fit[1], 'nm', 'eV')
        yi = fit[9]
        
        plt.subplot(axes[i])
        plt.plot(xi, yi, c='k', lw=2)
        plt.yticks(yticks[i])
        plt.ylim(ylims[i])
        plt.ylabel(ylabels[i], fontsize=25)
        plt.grid()
        plt.xlim(1.62, 2.12)
        
        subplot = axes[i]
        props = dict(boxstyle='square', facecolor='white', alpha=0.5)
        subplot.text(0.5, 0.9, texts[i], transform=subplot.transAxes, fontsize=20, 
                     verticalalignment='top', horizontalalignment='center', bbox=props)
  
    
    plt.xlabel('$\mathsf{\hslash \omega \/ (eV)}$', fontsize=25)
        
        
    plt.savefig('zerotune_new_labels.tif', transparent=True, dpi=300)
    plt.close('all')
    
    
### etc #######################################################################



if False:
    
    artist = wt.artists.mpl_2D(data_at_fluence, xaxis = 'w1', yaxis = 'w2')
    artist.sideplot(abs_data, y = False)
    artist.plot(local = True, pixelated = True, xbin = True, ybin = True,
                output_folder = 'images', contours=0)



if False:
    #data workup example for grant
    data_raw = wt.data.from_pickle(data_pickle_path)

    plt.figure(figsize=(11, 10))
    gs0 = grd.GridSpec(1, 2, width_ratios=(30, 1))

    gs00 = grd.GridSpecFromSubplotSpec(2, 2, width_ratios=(1,1), height_ratios=(1,1), subplot_spec=gs0[0], hspace = 0.05, wspace = 0.05)
    ax1 = plt.subplot(gs00[0])
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2 = plt.subplot(gs00[1])
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax3 = plt.subplot(gs00[2])
    ax4 = plt.subplot(gs00[3])
    plt.setp(ax4.get_yticklabels(), visible=False)    
    subplots = [ax1, ax2, ax3, ax4]
    
    gs01 = grd.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1])
    ax_cb = plt.subplot(gs01[0])
    
    #get list of data objects  
    data_at_delay = data_raw.chop('w2', 'w1', {'d2': [0, 'fs'], 'fluence': [227, 'uJ per sq. cm']})[0]
    data_at_delay.convert('eV')
    
    zis = []
    #detector voltages
    zis.append(data_at_delay.channels[4].values.copy())
    #chopped portion
    zis.append(data_at_delay.channels[0].values.copy())
    #invarient portion
    points = data_raw.d2.points[-3:]    
    zis_to_average = []
    for point in points:
        zis_to_average.append(data_raw.chop('w2', 'w1', {'d2': [point, 'fs'], 'fluence': [227, 'uJ per sq. cm']})[0].channels[0].values)
    zis_to_average = np.array(zis_to_average)
    print zis_to_average.shape
    a = np.average(zis_to_average, axis = 0)
    #zis.append(a)
    #subtracted portion
    data_raw.level(0, 'd2', -3)
    data_at_delay = data_raw.chop('w2', 'w1', {'d2': [0, 'fs'], 'fluence': [227, 'uJ per sq. cm']})[0]  
    zis.append(data_at_delay.channels[0].values.copy())
    #smoothed
    data_at_delay.smooth(2)
    zis.append(data_at_delay.channels[0].values.copy())
    #amplitude
    data_at_delay.scale()
    zis.append(data_at_delay.channels[0].values.copy())
    
    for i in range(4):
        #get subplot
        subplot = subplots[i]
        #get data
        xi = data_at_delay.w1.points
        yi = data_at_delay.w2.points
        zi = zis[i]
        zi -= zi.min()
        zi /= zi.max()
        levels = np.linspace(data_min, data_max, 200)
        if True:
            X, Y, Z = wt.artists.pcolor_helper(xi, yi, zi)
            plt.subplot(subplot)
            cax = plt.pcolor(X, Y, Z, cmap = wright_cm, vmin = levels.min(), vmax = levels.max())
            plt.xlim(xi.min(), xi.max())
            plt.ylim(yi.min(), yi.max())
        else:
            plt.subplot(subplot)
            cax = subplot.contourf(xi, yi, zi, levels, cmap = wright_cm)
        diag_min = max(min(xi), min(yi))
        diag_max = min(max(xi), max(yi))
        plt.plot([diag_min, diag_max],[diag_min, diag_max], 'k:')
        plt.grid()
        props = dict(boxstyle='square', facecolor='white', alpha=0.5)
        labels = ['A', 'B', 'C', 'D', 'E', 'F']
        subplot.text(0.1, 0.875, str(labels[i]), transform=subplot.transAxes, fontsize=8, 
                     verticalalignment='bottom', horizontalalignment='right', bbox=props)
        #plt.grid()
        #vertical and horizontal lines
        line_alpha = 1
        #plt.axhline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        #plt.axhline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        #plt.axvline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        #plt.axvline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        #sideplots
        divider = make_axes_locatable(subplot)
        #x
        axCorrx = divider.append_axes('top', 0.75, pad=0.0, sharex=subplot)
        axCorrx.autoscale(False)
        axCorrx.set_adjustable('box-forced')
        plt.setp(axCorrx.get_xticklabels(), visible=False)
        axCorrx.get_yaxis().set_visible(False)
        #axCorrx.axhline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        #axCorrx.axhline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        #axCorrx.axvline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        #axCorrx.axvline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorrx.grid()
        x_ax_int = zi.sum(axis = 0)
        x_ax_int = (x_ax_int - min(x_ax_int)) / (max(x_ax_int) - min(x_ax_int))
        axCorrx.set_ylim([-0.1,1.1])
        axCorrx.plot(xi, x_ax_int, color = 'b', lw = 1.5)
        #x_abs
        s_xi = abs_data.axes[0].points
        s_zi = abs_data.channels[0].values
        xlim =  axCorrx.get_xlim()
        min_index = np.argmin(abs(s_xi - min(xlim)))
        max_index = np.argmin(abs(s_xi - max(xlim)))
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi - min(s_zi_in_range)
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi / max(s_zi_in_range) 
        axCorrx.plot(s_xi, s_zi, lw = 1.5, color = 'k')
        #y
        axCorry = divider.append_axes('right', 0.75, pad=0.0, sharey=subplot)
        axCorry.autoscale(False)
        axCorry.set_adjustable('box-forced')
        plt.setp(axCorry.get_yticklabels(), visible=False)
        axCorry.get_xaxis().set_visible(False)
        #axCorry.axhline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        #axCorry.axhline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        #axCorry.axvline(A_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        #axCorry.axvline(B_eV, color = 'k', alpha = line_alpha, lw = 1, ls = '--')
        axCorry.grid()      
        y_ax_int = zi.sum(axis = 1)
        y_ax_int = (y_ax_int - min(y_ax_int)) / (max(y_ax_int) - min(y_ax_int))
        axCorry.set_xlim([-0.1,1.1])
        axCorry.plot(y_ax_int, yi, color = 'b', lw = 1.5)
        #y_abs
        s_xi = abs_data.axes[0].points
        s_zi = abs_data.channels[0].values
        xlim =  axCorry.get_ylim()
        min_index = np.argmin(abs(s_xi - min(xlim)))
        max_index = np.argmin(abs(s_xi - max(xlim)))
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi - min(s_zi_in_range)
        s_zi_in_range = s_zi[min(min_index, max_index):max(min_index, max_index)]
        s_zi = s_zi / max(s_zi_in_range) 
        axCorry.plot(s_zi, s_xi, lw = 1.5, color = 'k')
    
    #colorbar
    plt.subplot(ax_cb)
    plt.colorbar(cax, cax = ax_cb, ticks=np.linspace(0, 1, 11))
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)
    
    plt.figtext(0.01, 0.5, '$\mathsf{\hslash \omega_2 \,\, (eV)}$', rotation = 90, size = 30)
    plt.figtext(0.4, 0.02, '$\mathsf{\hslash \omega_1 \, = \, \hslash \omega_m \,\, (eV)}$', size = 30)
    
    
    #plt.savefig('data workup grant.png', transparent = True, dpi = 300)
    #plt.close()







