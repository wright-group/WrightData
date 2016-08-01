### import ####################################################################


import os
import itertools
import collections

import matplotlib
matplotlib.rcParams['font.size'] = 14
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import matplotlib.patheffects as PathEffects
plt.close('all')

import numpy as np

import WrightTools as wt


### define ####################################################################


directory = os.path.dirname(__file__)

cmap = wt.artists.colormaps['default']


### material introduction #####################################################


output_path = os.path.join(directory, 'material introduction.png')

force_plotting = False

if not os.path.isfile(output_path) or force_plotting:
    pass


### normalization procedure ###################################################


# TODO:

if False:
    # -*- coding: utf-8 -*-
    """
    Created on Fri Jul 31 23:27:32 2015
    ​
    apply corrections to the 2d frequency scans and compose a summary image
    ​
    @author: Dan
    """

    import os
    import numpy as np
    import fscolors as f
    from scipy.interpolate import interp1d

    user = os.path.expanduser('~')
    ###     abs file import
    abs_file = '\\'.join([user,
               'Documents', 
               'Colloidal Synthesis Characterization',
               'DK-2014.02.05-PbSe 019'])
    abs_file += r'\fs19 2014.09.01.trimmed.txt'
    l,a = np.loadtxt(abs_file, unpack = True, skiprows=18)
    a -= a.min()
    A = interp1d(1e7/l, a)

    cmap = f.cm.chw2#'cubehelix' #f.Dat.mycm

    ### define M-factors
    def M(w1,w2,n):
        a1 = A(w1)
        a2 = A(w2)
        out = 10**(-a1*n/2) * (1-10**(-a2[:,None]))
        out /= a2[:,None]
        return out

    # same weights as the other workups:
    n = np.array([0.625, 0.347, 0.147, 0.0764, 0.0475])
    od = np.array([0.79, 0.43, 0.18, 0.1, 0.06])
    n = n / n[0]

    ### opa power smoothness import
    opa_folder = '\\'.join([user,
                              'Desktop',
                              'incoming from fs table',
                              'fsb 19',
                              '2014.09.14',
                              'calibrations'])
    opa1_file = opa_folder + r'\opa1 power smoothness (1) 4000 ms wait [1x51] Freq.dat'
    opa2_file = opa_folder + r'\opa2 power smoothness (1) 4000 ms wait [1x51] Freq.dat'
    opa1_dat = np.loadtxt(opa1_file)
    opa2_dat = np.loadtxt(opa2_file)
    f_opa1 = interp1d(1e7/opa1_dat[:,1], opa1_dat[:,16])
    f_opa2 = interp1d(1e7/opa1_dat[:,3], opa1_dat[:,16])
    I1 = f_opa1(np.linspace(6256,8699,num=51))
    I2 = f_opa1(np.linspace(6281,8699,num=42))


    ###     2d freq imports
    dat_folder = '\\'.join([user,
                              'Desktop',
                              'incoming from fs table',
                              'fsb 19',
                              '2014.09.19 dilution study'])
    # only four files this time; make it easy
    names = [r'\fsb 19.1 (1) 60 ms wait [42x49] Freq.dat',
             r'\fsb 19.2 (2) 60 ms wait [42x49] Freq.dat',
             r'\fsb 19.3 (1) 60 ms wait [42x49] Freq.dat',
             r'\fsb 19.4 (1) 60 ms wait [42x49] Freq.dat']
    fs = []
    import matplotlib.pyplot as plt
    plt.close('all')
    for i in range(len(names)):
        name = names[i]
        fi = f.Dat(filepath=dat_folder+name,
                   xvar='w1', yvar='w2')
        print fi.zi.shape
        fi.zi -= -0.069
        fi.znull=0.
        mi = M(fi.xi,fi.yi, n[i+1])
        fi.zi /= mi**2
        if fi.zi.shape[1] == 49:
            fi.zi /= I1[1:-1] * I2[:,None]**2
        else:
            fi.zi /= I1 * I2[:,None]**2
        plt.figure()
        fi.zfit()
        fi.smooth(x=2,y=2)
        fi.plot2d(alt_zi='amp', contour=True)
        fi.side_plots(fi.s1)
        fi.colorbar()
        fs.append(fi)

    if True:
        plt.close('all')
        plt.figure(figsize=(12,4))
        xticks = np.array([6500,7500,8500])
        yticks = xticks
        for i in range(len(fs)):
            #spi = 111
            spi = 141 + i
            plt.subplot(spi, aspect='equal')
            plt.title(r'$\mathsf{OD_{1S}=}$ ' + '{0}'.format(od[i+1]))
            fi = fs[i]
            plt.contourf(fi.xi, fi.yi, fi.zi_norm, 
                         levels=np.linspace(0,fi.zi_norm.max(),256), 
                         cmap=cmap)
            plt.contour(fi.xi, fi.yi, fi.zi_norm, 
                        levels=np.linspace(0., fi.zi_norm.max(), 11)[1:-1], 
                        colors='k', linewidths=1, alpha=0.1)#'k')
            plt.plot([min([fi.xi.max(),fi.yi.max()]),max([fi.xi.min(),fi.yi.min()])],
                     [min([fi.xi.max(),fi.yi.max()]),max([fi.xi.min(),fi.yi.min()])],
                     'k:')
            if int(str(spi)[-1]) != 1:
                plt.yticks(yticks, visible=False)
            else:
                plt.yticks(yticks)
                plt.ylabel(r'$\mathsf{\bar\nu_2 (cm^{-1})}$')
            plt.xticks(xticks)
            plt.xlabel(r'$\mathsf{\bar\nu_1 (cm^{-1})}$')
            plt.grid(b=True)

    #plt.savefig('dilution3.png', transparent=True, dpi=300)


### batch 19 ##################################################################


if False:
    # -*- coding: utf-8 -*-
    """
    Created on Wed Jul 08 14:59:28 2015
    
    movie of fsb 19 fwm
    
    @author: Dan
    """
    
    import os
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import fscolors as f
    files_folder = os.path.expanduser(r'~\Desktop\incoming from fs table\fsb 19\2014.09.14')
    
    #data = np.load(files_folder + r'\tot rank 81 NRMSD 0.00%.npz')
    data = np.load(files_folder + r'\tot rank 54 NRMSD 1 0.00%.npz')
    # d2, w1, w2, zis[d2, w1, w2]
    d2 = data['d2']
    w1 = data['w1']
    w2 = data['w2']
    zis = data['zi']
    
    mycm = f.cm.chw2
    log = False
    cutoff = -4.5
    
    
    fig1 = plt.figure(figsize=(2*3.5,2*3))
    gs = matplotlib.gridspec.GridSpec(3,4, width_ratios=[10,10,10,1.25])
    cb = plt.subplot(gs[:,-1])
    plt.rcParams.update({
        'font.size':12
    })
    
    xticks = np.linspace(6500, 8500, num=3)
    yticks = xticks
    
    zmax = zis.max()
    zis /= zmax
    if log:
        levels = np.linspace(cutoff, 0, num=256)
    else:
        levels = np.linspace(0, 1, num=256)
    
    iss = [8,10,12,
           14,16,18,
           20,22,-1
           ]
    dss = [75, 50, 25, 
           0, -25, -50, 
           -75, -100, -800
           ]
    
    for i in range(9):
        rowi = i / 3
        coli = i % 3
        di = iss[i]
        sij = plt.subplot(gs[rowi,coli], aspect='equal')
        if rowi==1 and coli==0:
            plt.ylabel(r'$\mathsf{\bar\nu_2 (cm^{-1})}$', fontsize=16)
        if rowi==2 and coli==1:
            plt.xlabel(r'$\mathsf{\bar\nu_1 (cm^{-1})}$', fontsize=16)
       # need to mask negatives
        if log:
            zi = np.ma.masked_less(zis[di].T, 10**cutoff)
            zi = np.log10(zi.filled(10**cutoff))
        else:
            zi = np.ma.masked_less(zis[di].T, 0.)
            zi = np.sqrt(zi.filled(0.))
        if log:
            cb_ax = plt.contourf(w1, w2, zi/2., levels/2., 
                                 cmap=mycm)
            plt.contour(w1, w2, zi/2., 
                        levels=np.linspace(levels.min()-1e-6, levels.max(),10)/2, 
                        colors='k', linewidths=1, alpha=0.2)
        else:
            cb_ax = plt.contourf(w1, w2, zi, levels, 
                                 cmap=mycm)
            plt.contour(w1, w2, zi, 
                        levels=np.linspace(levels.min()-1e-6, levels.max(),10), 
                        colors='k', linewidths=1, alpha=0.2)
        plt.title(r'$\mathsf{\tau_{21}=}$'+r' {0} fs'.format(int(dss[i])), 
                  fontsize=12)
        if coli > 0:
            plt.yticks(yticks, visible=False)
        else:
            plt.yticks(yticks)
        if rowi != 2:
            plt.xticks(xticks, visible=False)
        else:
            plt.xticks(xticks)
        # draw lines for expected resonances?
        plt.grid()
    if log:
        cbx = fig1.colorbar(cb_ax, cax = cb, 
                            ticks=np.linspace(levels.min(),levels.max(),num=10)/2.)
    else:
        cbx = fig1.colorbar(cb_ax, cax = cb, 
                            ticks=np.linspace(levels.min(),levels.max(),num=11))
    #cbx.ax.yaxis.set_ticks_position('left')
    plt.tight_layout()
    #plt.savefig('fsb 25.png', transparent=True, dpi=300)



### 2D Frequency vs OD ########################################################

if False:
    # TODO:
    
    output_path = os.path.join(directory, 'WMELs.png')
    
    force_plotting = False
    
    if not os.path.isfile(output_path) or force_plotting:
        # prepare for plot
        fig = plt.figure(figsize=[6.5, 6])
        gs = grd.GridSpec(3, 8, width_ratios=[1.5, 1, 1, 1, 1, 1, 1, 1])
        # large levels figure
        ax = plt.subplot(gs[:, 0])
        # plot energies
        energies = list(np.linspace(0, 0.85, 4)) + [1]
        print energies
        for energy in energies:
            ax.axhline(energy, color='k', linewidth=2, ls='-')
        # add arrow method for subplot a
        def add_arrow(between, x_pos, width):
            head_size = 0.05
            color = 'k'
            # calculate arrow length
            arrow_length = energies[between[1]] - energies[between[0]]
            arrow_end = energies[between[1]]
            width_compensation_factor = 2e-3
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
        add_arrow([1, 2], 0.25, 12)
        add_arrow([2, 3], 0.25, 10)
        add_arrow([3, 2], 0.75, 8)
        add_arrow([2, 1], 0.75, 4)
        add_arrow([1, 0], 0.75, 2)
        state_text_buffer = 0.25
        state_names = ['n=0', 'n=1', 'n=2', 'n=3', 'n=N']
        for i in range(len(energies)):
            ax.text(-state_text_buffer, energies[i], state_names[i], fontsize=18, verticalalignment='center', horizontalalignment ='right')
        ax.text(-state_text_buffer, 0.925, '...', rotation=-90, fontsize=18, verticalalignment='center', horizontalalignment ='right')
        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-0.025, 1.025)
        ax.axis('off')
        # highlights
        import matplotlib.patches as patches
        patch1 = patches.Rectangle((0.3525, 0.1775), 0.245, 0.175, transform=fig.transFigure, figure=fig, facecolor='y', alpha=0.25)
        patch2 = patches.Rectangle((0.689, 0.6476), 0.16, 0.175, transform=fig.transFigure, figure=fig, facecolor='y', alpha=0.25)
        fig.patches.extend([patch1, patch2])
        # pathway 1 ---------------------------------------------------------------
        energies = [0., 0.5, 1.]
        state_text_buffer = 0.25
        state_names = ['g', 'a', '2a']
        # pathway 1 alpha
        ax = plt.subplot(gs[0, 2])
        wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='I')
        wmel.add_arrow(0, [0, 1], 'ket', '1', color='r')
        wmel.add_arrow(1, [0, 1], 'bra', '2', color='b')
        wmel.add_arrow(2, [1, 0], 'bra', '2\'', color='b')
        wmel.add_arrow(3, [1, 0], 'out', '', color='r')
        for i in range(len(energies)):
            ax.text(-state_text_buffer, energies[i], state_names[i], fontsize=12, verticalalignment='center', horizontalalignment ='right')
        # pathway 1 beta
        ax = plt.subplot(gs[1, 2])
        wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
        wmel.add_arrow(0, [0, 1], 'ket', '1', color='r')
        wmel.add_arrow(1, [0, 1], 'bra', '2', color='b')
        wmel.add_arrow(2, [1, 2], 'bra', '2\'', color='b')
        wmel.add_arrow(3, [2, 1], 'out', '', color='r')
        for i in range(len(energies)):
            ax.text(-state_text_buffer, energies[i], state_names[i], fontsize=12, verticalalignment='center', horizontalalignment ='right')
        # pathway 1 gamma
        ax = plt.subplot(gs[2, 2])
        wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
        wmel.add_arrow(0, [0, 1], 'ket', '1', color='r')
        wmel.add_arrow(1, [1, 0], 'ket', '2', color='b')
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
        wmel.add_arrow(2, [2, 1], 'ket', '2', color='b')
        wmel.add_arrow(3, [1, 0], 'out', '', color='r')
        # pathway 2 beta
        ax = plt.subplot(gs[1, 3])
        wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
        wmel.add_arrow(0, [0, 1], 'ket', '1', color='r')
        wmel.add_arrow(1, [1, 2], 'ket', '2\'', color='b')
        wmel.add_arrow(2, [0, 1], 'bra', '2', color='b')
        wmel.add_arrow(3, [2, 1], 'out', '', color='r')
        # pathway 3 ---------------------------------------------------------------
        # pathway 3 alpha
        ax = plt.subplot(gs[0, 4])
        wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='III')
        wmel.add_arrow(0, [0, 1], 'bra', '2', color='b')
        wmel.add_arrow(1, [0, 1], 'ket', '1', color='r')
        wmel.add_arrow(2, [1, 0], 'bra', '2\'', color='b')
        wmel.add_arrow(3, [1, 0], 'out', '', color='r')
        # pathway 3 beta
        ax = plt.subplot(gs[1, 4])
        wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
        wmel.add_arrow(0, [0, 1], 'bra', '2', color='b')
        wmel.add_arrow(1, [0, 1], 'ket', '1', color='r')
        wmel.add_arrow(2, [1, 2], 'ket', '2\'', color='b')
        wmel.add_arrow(3, [2, 1], 'out', '', color='r')
        # pathway 3 gamma
        ax = plt.subplot(gs[2, 4])
        wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
        wmel.add_arrow(0, [0, 1], 'bra', '2', color='b')
        wmel.add_arrow(1, [1, 0], 'bra', '1', color='r')
        wmel.add_arrow(2, [0, 1], 'ket', '2\'', color='b')
        wmel.add_arrow(3, [1, 0], 'out', '', color='r')
        # pathway 4 ---------------------------------------------------------------
        # pathway 4 alpha
        ax = plt.subplot(gs[0, 5])
        wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='IV')
        wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
        wmel.add_arrow(1, [1, 2], 'ket', '1', color='r')
        wmel.add_arrow(2, [2, 1], 'ket', '2', color='b')
        wmel.add_arrow(3, [1, 0], 'out', '', color='r')
        # pathway 4 beta
        ax = plt.subplot(gs[1, 5])
        wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
        wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
        wmel.add_arrow(1, [1, 2], 'ket', '1', color='r')
        wmel.add_arrow(2, [0, 1], 'bra', '2', color='b')
        wmel.add_arrow(3, [2, 1], 'out', '', color='r')
        # pathway 5 ---------------------------------------------------------------
        # pathway 5 alpha
        ax = plt.subplot(gs[0, 6])
        wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='V')
        wmel.add_arrow(0, [0, 1], 'bra', '2', color='b')
        wmel.add_arrow(1, [1, 0], 'bra', '2\'', color='b')
        wmel.add_arrow(2, [0, 1], 'ket', '1', color='r')
        wmel.add_arrow(3, [1, 0], 'out', '', color='r')
        # pathway 5 beta
        ax = plt.subplot(gs[1, 6])
        wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
        wmel.add_arrow(0, [0, 1], 'bra', '2', color='b')
        wmel.add_arrow(1, [0, 1], 'ket', '2\'', color='b')
        wmel.add_arrow(2, [1, 2], 'ket', '1', color='r')
        wmel.add_arrow(3, [2, 1], 'out', '', color='r')
        # pathway 5 gamma
        ax = plt.subplot(gs[2, 6])
        wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
        wmel.add_arrow(0, [0, 1], 'bra', '2', color='b')
        wmel.add_arrow(1, [0, 1], 'ket', '2\'', color='b')
        wmel.add_arrow(2, [1, 0], 'bra', '1', color='r')
        wmel.add_arrow(3, [1, 0], 'out', '', color='r')
        # pathway 6 ---------------------------------------------------------------
        text_buffer = 1.3
        # pathway 6 alpha
        ax = plt.subplot(gs[0, 7])
        wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4, title='VI')
        wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
        wmel.add_arrow(1, [1, 0], 'ket', '2', color='b')
        wmel.add_arrow(2, [0, 1], 'ket', '1', color='r')
        wmel.add_arrow(3, [1, 0], 'out', '', color='r')
        ax.text(text_buffer, 0.5, r'$\mathrm{\alpha}$', fontsize=25, verticalalignment='center', horizontalalignment='left')
        # pathway 6 beta
        ax = plt.subplot(gs[1, 7])
        wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
        wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
        wmel.add_arrow(1, [0, 1], 'bra', '2', color='b')
        wmel.add_arrow(2, [1, 2], 'ket', '1', color='r')
        wmel.add_arrow(3, [2, 1], 'out', '', color='r')
        ax.text(text_buffer, 0.5, r'$\mathrm{\beta}$', fontsize=25, verticalalignment='center', horizontalalignment='left')
        # pathway 6 gamma
        ax = plt.subplot(gs[2, 7])
        wmel = wt.diagrams.WMEL.Subplot(ax=ax, energies=energies, interactions=4)
        wmel.add_arrow(0, [0, 1], 'ket', '2\'', color='b')
        wmel.add_arrow(1, [0, 1], 'bra', '2', color='b')
        wmel.add_arrow(2, [1, 0], 'bra', '1', color='r')
        wmel.add_arrow(3, [1, 0], 'out', '', color='r')
        ax.text(text_buffer, 0.5, r'$\mathrm{\gamma}$', fontsize=25, verticalalignment='center', horizontalalignment='left')
        # labels
        dist = 0.05
        fig.text(0.075, 1-dist, 'a', bbox=props, horizontalalignment='left', verticalalignment='top', fontsize=18, alpha=1)
        fig.text(0.33, 1-dist, 'b', bbox=props, horizontalalignment='left', verticalalignment='top', fontsize=18, alpha=1)    
        # finish
        wt.artists.subplots_adjust(fig)
        plt.savefig(output_path, dpi=300, transparent=True,  pad_inches=1.)
        plt.close(fig)
    
