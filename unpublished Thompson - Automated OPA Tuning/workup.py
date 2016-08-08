'''
First Created 2016/05/06 by Blaise Thompson

Last Edited 2016/08/08 by Blaise Thompson

Contributors: Blaise Thompson
'''


### import ####################################################################


import os
import sys
import importlib
import collections

import WrightTools as wt


### define ####################################################################


# paths
directory = os.path.dirname(__file__)
key = os.path.basename(directory)
package_folder = os.path.dirname(directory)

# shared module
spec = importlib.util.spec_from_file_location('shared', os.path.join(package_folder, 'shared.py'))
shared_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(shared_module)

# dictionaries to fill
raw_dictionary = collections.OrderedDict()
processed_dictionary = collections.OrderedDict()


### download ##################################################################


bypass_download = False

if __name__ == '__main__' and not bypass_download:
    shared_module.download(key, directory)


if False:
    ### TOPAS-C amplitude and center (preamp) #####################################
    
    
    # raw and processed are identical in this case
    out_path = 'TOPAS_C_full_preamp.p'
    
    force_workup = False
    
    def workup():
        # ensure that user wants to spend the time doing the workup
        if not force_workup:
            prompt = 'TOPAS-C amplitude and center (preamp) workup may take some time, proceed?'
            response = raw_input(prompt)
            proceed = util.strtobool(response)
            if not proceed:
                return None, None
        # get path
        motortune_path = os.path.join(directory, 'TOPAS-C', 'MOTORTUNE [w1_Crystal_1, w1_Delay_1, wa] 2016.01.13 19_00_00.data')
        # read information from headers
        headers = wt.kit.read_headers(motortune_path)
        wa_index = headers['name'].index('wa')
        zi_index = headers['name'].index('array')
        c1 = np.array(headers['w1_Crystal_1 points'])
        d1 = np.array(headers['w1_Delay_1 points'])
        # TODO: check if my line count is correct (am I thinking floats?)
        # this array is large (~2.6 billion lines)
        # it cannot be imported directly into memory
        # instead I load chunks and fit them to Gaussians as I go
        acqns = c1.size * d1.size
        outs = np.full((acqns, 4), np.nan)
        file_slicer = wt.kit.FileSlicer(motortune_path)
        function = wt.fit.Gaussian()
        for i in range(acqns):
            # get data from file
            lines = file_slicer.get(256)
            arr = np.array([np.fromstring(line, sep='\t') for line in lines]).T
            # fit data, record
            out = function.fit(arr[zi_index], arr[wa_index])
            outs[i] = out
            wt.kit.update_progress(100.*i/10201)
        outs.shape = (c1.size, d1.size, 4)
        cen = outs[..., 0].T
        wid = outs[..., 1].T
        amp = outs[..., 2].T
        bas = outs[..., 3].T
        # assemble data object
        c1_axis = wt.data.Axis(c1, units=None, name='c1')
        d1_axis = wt.data.Axis(d1, units=None, name='d1')
        axes = [c1_axis, d1_axis]
        cen_channel = wt.data.Channel(cen, units='nm', name='center')
        wid_channel = wt.data.Channel(wid, units='wn', name='width')
        amp_channel = wt.data.Channel(amp, units=None, name='amplitude')
        bas_channel = wt.data.Channel(bas, units=None, name='baseline')
        channels = [amp_channel, cen_channel, wid_channel, bas_channel]
        data = wt.data.Data(axes, channels, name='TOPAS-C preamp')
        data.save(os.path.join(directory, out_path))
        # finish
        return data, data.copy()
          
    # get from pickle or create
    if os.path.isfile(os.path.join(directory, out_path)) and not force_workup:
        raw_data = wt.data.from_pickle(os.path.join(directory, out_path), verbose=False)
        processed_data = raw_data.copy()
    else:
        raw_data, processed_data = workup()
    
    # check version
    if raw_data is None:
        pass
    elif not raw_data.__version__.split('.')[0] == wt.__version__.split('.')[0]:
        raw_data, processed_data = workup()
    
    # add to dictionaries
    raw_dictionary['TOPAS-C preamp'] = raw_data
    processed_dictionary['TOPAS-C preamp'] = processed_data
    
    
    ### TOPAS-C amplitude and center (poweramp) ###################################
    
    
    # raw and processed are identical in this case
    out_path = 'TOPAS_C_full_poweramp_guassian.p'
    
    force_workup = False
    
    def workup():
        # ensure that user wants to spend the time doing the workup
        if not force_workup:
            prompt = 'TOPAS-C amplitude and center (poweramp) workup may take some time, proceed?'
            response = raw_input(prompt)
            proceed = util.strtobool(response)
            if not proceed:
                return None, None
        # get path
        motortune_path = os.path.join(directory, 'TOPAS-C', 'MOTORTUNE [w1, w1_Crystal_2, w1_Delay_2, wa] 2016.01.25 16_56_06.data')
        # for some reason, read headers fails for this file...
        # define information that would normally be contained in headers manually
        wa_index = 28
        zi_index = 29
        w1 = np.linspace(1140, 1620, 25)
        c2 = np.linspace(-2.5, 2.5, 51)  # TODO: actual values
        d2 = np.linspace(-1.5, 1.5, 51)  # TODO: actual values
        # this array is large (~16 million lines)
        # it cannot be imported directly into memory
        # instead I load chunks and fit them to Gaussians as I go
        acqns = w1.size * c2.size * d2.size
        outs = np.full((acqns, 4), np.nan)
        file_slicer = wt.kit.FileSlicer(motortune_path)
        function = wt.fit.Gaussian()
        for i in range(acqns):
            # get data from file
            lines = file_slicer.get(256)
            arr = np.array([np.fromstring(line, sep='\t') for line in lines]).T
            # fit data, record
            out = function.fit(arr[zi_index], arr[wa_index])
            outs[i] = out
            wt.kit.update_progress(100.*i/acqns)
        outs.shape = (w1.size, c2.size, d2.size, 4)
        cen = outs[..., 0]
        wid = outs[..., 1]
        amp = outs[..., 2]
        bas = outs[..., 3]
        # assemble data object
        w1_axis = wt.data.Axis(w1, units='nm', name='w1')
        c1_axis = wt.data.Axis(c2, units=None, name='c2')
        d1_axis = wt.data.Axis(d2, units=None, name='d2')
        axes = [w1_axis, c1_axis, d1_axis]
        cen_channel = wt.data.Channel(cen, units='nm', name='center')
        wid_channel = wt.data.Channel(wid, units='wn', name='width')
        amp_channel = wt.data.Channel(amp, units=None, name='amplitude')
        bas_channel = wt.data.Channel(bas, units=None, name='baseline')
        channels = [amp_channel, cen_channel, wid_channel, bas_channel]
        data = wt.data.Data(axes, channels, name='TOPAS-C poweramp')
        data.save(os.path.join(directory, out_path))
        # finish
        return data, data.copy()
    
    # get from pickle or create
    if os.path.isfile(os.path.join(directory, out_path)) and not force_workup:
        raw_data = wt.data.from_pickle(os.path.join(directory, out_path), verbose=False)
        processed_data = raw_data.copy()
    else:
        raw_data, processed_data = workup()
    
    # check version
    if raw_data is None:
        pass
    elif not raw_data.__version__.split('.')[0] == wt.__version__.split('.')[0]:
        raw_data, processed_data = workup()
    
    # add to dictionaries
    raw_dictionary['TOPAS-C poweramp'] = raw_data
    processed_dictionary['TOPAS-C poweramp'] = processed_data
    
    
    ### TOPAS-C amplitude and center (poweramp, moments) ##########################
    
    
    # raw and processed are identical in this case
    out_path = 'TOPAS_C_full_poweramp_moments.p'
    
    force_workup = False
    
    def workup():
        # ensure that user wants to spend the time doing the workup
        if not force_workup:
            prompt = 'TOPAS-C amplitude and center (poweramp) workup may take some time, proceed?'
            response = raw_input(prompt)
            proceed = util.strtobool(response)
            if not proceed:
                return None, None
        # get path
        motortune_path = os.path.join(directory, 'TOPAS-C', 'MOTORTUNE [w1, w1_Crystal_2, w1_Delay_2, wa] 2016.01.25 16_56_06.data')
        # for some reason, read headers fails for this file...
        # define information that would normally be contained in headers manually
        wa_index = 28
        zi_index = 29
        w1 = np.linspace(1140, 1620, 25)
        c2 = np.linspace(-2.5, 2.5, 51)
        d2 = np.linspace(-1.5, 1.5, 51)
        # this array is large (~16 million lines)
        # it cannot be imported directly into memory
        # instead I load chunks and fit them to Gaussians as I go
        acqns = w1.size * c2.size * d2.size
        outs = np.full((acqns, 6), np.nan)
        file_slicer = wt.kit.FileSlicer(motortune_path)
        function = wt.fit.Moments()
        for i in range(acqns):
            # get data from file
            lines = file_slicer.get(256)
            arr = np.array([np.fromstring(line, sep='\t') for line in lines]).T
            # fit data, record
            out = function.fit(arr[zi_index], arr[wa_index])
            outs[i] = out
            wt.kit.update_progress(100.*i/acqns)
        file_slicer.close()
        outs.shape = (w1.size, c2.size, d2.size, 6)
        # assemble data object
        w1_axis = wt.data.Axis(w1, units='nm', name='w1')
        c1_axis = wt.data.Axis(c2, units=None, name='c2')
        d1_axis = wt.data.Axis(d2, units=None, name='d2')
        axes = [w1_axis, c1_axis, d1_axis]
        ch_0 = wt.data.Channel(outs[..., 0], units='nm', name='integral')
        ch_1 = wt.data.Channel(outs[..., 1], units='wn', name='one')
        ch_2 = wt.data.Channel(outs[..., 2], units=None, name='two')
        ch_3 = wt.data.Channel(outs[..., 3], units=None, name='three')
        ch_4 = wt.data.Channel(outs[..., 4], units=None, name='four')
        ch_5 = wt.data.Channel(outs[..., 5], units=None, name='baseline')
        channels = [ch_0, ch_1, ch_2, ch_3, ch_4, ch_5]
        data = wt.data.Data(axes, channels, name='TOPAS-C poweramp')
        data.save(os.path.join(directory, out_path))
        # finish
        return data, data.copy()
    
    # get from pickle or create
    if os.path.isfile(os.path.join(directory, out_path)) and not force_workup:
        raw_data = wt.data.from_pickle(os.path.join(directory, out_path), verbose=False)
        processed_data = raw_data.copy()
    else:
        raw_data, processed_data = workup()
    
    # check version
    if raw_data is None:
        pass
    elif not raw_data.__version__.split('.')[0] == wt.__version__.split('.')[0]:
        raw_data, processed_data = workup()
    
    # add to dictionaries
    raw_dictionary['TOPAS-C poweramp moments'] = raw_data
    processed_dictionary['TOPAS-C poweramp moments'] = processed_data
    
    
    ### OPA800 signal and idler motortune #########################################
    
    
    if False:
        # TODO: make and save data pickle here
        # TODO: seperate figure making from data import
        data_path = r'MOTORTUNE [w2_Grating, w2_BBO] 2015.10.15 17_38_22.data'
        # this data is sufficiently old that we have to process it manually :-(
        # get values from file
        headers = wt.kit.read_headers(data_path)
        arr = np.genfromtxt(data_path).T
        # extract arrays
        grating_index = headers['name'].index('w2_Grating')
        bbo_index = headers['name'].index('w2_BBO')
        signal_index = headers['name'].index('pyro2_mean')
        gra = arr[grating_index]
        gra.shape = (-1, 401)
        gra = gra[:, 0]
        bbo = arr[bbo_index]
        bbo.shape = (-1, 401)
        bbo = bbo[0]
        sig = arr[signal_index]
        sig.shape = (-1, 401)
        sig -= sig.min()
        sig /= sig.max()
        sig = sig.T
        # prepare plot
        fig = plt.figure(figsize=[8, 6])
        gs = grd.GridSpec(1, 2, hspace=0.05, wspace=0.05, width_ratios=[20, 1])
        # pcolor
        cmap = wt.artists.colormaps['default']
        ax = plt.subplot(gs[0])
        X, Y, Z = wt.artists.pcolor_helper(gra, bbo, sig)    
        cax = plt.pcolor(X, Y, Z, vmin=0, vmax=np.nanmax(Z), cmap=cmap)
        ax.set_xlim(22, gra.max())
        ax.set_ylim(bbo.min(), bbo.max())
        # labels
        ax.set_xlabel('Grating (mm)', fontsize=16)
        ax.set_ylabel('BBO (mm)', fontsize=16)
        ax.grid()
        ax.axvline(34.6, c='k', alpha=0.5, lw=2)
        ax.axhline(39.2, c='k', alpha=0.5, lw=2)
        # on-plot labels
        distance = 0.05
        wt.artists.corner_text('II signal', ax=ax, corner='UL', fontsize=16, distance=distance)
        wt.artists.corner_text('II idler', ax=ax, corner='UR', fontsize=16, distance=distance)
        wt.artists.corner_text('III signal', ax=ax, corner='LL', fontsize=16, distance=distance)
        wt.artists.corner_text('III idler', ax=ax, corner='LR', fontsize=16, distance=distance)
        # colorbar
        plt.colorbar(cax, cax=plt.subplot(gs[1]))
        # finish
        plt.savefig('signal_and_idler_motortune.png', dpi=300, transparent=True)
        
        
    ### OPA800 DFG Mixer Motortune ################################################
    
    
    if False:
        # TODO: make and save data pickle here
        # TODO: seperate figure making from data import
        data_path = r'MOTORTUNE [w2, w2_Mixer] 2015.10.16 17_59_58.data'
        # this data is sufficiently old that we have to process it manually :-(
        # get values from file
        headers = wt.kit.read_headers(data_path)
        arr = np.genfromtxt(data_path).T
        # extract arrays
        w2_index = headers['name'].index('w2')
        mixer_index = headers['name'].index('w2_Mixer')
        signal_index = headers['name'].index('pyro2_mean')
        w2 = arr[w2_index]
        w2.shape = (-1, 501)
        w2 = w2[:, 0]
        mix = arr[mixer_index]
        mix.shape = (-1, 501)
        mix = mix[0]
        sig = arr[signal_index]
        sig.shape = (-1, 501)
        sig -= sig.min()
        sig /= sig.max()
        sig = sig.T
        # prepare plot
        fig = plt.figure(figsize=[8, 6])
        gs = grd.GridSpec(1, 2, hspace=0.05, wspace=0.05, width_ratios=[20, 1])
        # pcolor
        cmap = wt.artists.colormaps['default']
        ax = plt.subplot(gs[0])
        X, Y, Z = wt.artists.pcolor_helper(w2, mix, sig)    
        cax = plt.pcolor(X, Y, Z, vmin=0, vmax=np.nanmax(Z), cmap=cmap)
        ax.set_xlim(w2.min(), w2.max())
        ax.set_ylim(mix.min(), mix.max())
        ax.grid()
        # axis labels
        ax.set_xlabel('w2 (wn)', fontsize=16)
        ax.set_ylabel('Grating (mm)', fontsize=16)
        # colorbar
        plt.colorbar(cax, cax=plt.subplot(gs[1]))
        # finish
        plt.savefig('DFG_mixer_motortune.png', dpi=300, transparent=True)
