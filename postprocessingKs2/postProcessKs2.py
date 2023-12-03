# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 19:24:33 2023

@author: saxena
"""
import os
import numpy as np
import pandas as pd
import utilsClusterQual as ucq
import utilsWaveformMetrics as uwm

dirnames = [r'.\Spikesorted-SWIL\SWIL8HC']
vertical_site_spacing = 25e-6 # in meters
bit_volts = 0.195
fs = 30000.0

params = {}
params['samples_per_spike'] = 82
params['pre_samples'] = 20
params['spikes_per_epoch'] = 1000
params['upsampling_factor'] = 200/82
params['spread_threshold'] = 0.12

# loop through all directories
for dname in dirnames:
    print('Processing ' + str(dname))
    # load kilosort data
    spike_times, spike_clusters, templates, channel_map, cluster_ids, cluster_info = \
        uwm.load_kilosort_data(dname, fs, convert_to_seconds = False)
    print('Finished loading kilosort data..')

    # load raw dat file data to load waveforms
    num_channels = len(channel_map)
    if os.path.isfile(os.path.join(dname, 'temp_wh.dat')):
        rawData = np.memmap(os.path.join(dname, 'temp_wh.dat'), dtype='int16', mode='r')
    else:
        bin_file = [file for file in os.listdir(dname) if file.endswith('.bin')]
        bin_file = bin_file[0]
        rawData = np.memmap(bin_file, dtype='int16', mode='r')
    data = np.reshape(rawData, (int(rawData.size/num_channels), num_channels))
    del rawData
    print('Finished loading raw binary file..')
    
    # extract waveforms metrics
    waveforms, wfmetrics = uwm.extract_waveforms(data, spike_times, spike_clusters, templates, cluster_info['cluster_id'], cluster_info['ch'], channel_map, bit_volts, fs, vertical_site_spacing, params)
    wfmetrics = wfmetrics.reset_index(drop=True)
    print('Finished loading waveform metrics..')

    # load cluster quality metrics
    cluMetrics, spiketimesGood, spikeclustersGood = ucq.calcQualityMetrics(dname, wfmetrics['amp'])
    print('Finished loading cluster quality metrics..')
    
    metrics = pd.merge(wfmetrics, cluMetrics, on="cluster_id", how="left")