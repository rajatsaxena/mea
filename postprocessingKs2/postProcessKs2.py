# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 19:24:33 2023

@author: saxena
"""
import os
import warnings
import numpy as np
import pandas as pd
import utilsClusterQual as ucq
import utilsWaveformMetrics as uwm
warnings.simplefilter('ignore')

# other init params
vertical_site_spacing = 25e-6 # in meters
bit_volts = 0.195
fs = 30000.0

# parameters for waveform metrics calculation
params = {}
params['samples_per_spike'] = 82
params['pre_samples'] = 20
params['spikes_per_epoch'] = 1000
params['upsampling_factor'] = 200/82
params['spread_threshold'] = 0.12

# parameters for cluster quality metrics calculation
qualParams = {} 
qualParams['isi_threshold']=0.002
qualParams['min_isi']=0.0001
qualParams['isi_viol_th'] = 0.25 # 25% contamination
qualParams['presence_ratio'] = 0.7 # firing in 70% of behavior epochs
qualParams['firing_rate_th'] = 0.005  # firing rate above 0.005 Hz
qualParams['amp_cutoff_th'] = 0.2 # missing only 20% of spikes below threshold
qualParams['amp_th'] = 30 # ampltiude greater than 30 uV 

# loop through all directories
dirnames = [r'.\Spikesorted-SWIL\SWIL8HC']
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
    waveforms, wfMetrics = uwm.extract_waveforms(data, spike_times, spike_clusters, templates, cluster_info['cluster_id'], cluster_info['ch'], channel_map, bit_volts, fs, vertical_site_spacing, params)
    wfMetrics = wfMetrics.reset_index(drop=True)
    print('Finished loading waveform metrics..')
    
    # load cluster quality metrics
    cluMetrics, spiketimesGood, spikeclustersGood = ucq.calcQualityMetrics(dname, wfMetrics['amp'], params=qualParams)
    print('Finished loading cluster quality metrics..')
    
    # concatenate data
    metrics = pd.merge(wfMetrics, cluMetrics, on="cluster_id", how="left")
    del cluMetrics, wfMetrics