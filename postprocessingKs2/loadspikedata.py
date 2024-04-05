#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  7 18:06:15 2019

@author: rajat
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# load the raw data
dirname = r'T:\SWIL-Rajat\Spikesorted-SWIL\SWIL25PPC'
spike_clusters = np.ravel(np.load(os.path.join(dirname,'spike_clusters.npy'),mmap_mode='r'))
spike_times = np.ravel(np.load(os.path.join(dirname,'spike_times.npy'),mmap_mode='r'))

# convert the channel position to microns
channel_positions = np.load(os.path.join(dirname,'channel_positions.npy'))
channel_positionsX = channel_positions[:,0]*0.001
channel_positionsY = channel_positions[:,1]
# load the channel map
channel_map = np.ravel(np.load(os.path.join(dirname,'channel_map.npy')))

# load cluster information
cluster_info = pd.read_csv(os.path.join(dirname,'cluster_info.tsv'),delimiter='\t')
cluster_id = np.array(cluster_info['cluster_id'])
cluster_cnum = np.array(cluster_info['ch'])
cluster_depth = np.array(cluster_info['depth'])*0.001
cluster_num_spike = np.array(cluster_info['n_spikes'])
cluster_amp = np.array(cluster_info['Amplitude'])

# load cluster quality data and final sorted cluster amplitude and firing rate
good_clusters_id = cluster_id[np.where(cluster_info['group']=='good')[0]]
good_cluster_info = cluster_info.loc[cluster_info['cluster_id'].isin(good_clusters_id)]
good_clusters_depth = np.array(good_cluster_info['depth']*0.001)
good_clusters_cnum = np.array(good_cluster_info['ch'])
good_cluster_amp = np.array(good_cluster_info['Amplitude'])
good_clusters_channel_index = []
for cnum in good_clusters_cnum:
    good_clusters_channel_index.append((np.where(cnum==channel_map)[0][0]))
good_clusters_xpos = np.array(channel_positionsX[good_clusters_channel_index])
good_clusters_ypos = np.array(channel_positionsY[good_clusters_channel_index]*0.001)
del cluster_info

# calculate metrics for all good units 
#unitmetrics = pd.read_csv(os.path.join(dirname, 'VR24VisUnitMetrics.csv'))
#good_cluster_id_metrics = np.array(unitmetrics['cluster_id'][unitmetrics['isGood']])

# load the spiking data and plotting it according to the depth
spiketimes = []
spikeclusterID = []
cluster_meanfr = []
cluster_xpos = []
cluster_ypos = []
cluster_amp = []
fs=30000.0
# cluPeakVoltageTH = 30.
total_time = (np.max(spike_times) - np.min(spike_times))/fs
for c, cluid in enumerate(good_clusters_id):
    spike_t = spike_times[np.where(spike_clusters==cluid)[0]]/fs
    time = spike_t[-1] - spike_t[0]
    cluster_meanfr.append(len(spike_t)/time)
    cluster_xpos.append(good_clusters_xpos[c])
    cluster_ypos.append(good_clusters_ypos[c])
    spiketimes.append(spike_t)
    spikeclusterID.append(good_clusters_id[c])
    # my_ts[c] = nap.Ts(t=spike_t, time_units='s')
spiketimes = np.array(spiketimes, dtype='object')
# tsgroup = nap.TsGroup(my_ts)
# np.save('spiketimes.npy', spiketimes)
# np.save('spikeClusterID.npy', np.array(spikeclusterID))

# ccg = nap.compute_crosscorrelogram(group=tsgroup, binsize=0.001, windowsize=100)

# bin spike count in 15ms to get the Q matrix
spike_time_bins = np.arange(np.min(spike_times)/fs,total_time, 0.05)
spike_counts = np.zeros((len(spiketimes), len(spike_time_bins)-1))
for ind, st in enumerate(spiketimes):
    spike_counts[ind,:], _ = np.histogram(st, bins=spike_time_bins)
sum_spike_counts = np.sum(spike_counts, axis=0)
# np.save('Qbehav.npy', spike_counts)
# np.save('QbehavT.npy', spike_time_bins)
del spike_times, spike_clusters

# plotting the position, depth and firing rate  
fig = plt.figure(figsize=(6,10))
plt.scatter(channel_positionsX, channel_positionsY*0.001, c='gray', marker='s', alpha=0.25)
plt.scatter(cluster_xpos+np.random.uniform(-0.01,0.01,len(cluster_xpos)), 
            cluster_ypos+np.random.uniform(-0.01,0.01,len(cluster_ypos)), 
            c=cluster_meanfr, s=18, cmap='plasma')
plt.xlabel('Channel X Position', fontsize=16)
plt.ylabel('Deep -> Sup layers (mm)', fontsize=16)
for pos in ['right','top']:
    plt.gca().spines[pos].set_visible(False)
clb = plt.colorbar()
clb.ax.set_title('Mean FR (Hz)', fontsize=12)
plt.title('N='+str(len(spiketimes)), fontsize=18)
plt.grid(False)
plt.tight_layout()
plt.show()


# plotting data
fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(12,6))
for i, spiketrain in enumerate(spiketimes):
    t = np.array(spiketimes[i])
    ax[0].scatter(t, i * np.ones_like(t), c='k', marker='|')
ax[0].set_ylabel('Cell ID', fontsize=16)
ax[0].set_title('Spike Raster', fontsize=16)
ax[1].plot(spike_time_bins[:-1], sum_spike_counts, c='k')
ax[1].set_ylabel('spike sum', fontsize=16)
ax[1].set_title('Sum Q matrix', fontsize=16)
plt.tight_layout()
plt.show()