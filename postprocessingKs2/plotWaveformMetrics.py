# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 06:25:03 2023

@author: jshobe
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
warnings.simplefilter('ignore')

# read csv with start and end time for each experimental animal
dirname = r'.\Spikesorted-SWIL'
epochsfname = 'swil-animals.csv'
epochsdf = pd.read_csv(os.path.join(dirname,epochsfname))
filename = epochsdf['file_name']
fs = 30000.0

# waveform properties
troughToPeak = []
acgRiseTime = []
cellType = []
firingRate= []
duration = []
waveforms = []

# loop through all directories
for i,fname in enumerate(filename):
    print('.........')
    print('Processing ' + str(fname))
    fname = os.path.join(dirname,fname)
    
    # load waveform properties
    df = pd.read_csv(os.path.join(fname,'proc-metrics.csv'))
    idx = np.where(df['isGood'])[0]
    df = df[df['isGood']]
    peak_channel = list(df['ch'])
    troughToPeak.extend(df['troughToPeak'])
    acgRiseTime.extend(df['acgRiseTime'])
    firingRate.extend(df['firing_rate'])
    ctype = np.array(df['cellType'])
    ctype[ctype=='Pyramidal Cell']=0
    ctype[ctype=='Wide Interneuron']=1
    ctype[ctype=='Narrow Interneuron']=2
    cellType.extend(list(ctype))
    
    # load raw waveforms
    channel_map = np.load(os.path.join(fname, 'channel_map.npy'), allow_pickle=True)[0]
    wf = np.load(os.path.join(fname, 'proc-waveforms.npy'))
    wf = wf[idx,peak_channel,:]
    if i==0:
        waveforms = wf
    else:
        waveforms = np.concatenate((waveforms, wf))
#    for j,w in enumerate(wf):
#        plt.plot(np.linspace(0,82000*1./fs, 82), w)
#        s = 'tp: ' + str(round(troughToPeak[j],1))
#        plt.title(s)
#        plt.waitforbuttonpress(0)
#        plt.close()
        
troughToPeak = np.array(troughToPeak)
acgRiseTime = np.array(acgRiseTime)
cellType = np.array(cellType)
firingRate = np.array(firingRate)
waveforms = np.array(waveforms)

# plot the data
# Mapping different cell types to markers and colors
markers = {2: 'o', 0: '^', 1: '^'}  # Marker types for each cell type
colors = {2: 'b', 0: 'r', 1: 'r'}   # Colors for each cell type

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i, cell_type in enumerate(cellType):
    ax.scatter(troughToPeak[i], acgRiseTime[i], firingRate[i], c=colors[cell_type], marker=markers[cell_type])
ax.set_xlabel('Trough To Peak (ms)')
ax.set_ylabel('ACG Rise Time (ms)')
ax.set_zlabel('Mean Firing Rate (Hz)')
plt.title('3D Scatter Plot')
plt.show()


