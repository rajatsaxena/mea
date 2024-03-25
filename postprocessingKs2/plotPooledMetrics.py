# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 06:25:03 2023

@author: jshobe
"""

import os
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
warnings.simplefilter('ignore')

# read csv with start and end time for each experimental animal
dirname = r'T:\SWIL-Rajat\Spikesorted-SWIL'
epochsfname = 'swil-animals.csv'
epochsdf = pd.read_csv(os.path.join(dirname,epochsfname))
filename = epochsdf['file_name']
fs = 30000.0

# output variables
troughToPeak = []
duration = []
waveforms = []

# loop through all directories
for i,fname in enumerate(filename):
    print('.........')
    print('Processing ' + str(fname))
    aname = fname
    
    # load waveform properties
    df = pd.read_csv(os.path.join(dirname,'analyzedMetrics',aname+'-UnitMetrics.csv'), index_col=0)
    print(aname, len(df))
    df.insert(0, 'aname', aname)
    if i==0:
        maindf = df
    else:
        maindf = pd.concat([maindf, df])
    peak_channel = list(np.array(df['ch'], dtype='int'))
    duration.extend(df['duration'])
    troughToPeak.extend(df['wftroughToPeak'])
    
    # load and save raw waveforms
    wf = np.load(os.path.join(dirname, fname, 'proc-waveforms.npy'))
    wfmean = []
    for w, p in zip(wf,peak_channel):
        wfmean.append(w[p,:])
    wfmean = np.array(wfmean)
    if i==0:
        waveforms = wfmean
    else:
        waveforms = np.concatenate((waveforms, wfmean))
duration = np.array(duration)
troughToPeak = np.array(troughToPeak)
waveforms = np.array(waveforms)
maindf = maindf.reset_index(drop=True)
maindf = maindf.loc[:, ~maindf.columns.str.contains('^Unnamed')]
del wf, wfmean, w

# detect positive waveforms and update cell types
max_positive = np.max(waveforms*(waveforms>0), axis=1, initial=0)  # Maximum positive value in each row
min_negative = np.min(waveforms*(waveforms<0), axis=1, initial=0)  # Maximum negative value in each row
polarity = np.where(max_positive > -min_negative, np.abs(min_negative)/max_positive, -np.abs(min_negative)/max_positive)
maindf['polarityR'] = polarity
maindf.loc[(maindf.polarityR>0) & (maindf.polarityR<0.5), 'cellType'] = 'Positive'

maindf.loc[(maindf.cellType!='Positive') & (maindf.cellType!='Biphasic') & (maindf.wftroughToPeak<=0.425), 'cellType'] = 'Narrow Interneuron'
maindf.loc[(maindf.cellType!='Positive') & (maindf.cellType!='Biphasic') & (maindf.wftroughToPeak>0.425) &  (maindf.acgRiseTime>6), 'cellType'] = 'Wide Interneuron'
maindf.loc[(maindf.cellType!='Positive') & (maindf.cellType!='Biphasic') & (maindf.wftroughToPeak>0.425) &  (maindf.acgRiseTime<=6), 'cellType'] = 'Pyramidal Cell'

# save metrics pooled across animals after rearranging
maindf = maindf[['aname', 'cluster_id', 'ch', 'depth', 'cellType', 
                 'num_spikes', 'firing_rate', 'amp', 'amp_cutoff', 
                 'isi_viol', 'presence_ratio', 'duration', 'halfwidth', 
                 'PT_ratio', 'repolarization_slope', 'recovery_slope',
                 'wftroughToPeak', 'acgRiseTime', 'wfampCE',  'wfpeakToTrough', 
                 'wfabRatio', 'wfpolarity', 'wfderivative', 'wfpeakA', 'wfpeakB', 
                 'wftrough', 'polarityR', 'tmi', 'burstRoyer', 'burstDoublets']]
maindf.to_csv(os.path.join(dirname,'analyzedMetrics','pooledMetricsAll.csv'))
print(maindf.cellType.value_counts())

# plot data with waveform properties
figname = os.path.join(dirname,'analyzedMetrics','Fig2a.pdf')
colors = {'Narrow Interneuron': 'green', 'Pyramidal Cell': 'purple', 'Wide Interneuron': 'purple'}   # Colors for each cell type
df = maindf[(maindf.cellType!='Biphasic') & (maindf.cellType!='Positive')]
df.loc[df.cellType=='Wide Interneuron', 'cellType'] = 'Pyramidal Cell'
fig = plt.figure(figsize=(10,6))
sns.histplot(data=df, x="wftroughToPeak", hue='cellType', palette=colors, 
             stat='probability', bins=np.arange(0.075, 1.05, 0.025))
plt.axvline(x=0.425, c='k', linewidth=3, linestyle='--')
#sns.scatterplot(x='wftroughToPeak', y='firing_rate', data=df, hue='cellType', 
#                palette=colors, style='cellType', markers=['^', 'o'], 
#                legend=None)
plt.xlabel('Trough To Peak (ms)', fontsize=22)
plt.ylabel('Probability', fontsize=22)
plt.yticks(fontsize=18)
plt.xticks(fontsize=18)
sns.despine()
plt.tight_layout()
plt.show()