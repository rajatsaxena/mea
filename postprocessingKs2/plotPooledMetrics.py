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
dirname = r'.\Spikesorted-SWIL'
epochsfname = 'swil-animals.csv'
epochsdf = pd.read_csv(os.path.join(dirname,epochsfname))
filename = epochsdf['file_name']
fs = 30000.0

# waveform properties
troughToPeak = []
cellType = []
duration = []
waveforms = []

# loop through all directories
for i,fname in enumerate(filename):
    print('.........')
    print('Processing ' + str(fname))
    aname = fname
    
    # load waveform properties
    df = pd.read_csv(os.path.join(dirname,'analyzedMetrics',aname+'-UnitMetrics.csv'))
    df.insert(0, 'aname', aname)
    if i==0:
        maindf = df
    else:
        maindf = pd.concat([maindf, df])
    peak_channel = list(np.array(df['ch'], dtype='int'))
    duration.extend(df['duration'])
    troughToPeak.extend(df['troughToPeak'])
    ctype = np.array(df['cellType'])
    ctype[ctype=='Pyramidal Cell']=0
    ctype[ctype=='Wide Interneuron']=1
    ctype[ctype=='Narrow Interneuron']=2
    cellType.extend(list(ctype))
    
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
cellType = np.array(cellType)
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


# detect triphasic/ biphasic waveforms and update cell types
#peakA = np.max(waveforms[:,:20],1)
#trough = np.min(waveforms[:,14:58],1)
#ratio = peakA/np.abs(trough)
#idx = np.where((maindf['cellType']!='Positive') & (ratio>0.95) & (ratio<1.85))[0]
#maindf.loc[idx, 'cellType'] = 'Biphasic'
#
#
## MANUAL CURATION since both my script and CellExplorer ain't perfect
#r = maindf['duration']/maindf['troughToPeak']
#r[(r<0.77) | (r>1.5)] = 10
#r[maindf.cellType=='Positive'] = 1
#r[maindf.cellType=='Biphasic'] = 1
#idx = np.where(r==10)[0]
#indices = np.array([2732,2707,2616,2580,2546,2535,2509,2467,2464,2457,2396,
#                    2387,2375,2362,2289,2283,2232,2224,2206,2033,2031,1769,
#                    1690,1614,1433,1431,1421,1415,1404,1403,1402,1396,1393,
#                    1307,1285,1284,1282,1281,1279,1274,1272,1267,1261,1260,
#                    1254,1251,1250,1248,1243,1230,1227,1224,1222,1220,1217,
#                    1215,1213,1210,1206,1202,1201,1199,1198,1189,1186,1184,
#                    1183,1088,1001,972,920,829,775,761,753,718,716,688,686,
#                    665,662,648,642,474,473,472,470,466,464,463,459,458,455,
#                    450,446,445,443,441,440,434,431,426,425,422,415,410,409,
#                    406,405,404,399,398,395,394,393,391,390,387,386,385,384,
#                    383,382,380,378,377,375,371,370,369,368,367,365,356,350,
#                    348,342,340,338,337,335,334,332,327,325,319,299,286,215,
#                    214,148,136,129,3])
#maindf.loc[indices, 'troughToPeak'] = maindf.loc[indices, 'duration']
#maindf.loc[2231, 'troughToPeak'] = 0.401
#maindf.loc[1956, 'troughToPeak'] = 0.42
#maindf.loc[1734, 'troughToPeak'] = 0.63
#maindf.loc[1725, 'troughToPeak'] = 0.62
#maindf.loc[1690, 'troughToPeak'] = 0.42
#maindf.loc[1258, 'troughToPeak'] = 0.41
#maindf.loc[1229, 'troughToPeak'] = 0.41
#maindf.loc[1197, 'troughToPeak'] = 0.67
#maindf.loc[1168, 'troughToPeak'] = 0.72
#maindf.loc[1087, 'troughToPeak'] = 0.68
#maindf.loc[1078, 'troughToPeak'] = 0.64
#maindf.loc[952, 'troughToPeak'] = 0.7
#maindf.loc[813, 'troughToPeak'] = 0.35
#maindf.loc[798, 'troughToPeak'] = 0.37
#maindf.loc[598, 'troughToPeak'] = 0.37
#maindf.loc[583, 'troughToPeak'] = 0.42
#maindf.loc[499, 'troughToPeak'] = 0.41
#maindf.loc[461, 'troughToPeak'] = 0.55
#maindf.loc[454, 'troughToPeak'] = 0.72
#maindf.loc[442, 'troughToPeak'] = 0.46
#maindf.loc[403, 'troughToPeak'] = 0.49
#maindf.loc[388, 'troughToPeak'] = 0.62
#maindf.loc[372, 'troughToPeak'] = 0.58
#maindf.loc[343, 'troughToPeak'] = 0.51
#maindf.loc[326, 'troughToPeak'] = 0.35
#
## biphasic waveforms 
#indices = [2743,2690,2543,2524,2372,2339,2248,2229,2201,2193,2183,2154,2115,
#           2075,1977,1761,1735,1722,1692,1663,1649,1635,1634,1581,1526,1506,
#           1406,1388,1303,1273,1241,1236,1218,1188,1163,1068,1004,996,951,
#           937,847,827,804,739,738,685,659,630,609,509,438,437,435,432,417,
#           381,361,333,290,287,172]
#maindf.loc[indices, 'cellType'] = 'Biphasic'
maindf.loc[(maindf.cellType!='Positive') & (maindf.cellType!='Biphasic') & (maindf.troughToPeak<=0.425), 'cellType'] = 'FS'
maindf.loc[(maindf.cellType!='Positive') & (maindf.cellType!='Biphasic') & (maindf.troughToPeak>0.425) &  (maindf.acgRiseTime>6), 'cellType'] = 'WI'
maindf.loc[(maindf.cellType!='Positive') & (maindf.cellType!='Biphasic') & (maindf.troughToPeak>0.425) &  (maindf.acgRiseTime<=6), 'cellType'] = 'E'

# save metrics pooled across animals
maindf.to_csv(os.path.join(dirname,'analyzedMetrics','pooledMetricsAllAnimals.csv'))
print(maindf.cellType.value_counts())

# plot data with waveform properties
figname = os.path.join(dirname,'analyzedMetrics','pooledWFAnimals.pdf')
colors = {'FS': 'b', 'E': 'r', 'WI': 'r'}   # Colors for each cell type
df = maindf[(maindf.cellType!='Biphasic') & (maindf.cellType!='Positive')]
fig = plt.figure(figsize=(12,8))
#sns.histplot(data=df, x="troughToPeak", hue='cellType', palette=colors)
sns.scatterplot(x='troughToPeak', y='acgRiseTime', data=df, hue='cellType', 
                palette=colors, style='cellType', markers=['^', '^', 'o'], 
                legend=None)
plt.xlabel('Trough To Peak (ms)', fontsize=22)
plt.ylabel('ACG Rise Time (ms)', fontsize=22)
plt.yscale('log')
plt.yticks(fontsize=18)
plt.xticks(fontsize=18)
sns.despine()
plt.tight_layout()
plt.savefig(figname, dpi=300)
plt.close()