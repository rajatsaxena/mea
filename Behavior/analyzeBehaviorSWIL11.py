# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 09:49:20 2023

@author: jshobe
"""

import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
from behavutils import *
import scipy.stats as spst
import matplotlib.pyplot as plt

# init variables
c = {'1':[0.2,0.65], '2':[0.45,0.83], '28':[0.35, 0.83]}
Halltitle = {'1': 'Classroom', '2':'Sunset', '28':'Universal'}   
BIN2CM = 3.14 # cm
BINWIDTH = 0.01 # 3cm bins
fs = 30 # behavior sampling rate
radiusVR = 50 # default ball radius
ENVLENGTH = 314 # total env length in cm


# filenames
adname = r'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound3\SWIL11'
aname = 'SWIL11'
startidx = None
endidx = None

# load intan timestamps for two different recordings
digIn1 = np.load('Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound3\SWIL11\SWIL11_digIn1.npy', mmap_mode='r')
ts1 = np.arange(0, len(digIn1)*1./30000.0, 1./30000.0)
del digIn1

intanfilename = os.path.join(adname,aname+'_digIn1.npy')
intants1 = loadintants(intanfilename, stidx=startidx, etidx=endidx)
intanfilename = os.path.join(adname,aname+'_digIn2.npy')
intants2 = loadintants(intanfilename, stidx=startidx, etidx=endidx)
intants2 = intants2[3:] + ts1[-1] + 1./30000.0
intants = np.concatenate((intants1, intants2))

print('Processing '+str(aname)+'...')
# load behavior data file
fname = glob.glob(os.path.join(adname,'Behavior','dis_vel_rew_*RECORDING.mat'))
dat1 = readBehavFile(fname[0])
dat1 = processBehavFile(dat1, samplerate=fs, VRradius=radiusVR, intants=None)
dat1 = dat1.iloc[:len(intants1)]
fname = glob.glob(os.path.join(adname,'Behavior','dis_vel_rew_*RECORDINGSESSION2.mat'))
dat2 = readBehavFile(fname[0])
dat2 = processBehavFile(dat2, samplerate=fs, VRradius=radiusVR, intants=None)
dat2['time'] = dat2['time'] + ts1[-1] + 1./30000.0
dat2['lapnum'] = dat2['lapnum'] + dat1['lapnum'].iloc[-1]
dat = pd.concat((dat1,dat2))
dat = dat.reset_index(drop=True)
del dat1, dat2, intants1, intants2

# re calc velocity
dat['time'] = intants
vel = np.diff(dat['poscm'])/np.diff(intants)
vel = np.insert(vel,0,np.nan)
dat['vel'] = vel

# get unique hallway number
hallnum = np.unique(np.array(dat['hallnum'], dtype=int))
hallnum = hallnum[~np.isnan(hallnum)]
# variable to store final op
allhallwaydata = {}
meanSpeed = []
# iterate over each halls
for hnum in hallnum:
    # sanity check to skip certain hallways 
    if sum(np.array(dat['hallnum'],dtype=int)==hnum)>5:
        # create hallway specific dataframe
        allhallwaydata[hnum] = calcHallwayData(dat, hnum, binwidth=BINWIDTH)
        normspeed = norm1D(np.nanmean(allhallwaydata[hnum]['binned_speed'],0)[1:])
        print(allhallwaydata[hnum]['binned_speed'].shape[0])
        meanSpeed.append(normspeed)
figname = os.path.join(adname,'Behavior',aname+'-behav.png')
plotSpeedAcrossHallways(allhallwaydata, figname, aname, bin2cm=BIN2CM)
# save the data
opfname = os.path.join(adname,'Behavior', aname+'-behavop.npy')
np.save(opfname, allhallwaydata)
