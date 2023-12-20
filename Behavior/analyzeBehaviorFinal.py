# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 09:49:20 2023

@author: jshobe
"""

import os, glob
import numpy as np
from behavutils import *
import matplotlib.pyplot as plt

# init variables
rewardLocs = {'1':[0.2,0.65], '2':[0.45,0.83], '28':[0.35, 0.83]}
Halltitle = {'1': 'Classroom', '2':'Sunset', '28':'Universal'}   
BIN2CM = 3.14 # cm
BINWIDTH = 0.01 # 3cm bins
fs = 30 # behavior sampling rate
radiusVR = 50 # default ball radius


# directories to work on
# NOTE the main behavior file had RECORDING.mat in the end
animalsdirname = [r'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound3\SWIL105',
                  r'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound3\SWIL11',
                  r'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound3\SWIL12',
                  r'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound3\SWIL13',
                  r'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound4\SWIL15',
                  r'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound4\SWIL18',
                  r'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound4\SWIL19',
                  r'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound4\SWIL20',
                  r'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound5\SWIL22',
                  r'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound5\SWIL23',
                  r'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound5\SWIL24',
                  r'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound5\SWIL25',
                  r'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound5\SWIL26']
animalsname = ['SWIL105','SWIL11','SWIL12','SWIL13','SWIL15','SWIL18','SWIL19','SWIL20',
               'SWIL22','SWIL23','SWIL24','SWIL25','SWIL26']


# iterate through all behavior files
for aname,adname in zip(animalsname,animalsdirname):
    # load behavior data file
    fname = glob.glob(os.path.join(adname,'Behavior','dis_vel_rew_*RECORDING.mat'))
    dat = readBehavFile(fname[0])
    dat = processBehavFile(dat, samplerate=fs, VRradius=radiusVR)
    # load intan timestamps
    # intants = np.load(os.path.join(adname,'Behavior','intanTs.npy'), mmap_mode='r')
    
    # get unique hallway number
    hallnum = np.unique(np.array(dat['hallnum'], dtype=int))
    hallnum = hallnum[~np.isnan(hallnum)]
    # variable to store final op
    allhallwaydata = {}
    # iterate over each halls
    for hnum in hallnum:
        # sanity check to skip certain hallways 
        if sum(np.array(dat['hallnum'],dtype=int)==hnum)>5:
            # create hallway specific dataframe
            allhallwaydata[hnum] = calcHallwayData(dat, hnum, binwidth=BINWIDTH)
    # plot all hallways data
    figname = os.path.join(adname,'Behavior',aname+'-behav.png')
    plotSpeedAcrossHallways(allhallwaydata, figname, aname, bin2cm=BIN2CM)
    # save the data
    opfname = os.path.join(adname,'Behavior', aname+'-behavop.npy')
    np.save(opfname, allhallwaydata)
    del allhallwaydata