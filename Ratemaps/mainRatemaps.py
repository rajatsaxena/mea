#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 19:03:43 2019

@author: rajat
"""

import os
import rmaputils 
import numpy as np
import pandas as pd
import pylab as plt
from tqdm import tqdm
import scipy.io as spio
from matplotlib.backends.backend_pdf import PdfPages

# THRESHOLD PARAMS
speedTh = 2 # cm/s
MAX_CHAN = 256

# HALLWAY information
HALLWAYS = [1,2,28] # hallways numbers, Classroom: 1, Sunset:2, Universal: 28
rewardLocs = {1:[0.2,0.65], 2:[0.45,0.83], 28:[0.35,0.83]}
colors = {1:'Blues', 2:'Oranges', 28:'Reds'}
linecolors = {1:'b', 2:'darkorange', 28:'r'}

# POS variables
fspos = 30.
dt = 1./fspos
posMin = 0
posMax = 314 #cm
binwidth = 314/100 #cm
gsmooth = 0.25
bin2cm = binwidth
xt = np.array([0,50,100])
xtcm = xt*binwidth

# LFP variables
fslfp = 1000.0
f_theta = [6,10]
amplfpth = 0.7
phasebins = np.arange(0,4*360+10,10)

# read csv with start and end time for each experimental animal
ephysdirname = r'T:\SWIL-Rajat\Spikesorted-SWIL'
behavdirname = r'T:\SWIL-Rajat\Behavior-SWIL'
lfpdirname = r'T:\SWIL-Rajat\LFP-SWIL'
rmapdirname = r'T:\SWIL-Rajat\Ratemaps-SWIL'
epochsfname = 'swil-animals.csv'
# load default channel map
rezdat = spio.loadmat(r'T:\SWIL-Rajat\analysis-scripts\postProcessingKS2\chanMap_intan256F_bottom_SM.mat')
xcoords = np.ravel(rezdat['xcoords'])
ycoords = np.ravel(rezdat['ycoords'])
channelpos = np.array(np.stack((xcoords,ycoords)).T, dtype='int16')
del xcoords,ycoords,rezdat
# load each recording epochs file
epochsdf = pd.read_csv(os.path.join(ephysdirname, epochsfname))
filename = epochsdf['file_name']
start_time, end_time, referenceHC = epochsdf['start_time'], epochsdf['end_time'], epochsdf['referenceHC']
# load pooled unit metrics for all animals
pooledMetricsdf = pd.read_csv(os.path.join(ephysdirname,'analyzedMetrics','pooledAllAnimals.csv'))
if 'Unnamed: 0' in pooledMetricsdf.columns:
    pooledMetricsdf = pooledMetricsdf.drop(['Unnamed: 0'],axis=1)

# loop through all recordings from each animal
for dname, st, et, refHCChnum in zip(filename, start_time, end_time, referenceHC):
    print('.........')
    print('Processing ' + str(dname))
    if 'PPC' in dname:
        aname = dname[:-3]
    else:
        aname = dname[:-2]
    adat = {} # save each animals data 
    
    # load behavior data
    behavdata = np.load(os.path.join(behavdirname,aname+'-behavop.npy'), allow_pickle=True).item()
    # iterate through the file to analyze data through each hallway
    for hallwaynum in HALLWAYS:
        behavd = behavdata[hallwaynum]
        
        # load the hallway trial data
        behavd = behavd['trial_data'] 
        posX, posT, posSpeed, trStart, trEnd, omaptr = rmaputils.loadOccupancy(behavd, posMin=posMin, 
                                                                                  posMax=posMax, binwidth=binwidth, 
                                                                                  speedThr=speedTh, dt=dt)
        
        # create occupancy map
        omaptrsm, omaptrnorm, omap1d, omap1dsm, omapbins = rmaputils.processOccupancy(omaptr, posX, 
                                                                                      posMin=posMin, posMax=posMax, 
                                                                                      binwidth=binwidth, smwindow=gsmooth)
        # save hallway occupancy data for each animal
        adat[hallwaynum] = {}
        adat[hallwaynum]['posX'] = posX
        adat[hallwaynum]['posT'] = posT
        adat[hallwaynum]['posSpeed'] = posSpeed
        adat[hallwaynum]['trStart'] = trStart
        adat[hallwaynum]['trEnd'] = trEnd
        adat[hallwaynum]['omaptr'] = omaptr
        adat[hallwaynum]['omaptrsm'] = omaptrsm
        adat[hallwaynum]['omaptrnorm'] = omaptrnorm
        adat[hallwaynum]['omap1d'] = omap1d
        adat[hallwaynum]['omap1dsm'] = omap1dsm
        adat[hallwaynum]['omapbins'] = omapbins
    del posX, posT, posSpeed, trStart, trEnd 
    del omaptr, omaptrsm, omaptrnorm, omap1d, omap1dsm, omapbins
    del behavd, behavdata
    
    # load spiketimes and cluster id for good unit for a given recording
    metricsAnimaldf = pooledMetricsdf[pooledMetricsdf.aname==dname]
    spiketimes, clusterId, clusterDepth, brainRegion = rmaputils.loadGoodSpikeTimes(os.path.join(ephysdirname,dname), metricsAnimaldf)
    
    # load lfp data
    lfpfilename = os.path.join(lfpdirname, aname+'_lfp.npy')
    if os.path.exists(os.path.join(lfpdirname, aname+'_lfpts.npy')):
        lfpts = np.load(os.path.join(lfpdirname, aname+'_lfpts.npy'), allow_pickle=True)
    else:
        lfpts = None
    # load theta peak times for reference hippocampus electrode
    # slm layer electrode pre-determines marked as reference
    refThetaPeaks, lfpts, lfpsig, thetapk = rmaputils.calcThetaPeaks(lfpfilename, refHCChnum, st, et, eeg_times=lfpts, f_range=f_theta, fsl=fslfp, ampth=amplfpth)
    
    # load channel indices for main map
    rezChnum, phyChnum = rmaputils.loadCorrectedChanMap(ephysdirname, dname, clusterId, metricsAnimaldf, channelpos)
        
    # iterate through spike timestamps
    spikedata = {}
    for cid, spikets, chnum, rezch in tqdm(zip(clusterId, spiketimes, metricsAnimaldf.ch, rezChnum)):
#        if 'PPC' in dname:
##            assert chnum==rezch
#            chnum = int(rezch + 256) # add to account for PPC channels being at the end
#        else:
#            chnum = int(rezch)
        # load LFP from reference hc + current unit channel
        spkPhaseRefHC = rmaputils.calcSpikePhase(refThetaPeaks, spikets)
        
        celldata = {} # dict to hold spike data for each neurons in each hallway
        for hallwaynum in HALLWAYS:
            celldata[hallwaynum] = {} 
            # calculate spike pos data
            spkPos, spkPosts, spkPosSpeed, spkPhaseHC, spkPosTrial, spkPostsTrial = rmaputils.processSpikePos(spikets, spkPhaseRefHC, adat[hallwaynum]['trStart'], 
                                                                                                              adat[hallwaynum]['trEnd'], adat[hallwaynum]['posX'], 
                                                                                                              adat[hallwaynum]['posT'], adat[hallwaynum]['posSpeed'], 
                                                                                                              lfpts)
            # calcualte spike map
            spkmaptr, spkmaptrsm, spkmaptrnorm, spkmap1d, spkmap1dsm = rmaputils.processSpikeMap(spkPosTrial, spkPos, adat[hallwaynum]['omapbins'], posMin=posMin, 
                                                                                                 posMax=posMax, binwidth=binwidth, smwindow=gsmooth)
            
            # calculate ratemaps
            rmaptr, rmaptrsm, rmaptrnorm, rmap1d, rmap1dsm, _ = rmaputils.processRateMap(spkmaptr, adat[hallwaynum]['omaptr'], spkmaptrsm, 
                                                                                                adat[hallwaynum]['omaptrsm'], spkmap1d, 
                                                                                                adat[hallwaynum]['omap1d'], spkmap1dsm, 
                                                                                                adat[hallwaynum]['omap1dsm'], fs=fspos)
            
            # firing statistics: spatial information, sparsity, stability index
            sinfo = rmaputils.calcSpatialInformationScore(rmap1dsm, adat[hallwaynum]['omap1dsm'])
            sparsity, _ = rmaputils.calcSparsity(rmap1dsm, rmaptrsm)
            stabilityInd = rmaputils.calcStabilityIndex(rmaptrsm)
            # spatialinfo significance
            sinfo_p, shuffledinfo = rmaputils.calcShuffleSpatialInfo(sinfo, spikets, adat[hallwaynum]['omap1dsm'], adat[hallwaynum]['posT'][-1], 
                                                                     adat[hallwaynum]['posT'][0], adat[hallwaynum]['trStart'], adat[hallwaynum]['trEnd'], 
                                                                     adat[hallwaynum]['posX'], adat[hallwaynum]['posT'], adat[hallwaynum]['posSpeed'],  
                                                                     adat[hallwaynum]['omapbins'], posMin=posMin, posMax=posMax, binwidth=binwidth, fs=fspos)
            # get placefield statistics: num field, field peak firing rate, size, center, dispersion
            if sinfo_p<0.01:
                pfmap, pfPeakFr, pfCenter, pfSize, pfNumFields = rmaputils.calcfieldSize(rmap1dsm, rmaptrsm, pixel2cm=bin2cm, L=posMax)
            else:
                pfmap, pfPeakFr, pfCenter, pfSize, pfNumFields = np.nan, np.nan, np.nan, np.nan, np.nan
            # calculate tmi 
            tmihc, tmicounthc, tmiedges = rmaputils.calcTMI(spkPhaseHC, bins=phasebins)
            
            # store data for each cell in each hallway
            celldata[hallwaynum]['spkPos'] = spkPos
            celldata[hallwaynum]['spkPosts'] = spkPosts
            celldata[hallwaynum]['spkPosSpeed'] = spkPosSpeed
            celldata[hallwaynum]['spkPhaseHC'] = spkPhaseHC
            celldata[hallwaynum]['spkPosTrial'] = spkPosTrial
            celldata[hallwaynum]['spkPostsTrial'] = spkPostsTrial
            celldata[hallwaynum]['spkmaptr'] = spkmaptr
            celldata[hallwaynum]['spkmaptrsm'] = spkmaptrsm
            celldata[hallwaynum]['spkmaptrnorm'] = spkmaptrnorm
            celldata[hallwaynum]['spkmap1d'] = spkmap1d
            celldata[hallwaynum]['spkmap1dsm'] = spkmap1dsm
            celldata[hallwaynum]['rmaptr'] = rmaptr
            celldata[hallwaynum]['rmaptrsm'] = rmaptrsm
            celldata[hallwaynum]['rmap1d'] = rmap1d
            celldata[hallwaynum]['rmap1dsm'] = rmap1dsm
            celldata[hallwaynum]['sinfo'] = sinfo
            celldata[hallwaynum]['sparsity'] = sparsity
            celldata[hallwaynum]['stabilityInd'] = stabilityInd
            celldata[hallwaynum]['sinfo_p'] = sinfo_p
            celldata[hallwaynum]['shuffledinfo'] = shuffledinfo
            celldata[hallwaynum]['pfPeakFr'] = pfPeakFr
            celldata[hallwaynum]['pfCenter'] = pfCenter
            celldata[hallwaynum]['pfSize'] = pfSize
            celldata[hallwaynum]['pfNumFields'] = pfNumFields
            celldata[hallwaynum]['tmihc'] = tmihc
            celldata[hallwaynum]['tmicounthc'] = tmicounthc
            celldata[hallwaynum]['phasebins'] = tmiedges
        # plot the data across hallways for this cell
        spikedata[cid] = celldata
        del celldata
    adat['spikedata'] = spikedata
    del spikedata
    # save the processed data
    opfname = os.path.join(rmapdirname, dname+'-op.npy')
    np.save(opfname, adat)
    # plot the output pdf
    pdfname = os.path.join(rmapdirname, dname+'-op.pdf')
    with PdfPages(pdfname) as pdf: 
        for cid, region, depth, spikets  in tqdm(zip(clusterId, brainRegion, clusterDepth, spiketimes)):
            fig = rmaputils.genAnalysisReport(adat, aname, spikets, HALLWAYS, cid, region, depth, rewardLocs, xt, xtcm, colors, linecolors, posMin=posMin, posMax=posMax, binwidth=binwidth)
            pdf.savefig(fig, dpi=100)
            plt.close()