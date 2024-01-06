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
from matplotlib.backends.backend_pdf import PdfPages

# THRESHOLD PARAMS
speedTh = 5 # cm/s

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

# read csv with start and end time for each experimental animal
ephysdirname = r'X:\SWIL-Exp-Rajat\Spikesorted-SWIL'
behavdirname = r'X:\SWIL-Exp-Rajat\Behavior-SWIL'
rmapdirname = r'X:\SWIL-Exp-Rajat\Ratemaps-SWIL'
epochsfname = 'swil-animals.csv'
# load each recording epochs file
epochsdf = pd.read_csv(os.path.join(ephysdirname, epochsfname))
filename = epochsdf['file_name']
start_time, end_time = epochsdf['start_time'], epochsdf['end_time']
# load pooled unit metrics for all animals
pooledMetricsdf = pd.read_csv(os.path.join(ephysdirname,'analyzedMetrics',
                                           'pooledMetricsAllAnimals.csv')).drop(['Unnamed: 0'],axis=1)

# loop through all recordings from each animal
for dname, st, et in zip(filename, start_time, end_time):
    print('.........')
    print('Processing ' + str(dname))
    if 'PPC' in dname:
        aname = dname[:-3]
    else:
        aname = dname[:-2]
    animaldat = {} # save each animals data 
    
    # load spiketimes and cluster id for good unit for a given recording
    metricsAnimaldf = pooledMetricsdf[pooledMetricsdf.aname==dname]
    spiketimes, clusterId = rmaputils.loadGoodSpikeTimes(os.path.join(ephysdirname,dname), metricsAnimaldf)
    
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
        animaldat[hallwaynum] = {}
        animaldat[hallwaynum]['posX'] = posX
        animaldat[hallwaynum]['posT'] = posT
        animaldat[hallwaynum]['posSpeed'] = posSpeed
        animaldat[hallwaynum]['trStart'] = trStart
        animaldat[hallwaynum]['trEnd'] = trEnd
        animaldat[hallwaynum]['omaptr'] = omaptr
        animaldat[hallwaynum]['omaptrsm'] = omaptrsm
        animaldat[hallwaynum]['omaptrnorm'] = omaptrnorm
        animaldat[hallwaynum]['omap1d'] = omap1d
        animaldat[hallwaynum]['omap1dsm'] = omap1dsm
        animaldat[hallwaynum]['omapbins'] = omapbins
    # iterate through spike timestamps
    spikedata = {}
    for cid, spikets in zip(clusterId[10:12], spiketimes[10:12]):
        celldata = {} # dict to hold spike data for each neurons in each hallway
        for hallwaynum in HALLWAYS:
            print(cid, hallwaynum)
            celldata[hallwaynum] = {}
            #spikephase = ratemaputils.loadSpikePhase(lfpfilename, chnum, spikeposts, posT[0], posT[-1], fsl=1500.)
            
            # calculate spike pos data
            spkPos, spkPosts, spkPosSpeed, spkPosTrial, spkPostsTrial = rmaputils.processSpikePos(spikets, animaldat[hallwaynum]['trStart'], 
                                                                                                  animaldat[hallwaynum]['trEnd'],  animaldat[hallwaynum]['posX'], 
                                                                                                  animaldat[hallwaynum]['posT'], animaldat[hallwaynum]['posSpeed'])
            # calcualte spike map
            spkmaptr, spkmaptrsm, spkmaptrnorm, spkmap1d, spkmap1dsm = rmaputils.processSpikeMap(spkPosTrial, spkPos, 
                                                                                                 animaldat[hallwaynum]['omapbins'], posMin=posMin, 
                                                                                                 posMax=posMax, binwidth=binwidth, 
                                                                                                 smwindow=gsmooth)
            # calculate ratemaps
            rmaptr, rmaptrsm, rmaptrnorm, rmap1d, rmap1dsm, _ = rmaputils.processRateMap(spkmaptr, animaldat[hallwaynum]['omaptr'], spkmaptrsm, 
                                                                                                animaldat[hallwaynum]['omaptrsm'], spkmap1d, 
                                                                                                animaldat[hallwaynum]['omap1d'], spkmap1dsm, 
                                                                                                animaldat[hallwaynum]['omap1dsm'], fs=fspos)
            # firing statistics: spatial information, sparsity, stability index
            sinfo = rmaputils.calcSpatialInformationScore(rmap1dsm, omap1dsm)
            sparsity, _ = rmaputils.calcSparsity(rmap1dsm, rmaptrsm)
            stabilityInd = rmaputils.calcStabilityIndex(rmaptrsm)
            # spatialinfo significance
            sinfo_p, shuffledinfo = rmaputils.calcShuffleSpatialInfo(sinfo, spikets, omap1dsm, posT[-1], posT[0], trStart, trEnd, 
                                                       posX, posT, posSpeed,  omapbins, posMin=posMin, posMax=posMax, 
                                                       binwidth=binwidth, fs=fspos)
            # get placefield statistics: num field, field peak firing rate, size, center, dispersion
            if sinfo_p<0.01:
                pfmap, pfPeakFr, pfCenter, pfSize, pfNumFields, pfDispersion = rmaputils.calcfieldSize(rmap1dsm, rmaptrsm, pixel2cm=bin2cm, L=posMax)
            else:
                pfmap, pfPeakFr, pfCenter, pfSize, pfNumFields, pfDispersion = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
            
            # store data for each cell in each hallway
            celldata[hallwaynum]['spkPos'] = spkPos
            celldata[hallwaynum]['spkPosts'] = spkPosts
            celldata[hallwaynum]['spkPosSpeed'] = spkPosSpeed
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
            celldata[hallwaynum]['pfmap'] = pfmap
            celldata[hallwaynum]['pfPeakFr'] = pfPeakFr
            celldata[hallwaynum]['pfCenter'] = pfCenter
            celldata[hallwaynum]['pfSize'] = pfSize
            celldata[hallwaynum]['pfNumFields'] = pfNumFields
            celldata[hallwaynum]['pfDispersion'] = pfDispersion
        # plot the data across hallways for this cell
        spikedata[cid] = celldata
        del celldata
    animaldat['spikedata'] = spikedata
    del spikedata
    opfname = os.path.join(rmapdirname, dname+'-op.npy')
    np.save(opfname, animaldat)
    """
    1 x 4: + cell statistics
    5. phase vs. pos + TMI
    """
    pdfname = os.path.join(rmapdirname, dname+'-op.pdf')
    with PdfPages(pdfname) as pdf: 
        for cid in clusterId[10:12]:
            fig, ax = plt.subplots(nrows=6, ncols=4, figsize=(18,10))
            for h,hnum in enumerate(HALLWAYS):
                # trajectory plot
                ax[0][h].plot(animaldat[hnum]['posX'], animaldat[hnum]['posT'], color='k', linewidth=1, rasterized=True)
                ax[0][h].scatter(animaldat['spikedata'][cid][hnum]['spkPos'], animaldat['spikedata'][cid][hnum]['spkPosts'], color=linecolors[hnum], s=10, rasterized=True)
                # occ map
                ax[1][h].imshow(animaldat[hnum]['omaptrsm'], aspect='auto', cmap=colors[hnum], rasterized=True)
                ax[1][h].axvline(x=rewardLocs[hnum][0]*posMax//binwidth, c='k', linestyle='--')
                ax[1][h].axvline(x=rewardLocs[hnum][1]*posMax//binwidth, c='k', linestyle='--')
                # spike map
                ax[2][h].imshow(animaldat['spikedata'][cid][hnum]['spkmaptrsm'], aspect='auto', cmap=colors[hnum], rasterized=True)
                ax[2][h].axvline(x=rewardLocs[hnum][0]*posMax//binwidth, c='k', linestyle='--')
                ax[3][h].axvline(x=rewardLocs[hnum][1]*posMax//binwidth, c='k', linestyle='--')
                # rate map
                ax[3][h].imshow(animaldat['spikedata'][cid][hnum]['rmaptrsm'], aspect='auto', cmap=colors[hnum], rasterized=True)
                ax[3][h].axvline(x=rewardLocs[hnum][0]*posMax//binwidth, c='k', linestyle='--')
                ax[3][h].axvline(x=rewardLocs[hnum][1]*posMax//binwidth, c='k', linestyle='--')
                # mean over occ, spike, and ratemap plot
                ax[1][3].plot(np.nanmean(animaldat[hnum]['omaptrsm'],0), color=linecolors[hnum], rasterized=True)
                ax[2][3].plot(np.nanmean(animaldat['spikedata'][cid][hnum]['spkmaptrsm'],0), color=linecolors[hnum])
                ax[3][3].plot(np.nanmean(animaldat['spikedata'][cid][hnum]['rmaptrsm'],0), color=linecolors[hnum])
                # theta phase precession plot
                ax[4][h].axvline(x=rewardLocs[hnum][0]*posMax//binwidth, c='k', linestyle='--')
                ax[4][h].axvline(x=rewardLocs[hnum][1]*posMax//binwidth, c='k', linestyle='--')
                # label settings
                if hnum==1:
                    ax[0][h].set_ylabel('Time (s)', fontsize=18)
                    ax[1][h].set_ylabel('Trial #', fontsize=18)
                    ax[2][h].set_ylabel('Trial #', fontsize=18)
                    ax[3][h].set_ylabel('Trial #', fontsize=18)
                    ax[4][h].set_ylabel('Phase (deg)', fontsize=18)
                # axis limit and tick labels settings
                ax[1][h].set_xticks(xt), ax[1][h].set_xticklabels(xtcm)
                ax[1][3].set_xticks(xt), ax[1][3].set_xticklabels(xtcm)
                ax[2][h].set_xticks(xt), ax[2][h].set_xticklabels(xtcm)
                ax[2][3].set_xticks(xt), ax[2][3].set_xticklabels(xtcm)
                ax[3][h].set_xticks(xt), ax[3][h].set_xticklabels(xtcm)
                ax[3][3].set_xticks(xt), ax[3][3].set_xticklabels(xtcm)
                ax[4][h].set_xticks(xt), ax[4][h].set_xticklabels(xtcm)
                ax[0][h].set_xlim([0,314])
                ax[1][h].set_xlim([0,100])
                ax[2][h].set_xlim([0,100])
                ax[3][h].set_xlim([0,100])
                ax[0][h].invert_yaxis()
            ax[1][3].set_xlim([0,100])
            ax[2][3].set_xlim([0,100])
            ax[3][3].set_xlim([0,100])
            # firing statistics for each cell and hallway
            sth1 = 'pk:'+str(round(np.max(np.nanmean(animaldat['spikedata'][cid][1]['rmaptrsm'],0)),1))+', m:'+str(round(np.nanmean(np.nanmean(animaldat['spikedata'][cid][1]['rmaptrsm'],0)),1))+', i:'+str(round(animaldat['spikedata'][cid][1]['sinfo'],1))
            sth2 = 'pk:'+str(round(np.max(np.nanmean(animaldat['spikedata'][cid][2]['rmaptrsm'],0)),1))+', m:'+str(round(np.nanmean(np.nanmean(animaldat['spikedata'][cid][2]['rmaptrsm'],0)),1))+', i:'+str(round(animaldat['spikedata'][cid][2]['sinfo'],1))
            sth28 = 'pk:'+str(round(np.max(np.nanmean(animaldat['spikedata'][cid][28]['rmaptrsm'],0)),1))+', m:'+str(round(np.nanmean(np.nanmean(animaldat['spikedata'][cid][28]['rmaptrsm'],0)),1))+', i:'+str(round(animaldat['spikedata'][cid][28]['sinfo'],1))
            ax[0][3].text(0.01,0.9,'Hall#1: ' + sth1)
            ax[0][3].text(0.01,0.6,'Hall#2: ' + sth2)
            ax[0][3].text(0.01,0.3,'Hall#3: ' + sth28)
            ax[0][3].set_axis_off()
            plt.suptitle('Cluster: ' + str(cid), fontsize=20)
            plt.tight_layout()
            pdf.savefig(fig, dpi=100)
            plt.close()
    dada
#    processed_data['spikephase'] = spikephase
#    matfilename = os.path.join(outputdir, 'ClustId' + str(good_cluster_id[s]) +'_hall' + str(hallnum) + '_processed.mat')
#    spio.savemat(matfilename, processed_data)
#    
#    # # plotting the entire dataset ***************************************************************************************
#    fig_name = os.path.join(outputdir, 'ClustId' + str(good_cluster_id[s]) +'_hall' + str(hallnum) + '_ratemap.png')
#    fig = plt.figure(figsize=(8,8))
#    ax1 = plt.subplot2grid((10, 2), (0, 0), rowspan=2)
#    
#    # plot phase precession
#    spikepos = np.array(spikepos)
#    ax8.scatter(spikepos, spikephase, 0.2, c='k')
#    ax8.scatter(spikepos, np.array(spikephase)+360, s=0.2, c='k')
#    ax7.set_xticks([])
#    ax8.set_ylabel('Phase (deg)')
#    ax8.set_xlabel('Position')
#    ax8.set_xlim([posMin,posMax])
#    plt.savefig(fig_name, dpi=200)
#    plt.close()
#            
#    # rate maps for all cells 
#    ratemap_allcell = np.array(ratemap_allcell)
#    ratemap_allcell_norm = np.array(ratemap_allcell_norm)
#    
#    # sort rows according to FR
#    cell_fr_order = np.argmax(ratemap_allcell_norm[:,:], axis=1)
#    cell_fr_order = np.argsort(cell_fr_order)
#    np.save(os.path.join(outputdir, 'cellorder_hall'+str(hallnum)+'.npy'), cell_fr_order)
#                 
#    plt.figure()
#    plt.imshow(ratemap_allcell_norm[cell_fr_order,:], cmap='jet', vmin=0, vmax=1)
#    if hallnum==1:
#        plt.axvline(x=(0.2*posMax)//binwidth, ymin=0, ymax=len(good_cluster_id), c='w')
#        plt.axvline(x=(0.65*posMax)//binwidth, ymin=0, ymax=len(good_cluster_id), c='w')
#    elif hallnum==2:
#        plt.axvline(x=(0.45*posMax)//binwidth, ymin=0, ymax=len(good_cluster_id), c='w')
#        plt.axvline(x=(0.83*posMax)//binwidth, ymin=0, ymax=len(good_cluster_id), c='w')
#    elif hallnum==28:
#        plt.axvline(x=(0.35*posMax)//binwidth, ymin=0, ymax=len(good_cluster_id), c='w')
#        plt.axvline(x=(0.83*posMax)//binwidth, ymin=0, ymax=len(good_cluster_id), c='w')
#    plt.xticks(np.linspace(posMin//binwidth,posMax//binwidth,8), np.linspace(posMin,posMax,8,dtype='int'))
#    plt.title("Hall: " + str(hallnum), fontsize=22)
#    plt.xlabel('Linear Position', fontsize=22)
#    plt.ylabel('Cell #', fontsize=22)
#    plt.colorbar()
#    plt.show()