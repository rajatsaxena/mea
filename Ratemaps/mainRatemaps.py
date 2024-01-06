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

# ---------------------------MAIN----------------------------------------------------------------------------

# THRESHOLD PARAMS
speedTh = 2.5 # cm/s

# HALLWAY information
HALLWAYS = [1,2,28] # hallways numbers, Classroom: 1, Sunset:2, Universal: 28
rewardLocs = {'1':[0.2,0.65], '2':[0.45,0.83], '28':[0.35,0.83]}

# POS variables
fspos = 30.
dt = 1./fspos
posMin = 0
posMax = 314 #cm
binwidth = 314/100 #cm
gsmooth = 0.35
bin2cm = binwidth
halllength = 314.0

# read csv with start and end time for each experimental animal
ephysdirname = r'X:\SWIL-Exp-Rajat\Spikesorted-SWIL'
behavdirname = r'X:\SWIL-Exp-Rajat\Behavior-SWIL'
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
    for cid, spikets in zip(clusterId, spiketimes):
        celldata = {} # dict to hold spike data for each neurons in each hallway
        for hallwaynum in HALLWAYS:
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
            if sinfo_p<0.05:
                pfmap, pfPeakFr, pfCenter, pfSize, pfNumFields, pfDispersion = rmaputils.calcfieldSize(rmap1dsm, rmaptrsm, pixel2cm=bin2cm, L=halllength)
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
        """
        6 x 4: each hallway + mean
        1. trajectory plot for each hallway 
        2. occupancy trial for each hallway
        3. spikemap trial for each hallway
        4. ratemap trial for each hallway
        5. 1d occ map, 1d spike map, 1d rate map
        6. phase vs. pos + cell statistics
        """
        del celldata
    animaldat['spikedata'] = spikedata
    del spikedata
    dada
#            processed_data['spikephase'] = spikephase
#            matfilename = os.path.join(outputdir, 'ClustId' + str(good_cluster_id[s]) +'_hall' + str(hallnum) + '_processed.mat')
#            spio.savemat(matfilename, processed_data)
#            
#            # # plotting the entire dataset ***************************************************************************************
#            fig_name = os.path.join(outputdir, 'ClustId' + str(good_cluster_id[s]) +'_hall' + str(hallnum) + '_ratemap.png')
#            fig = plt.figure(figsize=(8,8))
#            ax1 = plt.subplot2grid((10, 2), (0, 0), rowspan=2)
#            ax2 = plt.subplot2grid((10, 2), (2, 0), rowspan=2)
#            ax3 = plt.subplot2grid((10, 2), (4, 0), rowspan=2)
#            ax4 = plt.subplot2grid((10, 2), (6, 0), rowspan=2)
#            ax5 = plt.subplot2grid((10, 2), (8, 0), rowspan=2)
#            ax6 = plt.subplot2grid((10, 2), (0, 1), rowspan=2)
#            ax7 = plt.subplot2grid((10, 2), (2, 1), rowspan=2)
#            ax8 = plt.subplot2grid((10, 2), (4, 1), rowspan=2)
#            ax9 = plt.subplot2grid((10, 2), (7, 1), rowspan=4)
#            
#            # trajectory plot
#            ax1.scatter(posX,posT, s=1)
#            ax1.scatter(spikepos, spikeposts, s=0.5, c='r')
#            ax1.set_ylim([np.nanmin(posT),np.nanmax(posT)])
#            ax1.set_ylabel('Time (second)')
#            ax1.set_xlim([posMin,posMax])
#            ax1.set_xticks([])
#            
#            # occupancy map across trial
#            ax2.imshow(occmaptrial*dt, aspect='auto', origin='lower')
#            ax2.set_xticks(np.linspace(0,occmaptrial.shape[1],8)) 
#            ax2.set_xticklabels(np.linspace(posMin*bin2cm,posMax*bin2cm,8,dtype='int'))
#            ax2.set_xticks([])
#            ax2.set_ylabel('Trial Number')
#            ax2.set_title('Occupancy map trial')
#            ax2.axvline(x=rewardLocs[hallnum][0]*posMax//binwidth, ymin=0, ymax=len(good_cluster_id), c='w', linestyle='--')
#            ax2.axvline(x=rewardLocs[hallnum][1]*posMax//binwidth, ymin=0, ymax=len(good_cluster_id), c='w', linestyle='--')
#            
#            # occupancy map 1d smooth
#            ax3.plot(occmap1dsm*dt)
#            ax3.set_xticks(np.linspace(0,occmaptrial.shape[1],6)) 
#            ax3.set_xticklabels(np.linspace(posMin*bin2cm,occmaptrial.shape[1]*bin2cm,6,dtype='int'))
#            ax3.set_xticks([])
#            ax3.set_yticks([])
#            ax3.set_title('Occupancy map')
#            
#            # spike map across trials
#            ax4.imshow(spikemaptrial, aspect='auto', origin='lower')   
#            ax4.set_xticks(np.linspace(0,occmaptrial.shape[1],6)) 
#            ax4.set_xticklabels(np.linspace(posMin*bin2cm,occmaptrial.shape[1]*bin2cm,6,dtype='int'))
#            ax4.set_xticks([])
#            ax4.set_ylabel('Trial Number')
#            ax4.set_title('Spike map trial')
#            ax4.axvline(x=rewardLocs[hallnum][0]*posMax//binwidth, ymin=0, ymax=len(good_cluster_id), c='w', linestyle='--')
#            ax4.axvline(x=rewardLocs[hallnum][1]*posMax//binwidth, ymin=0, ymax=len(good_cluster_id), c='w', linestyle='--')
#                
#            # plot spikemap1d smooth
#            ax5.plot(spikemap1dsm)  
#            ax5.set_xticks(np.linspace(0,occmaptrial.shape[1],6)) 
#            ax5.set_xticklabels(np.linspace(posMin*bin2cm,occmaptrial.shape[1]*bin2cm,6,dtype='int'), fontsize=10)
#            ax5.set_yticks([])
#            ax5.set_title('Spike map')
#            
#            # plot rate map across trials
#            ax6.imshow(ratemaptrial, aspect='auto', origin='lower')
#            ax6.set_xticks([])
#            ax6.set_ylabel('Trial Number')
#            ax6.set_title('Rate map trial')
#            ax6.axvline(x=rewardLocs[hallnum][0]*posMax//binwidth, ymin=0, ymax=len(good_cluster_id), c='w', linestyle='--')
#            ax6.axvline(x=rewardLocs[hallnum][1]*posMax//binwidth, ymin=0, ymax=len(good_cluster_id), c='w', linestyle='--')
#            
#            # plot rate map smoothed 1d
#            ax7.plot(ratemap1dsm)  
#            ax7.set_xticks(np.linspace(0,occmaptrial.shape[1],6)) 
#            ax7.set_xticklabels(np.linspace(posMin*bin2cm,occmaptrial.shape[1]*bin2cm,6,dtype='int'), fontsize=12)
#            ax7.set_yticks([])
#            ax7.set_xticks([])
#            maxfr = np.round(np.nanmax(ratemap1dsm),2)
#            meanfr = np.round(np.nanmean(ratemap1dsm),2)
#            ax7.set_title('Rate map   p:' + str(maxfr) + ' m: ' + str(meanfr))
#            
#            # plot phase precession
#            spikepos = np.array(spikepos)
#            ax8.scatter(spikepos, spikephase, 0.2, c='k')
#            ax8.scatter(spikepos, np.array(spikephase)+360, s=0.2, c='k')
#            ax7.set_xticks([])
#            ax8.set_ylabel('Phase (deg)')
#            ax8.set_xlabel('Position')
#            ax8.set_xlim([posMin,posMax])
#            
#            # firing statistics
#            ax9.text(0.1,0.65,'Mean FR =  {:.2f}'.format(meanfr), fontsize=10)
#            ax9.text(0.1,0.55,'Max FR =  {:.2f}'.format(maxfr), fontsize=10)
#            ax9.text(0.1,0.45,'# Spikes =  {:}'.format(str(len(spikeposts))), fontsize=10)
#            ax9.text(0.1,0.35,'Spatial Info =  {:.2f}'.format(spatialInfo), fontsize=10)
#            ax9.text(0.1,0.25,'Stability Index =  {:.2f}'.format(stabilityInd), fontsize=10)
#            ax9.text(0.1,0.15,'Out-In ratio =  {:.2f}'.format(outinratio), fontsize=10)
#            ax9.text(0.1,0.05,'Fields Size =  ' + str(placeFieldsSize), fontsize=10)
#            ax9.set_axis_off()
#
#            plt.suptitle('ClustId: ' + str(good_cluster_id[s]) + ' Hall: ' + str(hallnum), fontsize=20)
#            plt.tight_layout()
#            plt.savefig(fig_name, dpi=200)
#            plt.close()
#
#            
#            #************ end of plotting ***********************************************************************
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
#    
#    np.save(os.path.join(outputdir, 'hall'+str(hallnum)+'_ratemap_norm.npy'), ratemap_allcell_norm)
#    np.save(os.path.join(outputdir, 'hall'+str(hallnum)+'_ratemap.npy'), ratemap_allcell)