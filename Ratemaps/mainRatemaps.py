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
gsmooth = 0.25

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
    
    # load spiketimes and cluster id for good unit for a given recording
    metricsAnimaldf = pooledMetricsdf[pooledMetricsdf.aname==dname]
    spiketimes, clusterId = rmaputils.loadGoodSpikeTimes(os.path.join(ephysdirname,dname), metricsAnimaldf)
    
    # load behavior data
    behavdata = np.load(os.path.join(behavdirname,aname+'-behavop.npy'), allow_pickle=True).item()
    
    # iterate through the file to analyze data through each hallway
    for spikets in spiketimes:
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
            # calculate spike pos data
            spkPos, spkPosts, spkPosSpeed, spkPosTrial, spkPostsTrial = rmaputils.processSpikePos(spikets, trStart, 
                                                                                                  trEnd, posX, posT, 
                                                                                                  posSpeed)
            # calcualte spike map
            spkmaptr, spkmaptrsm, spkmaptrnorm, spkmap1d, spkmap1dsm = rmaputils.processSpikeMap(spkPosTrial, spkPos, 
                                                                                                 omapbins, posMin=posMin, 
                                                                                                 posMax=posMax, binwidth=binwidth, 
                                                                                                 smwindow=gsmooth)
            # calculate ratemaps
            rmaptr, rmaptrsm, rmaptrnorm, rmap1d, rmap1dsm, rmapnorm = rmaputils.processRateMap(spkmaptr, omaptr, spkmaptrsm, 
                                                                                                omaptrsm, spkmap1d, omap1d, 
                                                                                                spkmap1dsm, omap1dsm, fs=fspos)
            
            # firing statistics and properties
            spatialInfo = rmaputils.calcSpatialInformationScore(np.nanmean(rmaptrsm,0), np.nanmean(omaptrsm,0))
            sparsity, _ = rmaputils.calcSparsity(np.nanmean(rmaptrsm,0), rmaptrsm)
            stabilityInd = rmaputils.calcStabilityIndex(rmaptrsm)
            
            dadada
            
            # get placefield statistics
            placemap, pfPeakFr, pfCenter, pfSize, numFields = rmaputils.calcfieldSize(rmap1dsm, rmaptrsm, pixel2cm=bin2cm)
            # get field dispersion
            pfDispersion = rmaputils.getFieldDispersion(rmaptrsm, placemap, pfCenter, numFields, cmconversion=bin2cm, L=halllength)
            # get outinfield firing rate ratio
            outinratio = rmaputils.calcOutInFieldRatio(rmap1dsm, placemap)
        
        
#                # get phase of the individual spikes from the lfp
#                chnum = lfp_channel_num #good_clusters_ch[s]
#                spikephase = ratemaputils.loadSpikePhase(lfpfilename, chnum, spikeposts, posT[0], posT[-1], fsl=1500.)


#                ##################################################################################
#                # spatialinfo significance
#                # pval, _ = ratemaputils.calcShuffleSpatialInfo(spatialInfo, spikets, occmap1dsm, posT[-1], posT[0], trialStart, trialEnd, posX, posT, posSpeed,  occmapbins, posMin=posMin, posMax=posMax, binwidth=binwidth, fs=fs)
#                # might tackle spatial information shuffling later
#
#                # check for spatial info threshold                
#                if spatialInfo>spatialInfoTh: #and pval<0.05:
#                    # stability index
#                    stabilityInd = ratemaputils.calcStabilityIndex(ratemaptrialsm)
#                    # place field calculation
#                    placemap, placeFieldsPeakFr, placeFieldsCenter, placeFieldsSize, placeFieldsEdges, numFields = ratemaputils.calcfieldSizev2(ratemap1dsm, ratemaptrialsm, pixel2cm=binwidth)
#                    # get field dispersion
#                    placeFieldDispersion = ratemaputils.getFieldDispersion(ratemaptrialsm, placemap, placeFieldsCenter, numFields, cmconversion=posMax/binwidth, L=posMax-posMin)
#                    # get outinfield firing rate ratio
#                    outinratio, outinratioTrial = ratemaputils.calcOutInFieldRatio(ratemap1dsm, placemap, ratemaptrialsm)
#                    # calculate sparsity
##                    sparsity, spartsityTrial = ratemaputils.calcSparsity(ratemap1dsm, placemap, ratemaptrialsm)
#                else:
#                    placemap, placeFieldsPeakFr, placeFieldsCenter, placeFieldsSize, placeFieldsEdges, numFields = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
#                    outinratio = np.nan
#                    outinratioTrial = np.nan
#                    placeFieldDispersion = np.nan
#                    stabilityInd = np.nan
##                    sparsity = np.nan
##                    spartsityTrials = np.nan
#                                    
#                processed_data = {}
#                processed_data['posX'] = posX
#                processed_data['posT'] = posT
#                processed_data['posSpeed'] = posSpeed
#                processed_data['spikeposX'] = spikepos
#                processed_data['spikeposT'] = spikeposts
#                processed_data['spikeposSpd'] = spikeposspeed
#                processed_data['trialspikeposX'] = spikepostrial
#                processed_data['trialspikeposT'] = spikepoststrial
#                processed_data['rawoccmap1d'] = occmap1d
#                processed_data['smoothoccmap1d'] = occmap1dsm
#                processed_data['rawoccmaptrial'] = occmaptrial
#                processed_data['smoothoccmaptrial'] = occmaptrialsm
#                processed_data['normoccmap'] = occmaptrialnorm
#                processed_data['occmapbins'] = occmapbins
#                processed_data['rawspikemap1d'] = spikemap1d
#                processed_data['smoothspikemap1d'] = spikemap1dsm
#                processed_data['rawspikemaptrial'] = spikemaptrial
#                processed_data['smoothspikemaptrial'] = spikemaptrialsm
#                processed_data['normspikemaptrial'] = spikemaptrialnorm
#                processed_data['rawratemap1d'] = ratemap1d
#                processed_data['smoothratemap1d'] = ratemap1dsm
#                processed_data['rawratemaptrial'] = ratemaptrial
#                processed_data['smoothratemaptrial'] = ratemaptrialsm
#                processed_data['normratemaptrial'] = ratemaptrialsnorm
#                processed_data['spikephase'] = spikephase
#                processed_data['spatialInfo'] = spatialInfo
#                processed_data['stabilityIndex'] = stabilityInd
#                processed_data['placemap'] = placemap
#                processed_data['placeFieldPeakFr'] = placeFieldsPeakFr
#                processed_data['placeFieldCenter'] = placeFieldsCenter
#                processed_data['placeFieldSize'] = placeFieldsSize
#                processed_data['placeFieldEdges'] = placeFieldsEdges
#                processed_data['numFields'] = numFields
#                processed_data['placeFieldDispersion'] = placeFieldDispersion
#                processed_data['outinratio'] = outinratio
#                processed_data['outinratioTrial'] = outinratioTrial
##                processed_data['sparsity'] = sparsity
##                processed_data['sparsityTrials'] = spartsityTrial
#                processed_data['spikephase'] = spikephase
#                matfilename = os.path.join(outputdir, 'ClustId' + str(good_cluster_id[s]) +'_hall' + str(hallnum) + '_processed.mat')
#                spio.savemat(matfilename, processed_data)
#                
#                # # plotting the entire dataset ***************************************************************************************
#                fig_name = os.path.join(outputdir, 'ClustId' + str(good_cluster_id[s]) +'_hall' + str(hallnum) + '_ratemap.png')
#                fig = plt.figure(figsize=(8,8))
#                ax1 = plt.subplot2grid((10, 2), (0, 0), rowspan=2)
#                ax2 = plt.subplot2grid((10, 2), (2, 0), rowspan=2)
#                ax3 = plt.subplot2grid((10, 2), (4, 0), rowspan=2)
#                ax4 = plt.subplot2grid((10, 2), (6, 0), rowspan=2)
#                ax5 = plt.subplot2grid((10, 2), (8, 0), rowspan=2)
#                ax6 = plt.subplot2grid((10, 2), (0, 1), rowspan=2)
#                ax7 = plt.subplot2grid((10, 2), (2, 1), rowspan=2)
#                ax8 = plt.subplot2grid((10, 2), (4, 1), rowspan=2)
#                ax9 = plt.subplot2grid((10, 2), (7, 1), rowspan=4)
#                
#                # trajectory plot
#                ax1.scatter(posX,posT, s=1)
#                ax1.scatter(spikepos, spikeposts, s=0.5, c='r')
#                ax1.set_ylim([np.nanmin(posT),np.nanmax(posT)])
#                ax1.set_ylabel('Time (second)')
#                ax1.set_xlim([posMin,posMax])
#                ax1.set_xticks([])
#                
#                # occupancy map across trial
#                ax2.imshow(occmaptrial*dt, aspect='auto', origin='lower')
#                ax2.set_xticks(np.linspace(0,occmaptrial.shape[1],8)) 
#                ax2.set_xticklabels(np.linspace(posMin*bin2cm,posMax*bin2cm,8,dtype='int'))
#                ax2.set_xticks([])
#                ax2.set_ylabel('Trial Number')
#                ax2.set_title('Occupancy map trial')
#                ax2.axvline(x=rewardLocs[hallnum][0]*posMax//binwidth, ymin=0, ymax=len(good_cluster_id), c='w', linestyle='--')
#                ax2.axvline(x=rewardLocs[hallnum][1]*posMax//binwidth, ymin=0, ymax=len(good_cluster_id), c='w', linestyle='--')
#                
#                # occupancy map 1d smooth
#                ax3.plot(occmap1dsm*dt)
#                ax3.set_xticks(np.linspace(0,occmaptrial.shape[1],6)) 
#                ax3.set_xticklabels(np.linspace(posMin*bin2cm,occmaptrial.shape[1]*bin2cm,6,dtype='int'))
#                ax3.set_xticks([])
#                ax3.set_yticks([])
#                ax3.set_title('Occupancy map')
#                
#                # spike map across trials
#                ax4.imshow(spikemaptrial, aspect='auto', origin='lower')   
#                ax4.set_xticks(np.linspace(0,occmaptrial.shape[1],6)) 
#                ax4.set_xticklabels(np.linspace(posMin*bin2cm,occmaptrial.shape[1]*bin2cm,6,dtype='int'))
#                ax4.set_xticks([])
#                ax4.set_ylabel('Trial Number')
#                ax4.set_title('Spike map trial')
#                ax4.axvline(x=rewardLocs[hallnum][0]*posMax//binwidth, ymin=0, ymax=len(good_cluster_id), c='w', linestyle='--')
#                ax4.axvline(x=rewardLocs[hallnum][1]*posMax//binwidth, ymin=0, ymax=len(good_cluster_id), c='w', linestyle='--')
#                    
#                # plot spikemap1d smooth
#                ax5.plot(spikemap1dsm)  
#                ax5.set_xticks(np.linspace(0,occmaptrial.shape[1],6)) 
#                ax5.set_xticklabels(np.linspace(posMin*bin2cm,occmaptrial.shape[1]*bin2cm,6,dtype='int'), fontsize=10)
#                ax5.set_yticks([])
#                ax5.set_title('Spike map')
#                
#                # plot rate map across trials
#                ax6.imshow(ratemaptrial, aspect='auto', origin='lower')
#                ax6.set_xticks([])
#                ax6.set_ylabel('Trial Number')
#                ax6.set_title('Rate map trial')
#                ax6.axvline(x=rewardLocs[hallnum][0]*posMax//binwidth, ymin=0, ymax=len(good_cluster_id), c='w', linestyle='--')
#                ax6.axvline(x=rewardLocs[hallnum][1]*posMax//binwidth, ymin=0, ymax=len(good_cluster_id), c='w', linestyle='--')
#                
#                # plot rate map smoothed 1d
#                ax7.plot(ratemap1dsm)  
#                ax7.set_xticks(np.linspace(0,occmaptrial.shape[1],6)) 
#                ax7.set_xticklabels(np.linspace(posMin*bin2cm,occmaptrial.shape[1]*bin2cm,6,dtype='int'), fontsize=12)
#                ax7.set_yticks([])
#                ax7.set_xticks([])
#                maxfr = np.round(np.nanmax(ratemap1dsm),2)
#                meanfr = np.round(np.nanmean(ratemap1dsm),2)
#                ax7.set_title('Rate map   p:' + str(maxfr) + ' m: ' + str(meanfr))
#                
#                # plot phase precession
#                spikepos = np.array(spikepos)
#                ax8.scatter(spikepos, spikephase, 0.2, c='k')
#                ax8.scatter(spikepos, np.array(spikephase)+360, s=0.2, c='k')
#                ax7.set_xticks([])
#                ax8.set_ylabel('Phase (deg)')
#                ax8.set_xlabel('Position')
#                ax8.set_xlim([posMin,posMax])
#                
#                # firing statistics
#                ax9.text(0.1,0.65,'Mean FR =  {:.2f}'.format(meanfr), fontsize=10)
#                ax9.text(0.1,0.55,'Max FR =  {:.2f}'.format(maxfr), fontsize=10)
#                ax9.text(0.1,0.45,'# Spikes =  {:}'.format(str(len(spikeposts))), fontsize=10)
#                ax9.text(0.1,0.35,'Spatial Info =  {:.2f}'.format(spatialInfo), fontsize=10)
#                ax9.text(0.1,0.25,'Stability Index =  {:.2f}'.format(stabilityInd), fontsize=10)
#                ax9.text(0.1,0.15,'Out-In ratio =  {:.2f}'.format(outinratio), fontsize=10)
#                ax9.text(0.1,0.05,'Fields Size =  ' + str(placeFieldsSize), fontsize=10)
#                ax9.set_axis_off()
#
#                plt.suptitle('ClustId: ' + str(good_cluster_id[s]) + ' Hall: ' + str(hallnum), fontsize=20)
#                plt.tight_layout()
#                plt.savefig(fig_name, dpi=200)
#                plt.close()
#        else:
#                # save the processed data
#                processed_data = {}
#                processed_data['posX'] = np.nan
#                processed_data['posT'] = np.nan
#                processed_data['posSpeed'] = np.nan
#                processed_data['spikeposX'] = np.nan
#                processed_data['spikeposT'] = np.nan
#                processed_data['spikeposSpd'] = np.nan
#                processed_data['trialspikeposX'] = np.nan
#                processed_data['trialspikeposT'] = np.nan
#                processed_data['rawoccmap1d'] = np.nan
#                processed_data['smoothoccmap1d'] = np.nan
#                processed_data['rawoccmaptrial'] = np.nan
#                processed_data['smoothoccmaptrial'] = np.nan
#                processed_data['normoccmap'] = np.nan
#                processed_data['occmapbins'] = np.nan
#                processed_data['rawspikemap1d'] = np.nan
#                processed_data['smoothspikemap1d'] = np.nan
#                processed_data['rawspikemaptrial'] = np.nan
#                processed_data['smoothspikemaptrial'] = np.nan
#                processed_data['normspikemaptrial'] = np.nan
#                processed_data['rawratemap1d'] = np.nan
#                processed_data['smoothratemap1d'] = np.nan
#                processed_data['rawratemaptrial'] = np.nan
#                processed_data['smoothratemaptrial'] = np.nan
#                processed_data['normratemaptrial'] = np.nan
#                processed_data['spikephase'] = np.nan
#                processed_data['spatialInfo'] = np.nan
#                processed_data['stabilityIndex'] = np.nan
#                processed_data['placemap'] = np.nan
#                processed_data['placeFieldPeakFr'] = np.nan
#                processed_data['placeFieldCenter'] = np.nan
#                processed_data['placeFieldSize'] = np.nan
#                processed_data['placeFieldEdges'] = np.nan
#                processed_data['numFields'] = np.nan
#                processed_data['placeFieldDispersion'] = np.nan
#                processed_data['outinratio'] = np.nan
#                processed_data['outinratioTrial'] = np.nan
##                processed_data['sparsity'] = np.nan
##                processed_data['sparsityTrials'] = np.nan
#                matfilename =os.path.join(outputdir, 'ClustId' + str(good_cluster_id[s]) +'_hall' + str(hallnum) + '_processed.mat')
#                spio.savemat(matfilename, processed_data)
#                
#                #************ end of plotting ***********************************************************************
#                
#        # rate maps for all cells 
#        ratemap_allcell = np.array(ratemap_allcell)
#        ratemap_allcell_norm = np.array(ratemap_allcell_norm)
#        
#        # sort rows according to FR
#        cell_fr_order = np.argmax(ratemap_allcell_norm[:,:], axis=1)
#        cell_fr_order = np.argsort(cell_fr_order)
#        np.save(os.path.join(outputdir, 'cellorder_hall'+str(hallnum)+'.npy'), cell_fr_order)
#                     
#        plt.figure()
#        plt.imshow(ratemap_allcell_norm[cell_fr_order,:], cmap='jet', vmin=0, vmax=1)
#        if hallnum==1:
#            plt.axvline(x=(0.2*posMax)//binwidth, ymin=0, ymax=len(good_cluster_id), c='w')
#            plt.axvline(x=(0.65*posMax)//binwidth, ymin=0, ymax=len(good_cluster_id), c='w')
#        elif hallnum==2:
#            plt.axvline(x=(0.45*posMax)//binwidth, ymin=0, ymax=len(good_cluster_id), c='w')
#            plt.axvline(x=(0.83*posMax)//binwidth, ymin=0, ymax=len(good_cluster_id), c='w')
#        elif hallnum==28:
#            plt.axvline(x=(0.35*posMax)//binwidth, ymin=0, ymax=len(good_cluster_id), c='w')
#            plt.axvline(x=(0.83*posMax)//binwidth, ymin=0, ymax=len(good_cluster_id), c='w')
#        plt.xticks(np.linspace(posMin//binwidth,posMax//binwidth,8), np.linspace(posMin,posMax,8,dtype='int'))
#        plt.title("Hall: " + str(hallnum), fontsize=22)
#        plt.xlabel('Linear Position', fontsize=22)
#        plt.ylabel('Cell #', fontsize=22)
#        plt.colorbar()
#        plt.show()
#        
#        np.save(os.path.join(outputdir, 'hall'+str(hallnum)+'_ratemap_norm.npy'), ratemap_allcell_norm)
#        np.save(os.path.join(outputdir, 'hall'+str(hallnum)+'_ratemap.npy'), ratemap_allcell)