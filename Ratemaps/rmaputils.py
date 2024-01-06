#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 02:37:57 2020

@author: mcnlab
"""

import os
import bisect, mea
import numpy as np
import scipy.ndimage as scnd
import skimage.measure as measure
from detect_peaks import detect_peaks

#Find rightmost value less than or equal to x
def find_le(a, x):
    i = bisect.bisect_right(a, x)
    if i:
        return i-1, a[i-1]

#Find leftmost item greater than or equal to x
def find_ge(a, x):
    i = bisect.bisect_left(a, x)
    if i != len(a):
        return i, a[i]
    
# function to generate normalized 1d array
def norm1d(arr):
    return arr/np.nanmax(arr)

# function to generate 1d histogram
def hist_1d(a):
    return np.histogram(a, bins=np.arange(20,92,1))[0]

# function to run gaussian smoothing
def filter_1d(a):
    return scnd.filters.gaussian_filter(a, 1.0)
    
# binary search function
def binarySearch(data, val):
    lo, hi = 0, len(data) - 1
    best_ind = lo
    while lo <= hi:
        mid = lo + (hi - lo) // 2
        if data[mid] < val:
            lo = mid + 1
        elif data[mid] > val:
            hi = mid - 1
        else:
            best_ind = mid
            break
        # check if data[mid] is closer to val than data[best_ind] 
        if abs(data[mid] - val) < abs(data[best_ind] - val):
            best_ind = mid
    return best_ind, data[best_ind]

# load spike times and cluster id for good units
def loadGoodSpikeTimes(dirname, metricsdf):
    clusterId = np.array(metricsdf.cluster_id)
    if os.path.exists(os.path.join(dirname, 'proc-spiketimes.npy')) and os.path.exists(os.path.join(dirname, 'proc-spikeclusters.npy')):
        spikeclusters = np.load(os.path.join(dirname, 'proc-spikeclusters.npy'), allow_pickle=True)
        if (spikeclusters==clusterId).all():
            spiketimes = np.load(os.path.join(dirname, 'proc-spiketimes.npy'), allow_pickle=True)
            return spiketimes, clusterId
    return None

# function to load the trial by trial data ocupancy data
# posX - linear position
# posT - time for each position point
# posSpeed - speed for each position point
# trialStart - trial start timestamps
# trialEnd - trial end timestamps
def loadOccupancy(hallway_trial_data, posMin=0, posMax=314, binwidth=4, speedThr=4, dt=1./30):
    # variables to hold the position data
    posx = []
    post = []
    posSpd = []
    trialSt = []
    trialEt = []
    # occupancy mnap for each trials
    omaptrial = []
    # iterate thorugh each trial and interpolate the time 
    # values to ensure uniform sampling rate (fs=30 Hz)
    for x in hallway_trial_data.values():
        pos = x['pos']*posMax # position
        time = x['time'] # intan time
        if len(pos)>5:
            # interpolate so that sampling rate = 30hz
            new_time = np.arange(np.nanmin(time), np.nanmax(time)-dt, dt)
            new_pos = np.interp(new_time, time, pos)
            speed = np.ediff1d(new_pos)/np.ediff1d(new_time)
            speed = np.insert(speed,0,np.nan)
            if new_time[0]>=0:
                trialSt.append(np.round(new_time[0],2))
                trialEt.append(np.round(new_time[-1],2))
                # remove low and very high speed values
                rmv_idx = np.where((speed<=speedThr) | (speed>=120))[0]
                new_pos[rmv_idx] = np.nan
                # trial by trial occupancy map
                omap, omapbins = np.histogram(new_pos, np.arange(posMin,posMax+binwidth,binwidth))
                omaptrial.append(np.array(omap, dtype='float'))
                # note down start and end time of each trials 
                posx.extend(new_pos)
                post.extend(new_time)
                posSpd.extend(speed)
    posx = np.array(posx)
    post = np.array(post)
    posSpd = np.array(posSpd)
    trialSt = np.array(trialSt)
    trialEt = np.array(trialEt)
    omaptrial = np.array(omaptrial)
    return posx, post, posSpd, trialSt, trialEt, omaptrial

# function to process occupancy data
def processOccupancy(omaptrial, posx,  posMin=0, posMax=400, binwidth=4, smwindow=0.25):
    # smooth the occupancy map across trials
    omaptrialsm = np.apply_along_axis(filter_1d, 1, omaptrial)
    # normalized occupancy map across trials
    omaptrialnorm = []
    for omtn in omaptrial:
        omaptrialnorm.append(omtn/np.nanmax(omtn))
    omaptrialnorm = np.array(omaptrialnorm)
    # generate occupancy map
    omap1d, omapbins = np.histogram(posx, np.arange(posMin,posMax+binwidth,binwidth))
    omap1dsm = scnd.filters.gaussian_filter(omap1d, smwindow)
    # return the data
    return omaptrialsm, omaptrialnorm, omap1d, omap1dsm, omapbins

# function to process spike pos data
def processSpikePos(spikets, trialStart, trialEnd, posX, posT, posSpd):
    # variables to hold all the lap wise position and spike data
    # as well as overall spike pos, ts
    spkpos = []
    spkposts = []
    spkposspeed = []
    trspkpos = []
    trspkposts = []
    # find spike time and pos for each trial
    for tss, tse in zip(trialStart, trialEnd):
        # variables to hold trial by trial spike pos, spike ts
        tspkpos = []
        tspkposts = []
        # ensure that the spiketimestamps are after the start and end of trials
        ind = np.where((spikets>=posT[0]) & (spikets<=posT[-1]) & (spikets>=tss) & (spikets<=tse))[0]
        spkts = spikets[ind]
        for ts in spkts:
            if ts>=posT[0] and ts<=posT[-1]:
                ind, val = binarySearch(posT, ts)
                # add the timestamps to overall spike data per cell
                # and trial by trial spike data per cell
                spkpos.append(posX[ind])
                spkposts.append(val)
                spkposspeed.append(posSpd[ind])
                tspkpos.append(posX[ind])
                tspkposts.append(val)
        # add the individual trial data to overall trial by trial data
        trspkpos.append(tspkpos)
        trspkposts.append(tspkposts)
    spkpos = np.array(spkpos)
    spkposts = np.array(spkposts)
    spkposspeed = np.array(spkposspeed)
    trspkpos = np.array(trspkpos)
    trspkposts = np.array(trspkposts)
    # return spike data
    return spkpos, spkposts, spkposspeed, trspkpos, trspkposts

# fucntion process spikemap data
def processSpikeMap(spkpostrial, spkpos, omapbins, posMin=0, posMax=314, binwidth=4, smwindow=1.0):
    spkmaptrials = []
    for trspkpos in spkpostrial:
        spikecnt, _ = np.histogram(trspkpos, np.arange(posMin,posMax+binwidth,binwidth))
        spkmaptrials.append(spikecnt)
    spkmaptrials = np.array(spkmaptrials, dtype='float')
    spkmaptrials[np.isnan(spkmaptrials)] = 0.0
    # smoothed spike map across trials
    spkmaptrialsm = np.apply_along_axis(filter_1d, 1, spkmaptrials)
    spkmaptrialsm[np.isnan(spkmaptrialsm)] = 0.0
    # normalized spike map across trials
    spkmaptrialnorm = []
    for spktsm in spkmaptrials:
        spkmaptrialnorm.append(spktsm/np.nanmax(spktsm))
    spkmaptrialnorm = np.array(spkmaptrialnorm)
    # spikemap average over the session
    spkmap1d, _ = np.histogram(spkpos, bins=omapbins)
    spkmap1d = np.array(spkmap1d, dtype='float')
    # filter the spikemap average over the session
    spkmap1dsm = scnd.filters.gaussian_filter(spkmap1d, smwindow)
    return spkmaptrials, spkmaptrialsm, spkmaptrialnorm, spkmap1d, spkmap1dsm
    
# function to calculate ratemap
def processRateMap(spkmaptr, omaptr, spkmaptrsm, omaptrsm, spkmap1d, omap1d, spkmap1dsm, omap1dsm, fs=30.0):
    rmaptr = (spkmaptr/omaptr)*fs
    rmaptr[np.isinf(rmaptr)] = 0.0
    rmaptr[np.isnan(rmaptr)] = 0.0
    # smoothed ratemap across trials
    rmaptrsm = (spkmaptrsm/omaptrsm)*fs
    rmaptrsm[np.isinf(rmaptrsm)] = 0.0
    rmaptrsm[np.isnan(rmaptrsm)] = 0.0
    # normalized firing rate map across trials
    rmaptrnorm = []
    for rt in rmaptrsm:
        rmaptrnorm.append(rt/np.nanmax(rt))
    rmaptrnorm = np.array(rmaptrnorm)
    # raw ratemap for entire session
    rmap1d = (spkmap1d/omap1d)*fs
    # smooth ratemap for the entire session
    rmap1dsm = (spkmap1dsm/omap1dsm)*fs
    # normalized firing rate map
    rmapnorm = np.apply_along_axis(norm1d, axis=0, arr=rmap1dsm)                
    return rmaptr, rmaptrsm, rmaptrnorm, rmap1d, rmap1dsm, rmapnorm

# function to compute spatial information score
def calcSpatialInformationScore(rateMap, occmap):    
    si = 0
    occmap[np.isnan(occmap)] = 0.0
    if np.nanmax(rateMap)>=0:
        #firing rates in occupied pixels
        rocc = rateMap[occmap>=0]
        occ = occmap[occmap>=0]
        occ = occ/np.nansum(occ)
        #this threshold is used in Jim's program & oiginal Skaggs c program- you need to avoid pixels 
        #with rates = 0, since log(0) is undefined, but using 0.0001 as threshold instead of > 0 
        #doesn't really hurt the info score significantly
        mrate = np.mean(rocc)
        for r,p_i in zip(rocc, occ):
            if r> 0.0001:
                si = si + p_i*(r/mrate)*np.log2(r/mrate)
        if si<0:
            si=0
    else:
        si = np.nan
    return np.round(si,2)

# get sparsity calculation
def calcSparsity(rmap, rmaptrials):
    sparsity = np.nanmean(rmap)**2/ np.nanmean(rmap**2)
    sparsityTrials = []
    for rt in rmaptrials:
        sp = np.nanmean(rt)**2/ np.nanmean(rt**2)
        sparsityTrials.append(sp)
    return np.round(sparsity,3), np.round(np.array(sparsityTrials),3)
    
# get field stability
def calcStabilityIndex(rmaptrials):
    corr_ = np.corrcoef(rmaptrials)
    corr_ = np.triu(corr_)
    np.fill_diagonal(corr_,np.nan)
    corr_[np.tril_indices(corr_.shape[0], -1)] = np.nan
    stabilityInd = np.nanmean(corr_)
    return np.round(stabilityInd,4)

#function to assign phase to spike timestamps
def phase_assignment(thetapeaktime,spiketimestamps,lfpsignal,fs):
    #variable to hold spike phase
    spikephase = [] 
    #spectrogram to calculate power at different epochs for given lfp signal
#    frq,T,power = scsig.spectrogram(np.array(lfpsignal),fs=int(fs),window='hanning',nperseg=int(fs),noverlap=int(fs//2),mode='psd') 
#    # plot_only_specgram(frq, T, power)
#    #theta range used: 6-10Hz
#    frq_index_6hz = max(max(np.where(frq<=6)))
#    frq_index_10hz = min(min(np.where(frq>=10))) 
#    #calculate theta power (6-10Hz) and suprathetapower (>10Hz)
#    thetapower = sum(power[frq_index_6hz:frq_index_10hz,:],0)
#    suprathetapower = sum(power[frq_index_10hz+1:,:],0) 
#    #calculate theta supra theta power ratio, used as threshold later 
#    thetasuprathetaPowerRatio = np.divide(thetapower,suprathetapower)
#    thetasuprathetaPowerRatio = np.nan_to_num(thetasuprathetaPowerRatio)
#    #threshold settind for theta supra theta power ratio 
#    #calibrate it according to your own data
#    thetasuprathetaPowerRatiothreshold = np.nanmean(thetasuprathetaPowerRatio) - 0.8*np.nanstd(thetasuprathetaPowerRatio) 
    #iterate over each spike
    for i in range(0,len(spiketimestamps)):
        #find timestamps in spectrogram closest to ith spike
#        index, _ = binarySearch(T,spiketimestamps[i])
        #check if the theta supra theta ratio crosses the threshold (to confirm that there was actually theta modulation)
        #assign phase=nan if it does not cross the threshold
        #the second and condition is to account for that the spike occured after the 0th theta peak timestamps or before the last peak timestamps
        if spiketimestamps[i]>=thetapeaktime[0] and spiketimestamps[i]<=thetapeaktime[-1]: # and thetasuprathetaPowerRatio[index]>=thetasuprathetaPowerRatiothreshold and 
            #find before and after thetapeak timestamps nearest to the ith spike
            _, thetaPeakBefore = find_le(thetapeaktime,spiketimestamps[i])
            _, thetaPeakAfter = find_ge(thetapeaktime,spiketimestamps[i])
            #calculate the interpeak diff
            interpeak_diff = thetaPeakAfter - thetaPeakBefore
            #if theta peak after and theta peak before are same assign spike phase=0
            if thetaPeakAfter==thetaPeakBefore:
                spikephase.append(0.0)
            #this is to account for the fact that the theta peak interval is within 1/6hz to 1/12 i.e. 0.16 to 0.08
            #added a epsilon for now 0.2 to 0.05
            #everyone has to play with this as per their data
            elif 0.05 < interpeak_diff < 0.2:
                #calculate by linear interpolation
                # phase = ((spiketime - thetapeakbefore)/inter_peak_interval)*360
                phase = round((float(spiketimestamps[i]-thetaPeakBefore)/interpeak_diff)*360,3)
            #if the inter peak difference is outside assign nan
            else:
                phase = np.nan
        #phase=nan did not cross threshold
        else:
            phase = np.nan
        spikephase.append(phase)
    return spikephase

# function to spike phase data 
def loadSpikePhase(filename, chnum, spikets, sessionStart, sessionEnd, fsl=3000.):
    lfpsig = np.load(filename, mmap_mode='r')
    eeg_sig = lfpsig[chnum,:]
    del lfpsig
    dtl=1./fsl
    eeg_times = np.array(np.arange(0, len(eeg_sig)*dtl, dtl), dtype='float32')
    behav_end_idx = np.where(eeg_times<=sessionEnd)[0][-1]
    behav_start_idx = np.where(eeg_times>=sessionStart)[0][0]
    eeg_sig = eeg_sig[behav_start_idx:behav_end_idx]
    eeg_times = eeg_times[behav_start_idx:behav_end_idx]
    
    # filter in theta range 6-10hz
    theta_eeg_sig = mea.get_bandpass_filter_signal(eeg_sig,fsl,(6,10))
    #amplitude threshold for peak detection. taken directly from Sachin's code
    ampthresh = np.mean(theta_eeg_sig) + 0.7*np.std(theta_eeg_sig)
    #detect peak indices. min peak height = ampthresh, min peak distance = fs/10
    theta_peakindices = detect_peaks(theta_eeg_sig, mph=ampthresh, mpd=int(fsl/10.), show=False)
    #calculate the amplitude and time for theta peak
    thetapeaktime = eeg_times[theta_peakindices]
    spike_phase = phase_assignment(thetapeaktime,spikets,eeg_sig,fsl)
    return spike_phase

# def getautocorr(st1,st2):
#     corrbins = np.arange(-1000,1000)
#     autocorr = np.array(getBinnedISI(st1*1000, len(st1), st2*1000, len(st2)))
#     return autocorr, corrbins    

# function to get shuffled spatial information
def calcShuffleSpatialInfo(observedSI, spikets, occmap1d, endTime, startTime, trialStart, trialEnd, posX, posT, posSpeed,  occmapbins, posMin=0, posMax=314, binwidth=4, fs=30.):
    shuffledsi = []
    for i in range(1000):
        print(i)
        #get a lag from 30seconds to behavior duration - 15seconds
        lag = np.random.uniform(15, (endTime-startTime-15), 1)[0]
        #add the lag to generate timestamps
        newSpikeTimestamps = spikets + lag
        #remove spike timestamps which do not fall in maze time
        #find the indices where spiketimestamps exceed end maze time
        indices = np.where(newSpikeTimestamps>endTime)[0]
        if indices.size!=0:
            #get the modulo function to adjust for spike train exceeding end maze time
            newSpikeTimestamps[indices] = newSpikeTimestamps[indices] - endTime + startTime
            #do a cyclic shift to move the updated timestamps which earlier exceeded spike timestamps to the front 
            shuffledSpikeTimestamps = np.concatenate((newSpikeTimestamps[np.min(indices):],newSpikeTimestamps[:np.min(indices)]))
        else:
        	shuffledSpikeTimestamps = newSpikeTimestamps
        shuffledSpikeTimestamps = np.delete(shuffledSpikeTimestamps, np.where(shuffledSpikeTimestamps>endTime))
        shuffledSpikeTimestamps = np.delete(shuffledSpikeTimestamps, np.where(shuffledSpikeTimestamps<startTime))
        # function to process spike pos data
        spikepos, _, _, spikepostrial, _ = processSpikePos(shuffledSpikeTimestamps, trialStart, trialEnd, posX, posT, posSpeed)
        # generate spike map across trials
        _, _, _, _, spikemap1dsm = processSpikeMap(spikepostrial, spikepos, occmapbins, posMin=posMin, posMax=posMax, binwidth=binwidth)
        ratemap1dsm = (spikemap1dsm/occmap1d)*fs
        si = calcSpatialInformationScore(ratemap1dsm, occmap1d)
        shuffledsi.append(si)
    return 1 - np.nansum(observedSI>shuffledsi)/len(shuffledsi), np.array(shuffledsi)

# rotate place map
def rotatedPlacemap(rmap, rmaptrial):
    valleyind = detect_peaks(rmap, mph=-2, valley=True)
    if len(valleyind):
        if len(valleyind)>=2:
            valleyind = valleyind[1]
        else:
            valleyind = valleyind[0]
        rmap = np.roll(rmap,-valleyind)
        rmaptrial = np.roll(rmaptrial, -valleyind, axis=1) 
    return rmap, rmaptrial, valleyind

# field size determination
def calcfieldSize(rmap, rmaptrial, thresholdrate=0.5, peakfallTh=0.15, fieldpeakFr=0.5, fieldcutoff=8, pixel2cm=4):
    rmap = rmap - np.nanmin(rmap)
    placemap = np.zeros(len(rmap))
    placemap[np.where(rmap >= thresholdrate)] = 1 #same as mosers' 1hz threshold
    placemap = measure.label(placemap)
    placemap = measure.label(placemap, background=0)
    numfields = np.max(placemap) #find num of fields
    placeFieldsPeakFr = []
    placeFieldsCenter = []
    placeFieldsSize = []    
    ##lets try to find center of the peaks
    for i in range(1,numfields+1):
        fieldPixels = np.where(placemap==i)
        fieldpeak = max(rmap[fieldPixels])
        #fall off of firing rate to 15% PFR field
        fieldthreshold = max(peakfallTh*fieldpeak, thresholdrate/2.)
        #indixes where firing rate is less than 15%
        fieldsizecutoffind = np.where(rmap[fieldPixels]<fieldthreshold)
        #remove fields samller than 15cm and field peak firing rate less than 
        #1hz and field PFR less than 15% of the max firing rate of the neuron
        if len(fieldPixels[0])<fieldcutoff/pixel2cm or fieldpeak<fieldpeakFr:
            placemap[fieldPixels]=0
        #shorten the field size using the 15% fall off threshold
        if len(fieldsizecutoffind[0])>0:
            placemap[fieldPixels[0][fieldsizecutoffind[0]]]=0
            fieldsizecutoffind1 = np.where(rmap[fieldPixels]>=fieldthreshold)
            fieldPixels = np.array([fieldPixels[0][fieldsizecutoffind1[0]]]) 
    #retiteration after shortening things downs
    placemap = measure.label(placemap)
    placemap = measure.label(placemap, background=0)
    numFields = np.max(placemap)
    #save the field peak firing rate, peak indices, size after all the criteria 
    #this can be improved further
    for i in np.arange(1,numFields+1):
        fieldPixels = np.where(placemap==i)
        fieldpeak = max(rmap[fieldPixels])
        # check if the cell fire in more than 1/3 of the total laps
        # and is greater than 10 laps
        countlaps=0
        for rt in rmaptrial:
            lappeakfr = np.nanmax(rt[fieldPixels])
            if lappeakfr>=fieldpeak-0.5 or lappeakfr<=fieldpeak+0.5:
                countlaps+=1
        if countlaps>=rmaptrial.shape[0]//3 and countlaps>=15:      
            ind = np.where(rmap[fieldPixels]==fieldpeak)
            fieldpeakind = fieldPixels[0][ind[0]][0] 
            placeFieldsPeakFr.append(fieldpeak)
            placeFieldsCenter.append(fieldpeakind)
            placeFieldsSize.append(round(len(fieldPixels[0])*pixel2cm,2)) #converting to cm
    numFields = len(placeFieldsCenter)
    return placemap, np.array(placeFieldsPeakFr), np.array(placeFieldsCenter), np.array(placeFieldsSize), numFields   

# field size determination v2
def calcfieldSizev2(rmap, rmaptrial, pixel2cm=4.0):
    falloffTh = 0.2 # 20% of peak firing of the field
    contiguousTh = 3 # number of contiguous pixels
    pkfr = np.nanmax(rmap)
    minfieldpkfr = 0.75 # 1Hz
    shiftpixels = 3 # pixels to extend the field by
    
    # hack
    rmap = rmap - np.nanmin(rmap)
    locs = np.arange(len(rmap))
    
    ind = detect_peaks(rmap, mph=minfieldpkfr, mpd=10)
    ind = ind[np.argsort(rmap[ind])]
    ind = ind[::-1][:len(ind)]
    
    # final output variables
    fieldscenter = []
    fieldsedges = []
    fieldssize = []
    fieldspeakFr = []
    # create a zero-based place map
    placemap = np.zeros(rmap.shape)
    fieldcount=0
    
    if pkfr>=minfieldpkfr and len(ind):
        # start with max peak first
        for pi in ind:
            # get field peak rate and location
            pkloc = locs[pi]
            pk = rmap[pi]
            
            # check if the peak lies in the previously detected 
            # place fields only proceed if it does not
            if pk>=minfieldpkfr and not pkloc in np.where(placemap)[0]:
                # falloff threshold for the peak
                if pkfr>=2:
                    falloffpk = np.nanmax([falloffTh*pk, 0.75])
                else:
                    falloffpk = np.nanmin([falloffTh*pk, 0.5])
                
                # find left and right edge of the peak
                # first let's deal with right edge
                rightedgeind = np.where(rmap[pkloc:]<=falloffpk)[0]
                # if the peak lies in 80% range of the entire track length 
                # and there is a falloff value before 360
                if pkloc<=0.75*rmap.shape[0] and len(rightedgeind):
                    rightedge = pkloc+rightedgeind[0]
                else:
                    # circularly shift the map left by 100 values
                    rotatedmap = np.roll(rmap, -40)
                    rightedgeind = np.where(rotatedmap[pkloc-40:]<=falloffpk)[0]
                    if len(rightedgeind)==0:
                        rightedgeind = rmap.shape[0]
                    else:
                        rightedge = pkloc+rightedgeind[0] # right end in 360+".." values                    
                
                # now dealing with left edge
                leftedgeind = np.where(rmap[:pkloc]<=falloffpk)[0]
                # if the peak lies in 60:360 range and 
                # there is a falloff value before 0
                if pkloc>=0.4*rmap.shape[0] and len(leftedgeind):
                    leftedge = leftedgeind[-1]
                else:
                    # circularly shift the map right by 100 values
                    rotatedmap = np.roll(rmap, 30);
                    leftedgeind = np.where(rotatedmap[:pkloc+30]<=falloffpk)[0]
                    if len(leftedgeind)==0:
                        leftedge=0
                    else:
                        leftedge = leftedgeind[-1]-30; # left end in negative".." values


                # width calculation, spikepos and spikephase
                # extraction for each place field
                # check if field size is greater than continguous threshold 
                if (rightedge-leftedge)>=contiguousTh and not (rightedge in np.where(placemap)[0]) and not (leftedge in np.where(placemap)[0]):                            
                    leftedge=leftedge-shiftpixels
                    rightedge=rightedge+shiftpixels
                    
                    if rightedge>rmap.shape[0]:
#                        if not (rightedge-rmap.shape[0] in np.where(placemap)[0]):
                        fieldPixels = np.concatenate((np.arange(leftedge, rmap.shape[0]), np.arange(rightedge-rmap.shape[0])))
                    elif leftedge<0:
#                        if not (leftedge+rmap.shape[0] in np.where(placemap)[0]):
                        fieldPixels = np.concatenate((np.arange(rightedge), np.arange(rmap.shape[0]+leftedge,rmap.shape[0])))
                    else:
                        fieldPixels = np.arange(leftedge,rightedge)
                     
                    # check if the cell fire in more than 1/3 of the total laps
                    # and is greater than 10 laps
                    countlaps=0
                    for rt in rmaptrial:
                        lappeakfr = np.nanmax(rt[fieldPixels])
                        if lappeakfr>=pk-0.5 or lappeakfr<=pk+0.5:
                            countlaps+=1
                    if countlaps>=rmaptrial.shape[0]//3 and countlaps>=15:
                        if not any(x in np.where(placemap)[0] for x in fieldPixels):
                            placemap[fieldPixels] = fieldcount+1
                            fieldscenter.append(pkloc)
                            fieldspeakFr.append(pk)
                            fieldssize.append((rightedge-leftedge)*pixel2cm)
                            fieldsedges.append([leftedge,rightedge])
                            fieldcount = fieldcount+1
    numFields  = len(fieldscenter)
    return placemap, np.array(fieldspeakFr), np.array(fieldscenter), np.array(fieldssize), np.array(fieldsedges), numFields   

# calculate place field dispersion
def getFieldDispersion(rmaptrial, plmap, pfcenter, nfields, cmconversion=4, L=400):
    M = rmaptrial.shape[0]
    N = rmaptrial.shape[1]
    placeFieldDispersion = []
    for i in range(nfields):
        C = pfcenter[i]*cmconversion
        FieldDis = 0
        dispersionPF = 0
        for m in range(M):
            fieldIndex = np.where(plmap==i+1)[0]
            com = scnd.measurements.center_of_mass(rmaptrial[m,fieldIndex])
            com = com[0]+fieldIndex[0]
            C_i = com*cmconversion
            if not np.isnan(C_i):
                FieldDis += (C - C_i)**2
        dispersionPF = (L/N)*np.sqrt((1./M)*FieldDis)
        placeFieldDispersion.append(dispersionPF)
    return np.array(placeFieldDispersion)

# get out-field in-field ratio
def calcOutInFieldRatio(rmap, pmap, rmaptrials):
    inrate = np.nanmean(rmap[pmap!=0])
    outrate = np.nanmean(rmap[pmap==0])
    ratio = []
    for rt in rmaptrials:
        inrt = np.nanmean(rt[pmap!=0])
        outrt = np.nanmean(rt[pmap==0])
        ratio.append(outrt/inrt)
    return outrate/inrate, np.array(ratio)