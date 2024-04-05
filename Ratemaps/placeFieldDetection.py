# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 11:04:35 2024

@author: jshobe
"""

import os
import numpy as np
import pandas as pd
import pylab as plt
import scipy.ndimage as scnd
import skimage.measure as measure

ephysdirname = r'T:\SWIL-Rajat\Spikesorted-SWIL'
rmapdirname = r'T:\SWIL-Rajat\Ratemaps-SWIL'
epochsfname = 'swil-animals.csv'
# load each recording epochs file
epochsdf = pd.read_csv(os.path.join(ephysdirname, epochsfname))
filename = epochsdf['file_name']
# load pooled unit metrics for all animals
pooledMetricsdf = pd.read_csv(os.path.join(ephysdirname,'analyzedMetrics','pooledAllAnimals.csv'))
if 'Unnamed: 0' in pooledMetricsdf.columns:
    pooledMetricsdf = pooledMetricsdf.drop(['Unnamed: 0'],axis=1)

"""
Fields with fewer than 30 spikes in total were excluded from the analyses.
"""
# loop through all recordings from each animal
for dname in filename:
    adat = np.load(os.path.join(rmapdirname, dname+'-op.npy'), allow_pickle=True).item()
    adat = adat['spikedata']
    cluid = list(adat.keys())
    for c in cluid[4:]:
        for hnum in [1,2,28]:
            hdat = adat[c][hnum]
            rmaptr = hdat['rmaptrsm']
            rmap = hdat['rmap1dsm']
            rmap = scnd.gaussian_filter1d(rmap,1.5)
            rmaptrsm = scnd.gaussian_filter1d(rmaptr,1.5)
            
            peakfallTh = 0.15
            pixel2cm = 3.14
                
            if hdat['sinfo']>0.2 and hdat['sinfo_p']<0.05:
                # detect preliminary fields
                rmap = rmap - np.nanmin(rmap)
                peakfr = np.nanmax(rmap)
                placemap = np.zeros_like(rmap, dtype=int)
                placemap[rmap >= max(0.2*peakfr, 0.5)] = 1 
                placemap = measure.label(placemap)
                placemap = measure.label(placemap, background=0)
                numfields = np.max(placemap)  
                
                # Iterate over each field
                fieldpeaks = []
                pfnumFields, pfPeak, pfCenter, pfSize, pfEdges = 0, [], [], [], []
                for field_label in range(1, numfields + 1):  # Field labels start from 1
                    # Compute the peak firing rate within the field
                    field_mask = placemap == field_label
                    fieldpk = np.nanmax(rmap[field_mask])
                    fieldpeaks.append((field_label, fieldpk))
                # Sort the fields based on their peak firing rates in descending order
                fieldpeaks.sort(key=lambda x: x[1], reverse=True)
                # Iterate over the sorted fields
                for field_label, fieldpeak in fieldpeaks:
                    fieldPixels = np.where(placemap==field_label)[0]
                    #fall off of firing rate to 20% PFR field
                    fieldthreshold = peakfallTh*fieldpeak
                    # Start extending the preliminary field
                    start_index, end_index = min(fieldPixels), max(fieldPixels)
                    
                    # extend left edge
                    peak_reached = False
                    while start_index > 0 and not peak_reached:
                        start_index -= 1
                        if rmap[start_index] >= 2*rmap[start_index+1]:
                            peak_reached = True
                        elif rmap[start_index] <= fieldthreshold:
                            break
                    # extend right edge
                    peak_reached = False
                    while end_index < len(rmap) - 1 and not peak_reached:
                        end_index += 1
                        if rmap[end_index] >= 2*rmap[end_index - 1]:
                            peak_reached = True
                        elif rmap[end_index] <= fieldthreshold:
                            break
                    # new field calculation
                    newfield = rmap[start_index:end_index + 1]
                    fieldsize = (end_index - start_index)*pixel2cm
                    fieldCenter = np.where(rmap==scnd.center_of_mass(newfield))[0]
                    
                    print(start_index, end_index+1)
                    
                    # size threshold (min. field size > 10 cm)
                    if fieldsize>10:
                        lappeakfr = np.nanmax(rmaptrsm[:,start_index:end_index+1],1)
                        fieldth = max(fieldpeak/1.25, 2.0)
                        countlaps = sum(lappeakfr>fieldth)
                        if countlaps>15 and countlaps>rmaptr.shape[0]//2:
                            pfPeak.append(fieldpeak)
                            pfCenter.append(fieldCenter)
                            pfSize.append(fieldsize)
                            pfEdges.append([start_index, end_index+1])
                            pfnumFields = pfnumFields + 1
                dadada
                pfNumFields, pfCenter, pfSize = np.nan, np.nan, np.nan
                adat[c][hnum]['pfNumFields'] = pfNumFields
                adat[c][hnum]['pfCenter'] = pfCenter
                adat[c][hnum]['pfSize'] = pfSize