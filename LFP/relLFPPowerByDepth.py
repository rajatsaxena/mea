# -*- coding: utf-8 -*-
"""
Created on Wed May  8 08:07:02 2024

@author: Justin
"""

import os, sys
import numpy as np
from tqdm import tqdm
import scipy.io as spio
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2
sys.path.insert(1, '..\..\Ratemaps')
sys.path.append('..\..\LFP')
import mea

# lfp data init
fs = 1000.0
dt = 1./fs
lfpdir = r'T:\SWIL-Rajat\LFP-SWIL'
animals = ['SWIL105','SWIL11','SWIL12','SWIL13','SWIL18','SWIL19','SWIL22','SWIL23','SWIL24','SWIL26','SWIL25']

# load default channel map
rezdat = spio.loadmat(r'T:\SWIL-Rajat\analysis-scripts\postProcessingKS2\chanMap_intan256F_bottom_SM.mat')
xcoords = np.ravel(rezdat['xcoords'])
ycoords = np.ravel(rezdat['ycoords'])
del rezdat

# sort shanks according to depth
leftshank = ycoords[xcoords<100]
leftshankidx = np.argsort(leftshank)
rightshank = ycoords[xcoords>100]
rightshankidx = np.argsort(rightshank)
# same sorting for y for left and riht shank
ysorted = ycoords[xcoords<100][leftshankidx]


# load lfp data for each shank
relativePowerAll = []
for aname in animals:
    print(aname)
    lfp = np.load(os.path.join(lfpdir, aname+'_lfp.npy'), mmap_mode='r')
    lfpvis_sh1 = lfp[:128,int(fs*20*60):int(fs*22*60)] # load 10 mins data
    lfpvis_sh2 = lfp[128:256,int(fs*20*60):int(fs*22*60)] 
    lfpppc_sh1 = lfp[256:-128,int(fs*20*60):int(fs*22*60)] # load 10 mins data
    lfpppc_sh2 = lfp[-128:,int(fs*20*60):int(fs*22*60)] 
    del lfp
    
    lfpvis_sh1  = lfpvis_sh1[leftshankidx]
    lfpvis_sh2  = lfpvis_sh2[rightshankidx]
    lfpppc_sh1 = lfpppc_sh1[leftshankidx]
    lfpppc_sh2 = lfpppc_sh2[rightshankidx]
    
    # calculate relative LFP power by depth
    relativePower = []
    for l1,l2,l3,l4 in tqdm(zip(lfpvis_sh1, lfpvis_sh2, lfpppc_sh1, lfpppc_sh2)):
        Sxx, f, t = mea.get_spectrogram(l1, fs, freq_band=None, norm=None)
        p1 = np.nansum(Sxx[(f>300)])/np.nansum(Sxx)
        Sxx, f, t = mea.get_spectrogram(l2, fs, freq_band=None, norm=None)
        p2 = np.nansum(Sxx[(f>300)])/np.nansum(Sxx)
        Sxx, f, t = mea.get_spectrogram(l3, fs, freq_band=None, norm=None)
        p3 = np.nansum(Sxx[(f>300)])/np.nansum(Sxx)
        Sxx, f, t = mea.get_spectrogram(l4, fs, freq_band=None, norm=None)
        p4 = np.nansum(Sxx[(f>300)])/np.nansum(Sxx)
        relativePower.append([p1,p2,p3,p4])
    relativePower = np.array(relativePower)
    del lfpvis_sh1, lfpvis_sh2, lfpppc_sh1, lfpppc_sh2
    
    plt.figure(figsize=(15,10))
    for i in range(4):
        plt.subplot(1,4,i+1)
        plt.plot(relativePower[:,i], ysorted, c='k', lw=2)
    plt.ylim([0, 2125])
    plt.suptitle(aname)
    plt.tight_layout()
    plt.savefig(aname+'-lfpdepth.pdf', dpi=300)
    plt.show()
        
    relativePowerAll.append(relativePower)
relativePowerAll = np.array(relativePowerAll)

import scipy.ndimage as scnd
plt.figure(figsize=(15,10))
for a in range(len(animals)):
    for i in range(4):
        plt.subplot(1,4,i+1)
        x = scnd.gaussian_filter1d(relativePowerAll[a,:,i], 3)
        x = x/np.nanmax(x)
        plt.plot(x, ysorted, lw=2, label=animals[a])
        plt.ylim([0, 2125])
plt.legend()
plt.savefig('SWIL-lfpdepth.pdf', dpi=300)
plt.show()