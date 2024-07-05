# -*- coding: utf-8 -*-
"""
Created on Wed May  8 08:07:02 2024

@author: rajat
"""

import os, mea, sys
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2
sys.path.insert(1, '..\..\Ratemaps')
sys.path.append('..\..\LFP')


# calculate CSD
def calcCSD(lfp, w):
    lfp = -lfp
    CSD = np.zeros_like(lfp)
    for i in range(lfp.shape[0]):
        if i - 2 >= 0 and i + 2 < lfp.shape[0]:
            u1 = lfp[i - 2, :]
            u2 = lfp[i - 1, :]
            u3 = lfp[i, :]
            u4 = lfp[i + 1, :]
            u5 = lfp[i + 2, :]  
            CSD[i, :] = -(w[0] * u1 + w[1] * u2 + w[2] * u3 + w[3] * u4 + w[4] * u5) / (2 * dz * 2 * dz)
    # CSD = CSD[2:-2, :]  # Adjusted to avoid boundary effects
    return CSD


# calculate bandpass filter
def bpfilter(lfp, fs, frange):
    lfpfilt = []
    for l in lfp:
        lfpfilt.append(mea.butter_bandpass_filter(l, frange[0], frange[1], fs, order=2))
    return np.array(lfpfilt)

# plot CSD
def plotCSD(timestamps, CSD, lfp_frag, ax, tstr=None):
    lfp_frag = -lfp_frag.T
    cmax = np.max(np.abs(CSD))
    # plot CSD
    ax.contourf(np.arange(CSD.shape[1]), np.arange(CSD.shape[0]), CSD, 40, cmap='jet', vmin=-cmax, vmax=cmax)
    ax.set_xlabel('time (s)')
    ax.set_title(tstr)
    ax.invert_yaxis()
    # Plot LFP
    for ch in range(lfp_frag.shape[1]):
        ax.plot(np.arange(len(lfp_frag)), lfp_frag[:,ch]/0.3e3+ch-1, 'k', linewidth=1)
    
    plt.tight_layout()
    plt.show()

        
# lfp data init
fs = 1000.0
dt = 1./fs
dz = 50
lfpdir = r'T:\SWIL-Rajat\LFP-SWIL'
animals = ['SWIL22']
w = np.array([0.23, 0.08, -0.62, 0.08, 0.23])
frange = [0.1, 4]


# load default channel map
rezdat = spio.loadmat(r'T:\SWIL-Rajat\analysis-scripts\postProcessingKS2\chanMap_intan256F_bottom_SM.mat')
xcoords = np.ravel(rezdat['xcoords'])
ycoords = np.ravel(rezdat['ycoords'])
del rezdat

# sort shanks according to depth
ysh1 = ycoords[xcoords<100]
xsh1 = xcoords[xcoords<100]
leftshankidx = np.argsort(ysh1)
xsh1center = xsh1[leftshankidx]
ysh2 = ycoords[xcoords>100]
xsh2 = xcoords[xcoords>100]
rightshankidx = np.argsort(ysh2)
xsh2center = xsh2[rightshankidx]
# same sorting for y for left and riht shank
ysorted = ycoords[xcoords<100][leftshankidx]
ysh1center = ysorted[xsh1center==0]
ysh2center = ysorted[xsh2center==500]


# load lfp data for each shank
relativePowerAll = []
for aname in animals:
    lfp = np.load(os.path.join(lfpdir, aname+'_lfp.npy'), mmap_mode='r')
    lfpts = np.linspace(0, lfp.shape[1]*1./fs, lfp.shape[1])
    idx = (lfpts>6528) & (lfpts<6530)    
    lfpvis_sh1 = lfp[:128,idx]
    lfpvis_sh2 = lfp[128:256,idx] 
    lfpppc_sh1 = lfp[256:-128,idx]
    lfpppc_sh2 = lfp[-128:,idx]
    del lfp
    
    # subset for each vertical strip
    lfpvis_sh1  = lfpvis_sh1[leftshankidx][xsh1center==0]
    lfpvis_sh2  = lfpvis_sh2[rightshankidx][xsh2center==500]
    lfpppc_sh1 = lfpppc_sh1[leftshankidx][xsh1center==0]
    lfpppc_sh2 = lfpppc_sh2[rightshankidx][xsh2center==500]
    
    # bandpass filter in UP_DOWN states
    lfpvis_sh1 = bpfilter(lfpvis_sh1, fs, frange)
    lfpvis_sh2 = bpfilter(lfpvis_sh2, fs, frange)
    lfpppc_sh1 = bpfilter(lfpppc_sh1, fs, frange)
    lfpppc_sh2 = bpfilter(lfpppc_sh2, fs, frange)
    
    # run CSD analysis for only cortical electrodes
    CSDvis_sh1 = calcCSD(lfpvis_sh1[:30,:], w)
    CSDvis_sh2 = calcCSD(lfpvis_sh2[:30,:], w)    
    CSDppc_sh1 = calcCSD(lfpppc_sh1[:30,:], w)
    CSDppc_sh2 = calcCSD(lfpppc_sh2[:30,:], w)
    
    # plt.figure()
    # for l, lfp in enumerate(lfpvis_sh1):
    #     plt.plot(lfp+l*100)
    # plt.show()
    
    # plt.figure()
    # for l, lfp in enumerate(lfpppc_sh1):
    #     plt.plot(lfp+l*100)
    # plt.show()
    
    # plot the CSDs
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    plotCSD(lfpts[idx], CSDvis_sh1, lfpvis_sh1[:30,:], axs[0])
    plotCSD(lfpts[idx], CSDppc_sh1, lfpppc_sh1[:30,:], axs[1])
    
    fsda