#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 05:35:00 2023

@author: rajat
"""
import os
import numpy as np
import pandas as pd
from collections import OrderedDict
from scipy.ndimage.filters import gaussian_filter1d

# borrowed from Allen Institute Pipeline, Cortex-lab pipeline


def calculate_isi_violations(spike_times, spike_clusters, total_units, isi_threshold, min_isi):
    cluster_ids = np.unique(spike_clusters)
    viol_rates = np.zeros((total_units,))
    for idx, cluster_id in enumerate(cluster_ids):
        for_this_cluster = (spike_clusters == cluster_id)
        viol_rates[idx], num_violations = calcISIViolations(spike_times[for_this_cluster],
                                                       min_time = np.min(spike_times[for_this_cluster]),
                                                       max_time = np.max(spike_times[for_this_cluster]),
                                                       isi_threshold=isi_threshold,
                                                       min_isi = min_isi)
    return viol_rates


def calculate_firing_rate(spike_times, spike_clusters, total_units, epoch):
    cluster_ids = np.unique(spike_clusters)
    firing_rates = np.zeros((total_units,))
    for idx, cluster_id in enumerate(cluster_ids):
        for_this_cluster = (spike_clusters == cluster_id)
        firing_rates[idx] = calcFiringRate(spike_times[for_this_cluster],
                                        min_time = epoch[0],
                                        max_time = epoch[1])
    return firing_rates


def calculate_presence_ratio(spike_times, spike_clusters, total_units, epoch):
    cluster_ids = np.unique(spike_clusters)
    ratios = np.zeros((total_units,))
    for idx, cluster_id in enumerate(cluster_ids):
        for_this_cluster = (spike_clusters == cluster_id)
        ratios[idx] = calcPresenceRatio(spike_times[for_this_cluster],
                                                       min_time = epoch[0],
                                                       max_time = epoch[1])
    return ratios


def calculate_amplitude_cutoff(spike_clusters, amplitudes, total_units):
    cluster_ids = np.unique(spike_clusters)
    amplitude_cutoffs = np.zeros((total_units,))
    for idx, cluster_id in enumerate(cluster_ids):
        for_this_cluster = (spike_clusters == cluster_id)
        amplitude_cutoffs[idx] = calcAmpCutoff(amplitudes[for_this_cluster])
    return amplitude_cutoffs


def calcISIViolations(spike_train, min_time, max_time, isi_threshold, min_isi=0.0005):
    """Calculate ISI violations for a spike train.
    Based on metric described in Hill et al. (2011) J Neurosci 31: 8699-8705
    modified by Dan Denman from cortex-lab/sortingQuality GitHub by Nick Steinmetz
    Inputs:
    -------
    spike_train : array of spike times
    min_time : minimum time for potential spikes
    max_time : maximum time for potential spikes
    isi_threshold : threshold for isi violation
    min_isi : threshold for duplicate spikes
    Outputs:
    --------
    fpRate : rate of contaminating spikes as a fraction of overall rate
        A perfect unit has a fpRate = 0
        A unit with some contamination has a fpRate < 0.5
        A unit with lots of contamination has a fpRate > 1.0
    num_violations : total number of violations
    """

    duplicate_spikes = np.where(np.diff(spike_train) <= min_isi)[0]

    spike_train = np.delete(spike_train, duplicate_spikes + 1)
    isis = np.diff(spike_train)

    num_spikes = len(spike_train)
    num_violations = sum(isis < isi_threshold)
    violation_time = 2*num_spikes*(isi_threshold - min_isi)
    total_rate = calcFiringRate(spike_train, min_time, max_time)
    violation_rate = num_violations/violation_time
    fpRate = violation_rate/total_rate
    return fpRate, num_violations


def calcFiringRate(spike_train, min_time=None, max_time=None):
    """Calculate firing rate for a spike train.
    If no temporal bounds are specified, the first and last spike time are used.
    Inputs:
    -------
    spike_train : numpy.ndarray
        Array of spike times in seconds
    min_time : float
        Time of first possible spike (optional)
    max_time : float
        Time of last possible spike (optional)
    Outputs:
    --------
    fr : float
        Firing rate in Hz
    """
    if min_time is not None and max_time is not None:
        duration = max_time - min_time
    else:
        duration = np.max(spike_train) - np.min(spike_train)
    fr = spike_train.size / duration
    return fr


def calcPresenceRatio(spike_train, min_time, max_time, hist_win=10):
    """Calculate fraction of time the unit is present within an epoch.
    Inputs:
    -------
    spike_train : array of spike times
    min_time : minimum time for potential spikes
    max_time : maximum time for potential spikes
    hist_win: binwidth (10s)
    Outputs:
    --------
    presence_ratio : fraction of time bins in which this unit is spiking
    """
    bins = np.arange(min_time, max_time + hist_win, hist_win)
    spks_bins, _ = np.histogram(spike_train, bins)
    pr = len(np.where(spks_bins)[0]) / len(spks_bins)
    return pr

def calcAmpCutoff(amplitudes, num_histogram_bins = 500, histogram_smoothing_value = 3):
    """ Calculate approximate fraction of spikes missing from a distribution of amplitudes
    Assumes the amplitude histogram is symmetric (not valid in the presence of drift)
    Inspired by metric described in Hill et al. (2011) J Neurosci 31: 8699-8705
    Input:
    ------
    amplitudes : numpy.ndarray
        Array of amplitudes (don't need to be in physical units)
    Output:
    -------
    fraction_missing : float
        Fraction of missing spikes (0-0.5)
        If more than 50% of spikes are missing, an accurate estimate isn't possible
    """
    h,b = np.histogram(amplitudes, num_histogram_bins, density=True)
    pdf = gaussian_filter1d(h,histogram_smoothing_value)
    support = b[:-1]
    peak_index = np.argmax(pdf)
    G = np.argmin(np.abs(pdf[peak_index:] - pdf[0])) + peak_index
    bin_size = np.mean(np.diff(support))
    fraction_missing = np.sum(pdf[G:])*bin_size
    return fraction_missing


# function to calculate quality metrics
def calcQualityMetrics(dirname, fs=30000.0, params=None):
    ## Params for quality metrics
    if params is None:
        params = {}
        params['isi_threshold']=0.002
        params['min_isi']=0.0001
        params['isi_viol_th'] = 0.2
        params['presence_ratio'] = 0.75
        params['firing_rate_th'] = 0.01 
        params['amp_cutoff_th'] = 0.2 
        params['amp_th'] = 30 # uV 
    
    # load unit data
    spike_times = np.ravel(np.load(os.path.join(dirname,'spike_times.npy'), allow_pickle=True))/fs
    spike_clusters = np.ravel(np.load(os.path.join(dirname,'spike_clusters.npy'), allow_pickle=True))
    amplitudes = np.ravel(np.load(os.path.join(dirname,'amplitudes.npy'), allow_pickle=True))
    cluster_info = pd.read_csv(os.path.join(dirname,'cluster_info.tsv'), sep='\t')
    total_units = len(np.unique(spike_clusters))
    epoch = [0, spike_times[-1]]
    in_epoch = (spike_times > epoch[0]) * (spike_times < epoch[-1])
    
    # Calculate unit quality metrics
    metrics = pd.DataFrame()
    isi_viol = calculate_isi_violations(spike_times[in_epoch], spike_clusters[in_epoch], total_units, params['isi_threshold'], params['min_isi'])
    presence_ratio = calculate_presence_ratio(spike_times[in_epoch], spike_clusters[in_epoch], total_units, epoch)
    firing_rate = calculate_firing_rate(spike_times[in_epoch], spike_clusters[in_epoch], total_units, epoch)
    amplitude_cutoff = calculate_amplitude_cutoff(spike_clusters[in_epoch], amplitudes[in_epoch], total_units)
    cluster_ids = np.unique(spike_clusters)
    
    # finalize the metrics
    metrics = pd.concat((metrics, pd.DataFrame(data= OrderedDict((('cluster_id', cluster_ids),
                                    ('firing_rate' , firing_rate),
                                    ('presence_ratio' , presence_ratio),
                                    ('isi_viol' , isi_viol),
                                    ('amp_cutoff' , amplitude_cutoff),
                                    )))))
    metrics['group'] = cluster_info['group']
    metrics['depth'] = cluster_info['depth']
    metrics['ch'] = cluster_info['ch']
    metrics['num_spikes'] = cluster_info['n_spikes']
    
    # find good cell based on cutoff
    isiflag = (metrics['isi_viol']<=params['isi_viol_th']) 
    pratioflag = (metrics['presence_ratio']>=params['presence_ratio']) 
    frflag = (metrics['firing_rate']>=params['firing_rate_th']) 
    cluqualflag = (metrics['group']=='good')
    ampflag = (metrics['amp_cutoff']<=params['amp_cutoff_th']) 
    isGoodCluster = isiflag & pratioflag & frflag & ampflag & cluqualflag
    metrics['isGood'] = isGoodCluster
    print('Number of Good cluster: ' + np.sum(isGoodCluster))
    
    # save spiketimes and spikeclusters
    goodCluId = np.array(metrics['cluster_id'][metrics['isGood']])
    spiketimesGood = []
    spikeclustersGood = []
    for c, cluid in enumerate(goodCluId):
        spike_t = spike_times[np.where(spike_clusters==cluid)[0]]
        spiketimesGood.append(spike_t)
        spikeclustersGood.append(goodCluId[c])
    spiketimesGood = np.array(spiketimesGood)
    spikeclustersGood = np.array(spikeclustersGood)
    
    return metrics, spiketimesGood, spikeclustersGood
