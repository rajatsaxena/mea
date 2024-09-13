# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 03:22:23 2024

@author: saxena
"""

import numpy as np
from tqdm import tqdm
from CCGEngine import CCGEngine
from scipy.stats import poisson
from scipy.signal import convolve
from joblib import Parallel, delayed

# Helper function for smoothing
def Smooth(data, smooth):
    if smooth > 0:
        window = np.ones(smooth) / smooth
        return np.convolve(data, window, mode='same')
    return data

# CCG - Compute multiple cross- and auto-correlograms, or cross-covariances
# from Stark CCG.c mex scripts
def CCG(times, ids, binSize=0.001, duration=0.1, nBins=None, smooth=0, alpha=0.001):
    # Check parameters
    if len(times) != len(ids):
        raise ValueError("Parameters 'times' and 'id' have different lengths.")
    
    ids = np.array(ids, dtype=np.int32).reshape(-1)
    times = np.array(times, dtype=np.float64).reshape(-1)
    
    halfBins = np.ceil(duration/binSize/2.0)
    nBins = int(2*halfBins+1)
    t = np.arange(-halfBins, halfBins + 1) * binSize
    
    if len(times) <= 1:
        return np.array([]), t, np.array([]), np.array([])
    
    # Sort events in time and compute CCGs
    sort_indices = np.argsort(times)
    times = times[sort_indices]
    ids = ids[sort_indices]
    
    # Call CCGEngine (Cython function)
    counts = CCGEngine(times, ids, binSize, halfBins)
    
    # Reshape the results
    nIDs = max(ids)
    counts = counts.reshape((nIDs, nIDs, nBins)).transpose((2, 0, 1))
    
    ccg = np.flipud(counts)
    
    return counts, t


# cch 3d reshape
def CCH3D_reshape(ccg, nspks0, pruneAutoCCG=False):
    m, n, n2 = ccg.shape
    if n != n2:
        raise ValueError('Must be a 3D matrix with equal number of columns and sheets')

    # reshape cch and extract ach to prepare for deconvolution
    cch = ccg.reshape((m, n*n), order='F')  # columns are [ 11 12 ... 1n 21 22 ... 2n ... ]

    # Create index matrices
    aidx = np.tile(np.arange(n) * n + np.arange(n), (n, 1)).flatten()
    aidx1 = np.reshape(np.arange(n * n), (n, n)).T.flatten()
    aidx2 = np.arange(n * n)
    nidx = np.tile(np.arange(n).reshape(-1, 1), (1, n)).flatten()
  
    # Create ach1 and ach2
    ach1 = np.zeros((m, n * n))
    ach1[:, aidx1] = cch[:, aidx]
    ach2 = np.zeros((m, n * n))
    ach2[:, aidx2] = cch[:, aidx]
  
    # Calculate spike counts
    nspks1 = np.zeros(n * n)
    nspks2 = np.zeros(n * n)
    nspks1[aidx2] = nspks0[nidx]
    nspks2[aidx1] = nspks0[nidx]
  
    # Prune auto-CCHs
    if pruneAutoCCG:
        ridx = np.arange(0, n * n, n + 1)
        cch = np.delete(cch, ridx, axis=1)
        ach1 = np.delete(ach1, ridx, axis=1)
        ach2 = np.delete(ach2, ridx, axis=1)
        nspks1 = np.delete(nspks1, ridx)
        nspks2 = np.delete(nspks2, ridx)
  
    return cch, ach1, nspks1, ach2, nspks2

# deconvolved cch from Spivak & Stark, 2022
def cchdeconv(cch, ach1, nspks1, ach2, nspks2):
    # Argument handling
    if cch.ndim == 1:
      cch = cch[:, np.newaxis]
    if ach1.ndim == 1:
      ach1 = ach1[:, np.newaxis]
    if nspks1.ndim == 1:
      nspks1 = nspks1[np.newaxis, :]
    if ach2.ndim == 1:
      ach2 = ach2[:, np.newaxis]
    if nspks2.ndim == 1:
      nspks2 = nspks2[np.newaxis, :]
  
    m, n0 = cch.shape
    ma1, na1 = ach1.shape
    ma2, na2 = ach2.shape
  
    # Preparations
    hw = (m - 1) // 2
  
    # Scale ACHs
    ach1k = ach1 - ach1.mean(axis=0)
    ach1k /= nspks1
    hidx = np.r_[np.arange(hw), np.arange(hw + 2, m)]
    ach1k[hw, :] = 1 - ach1k[hidx, :].sum(axis=0)
  
    ach2k = ach2 - ach2.mean(axis=0)
    ach2k /= nspks2
    ach2k[hw, :] = 1 - ach2k[hidx, :].sum(axis=0)
  
    # Deconvolve ACHs from the CCH
    den = np.fft.fft(ach1k, m, axis=0) * np.fft.fft(ach2k, m, axis=0)
    dccch = np.real(np.fft.ifft(np.fft.fft(cch, m, axis=0) / den, m, axis=0))
  
    # Organize output
    dccch = np.roll(dccch, -1, axis=0)
    dccch[dccch < 0] = 0
  
    return dccch

# jitter timestamps
def jitter_timestamps(timestamps, jitter_ms=2):
    jitter_s = jitter_ms / 1000.0
    jitter_values = np.random.uniform(-jitter_s, jitter_s, size=timestamps.shape)
    jittered_timestamps = timestamps + jitter_values
    
    # Determine the range of timestamps
    min_timestamp = np.min(timestamps)
    max_timestamp = np.max(timestamps)
    timestamp_range = max_timestamp - min_timestamp

    # Circularly rotate the jittered timestamps
    jittered_timestamps = (jittered_timestamps - min_timestamp) % timestamp_range + min_timestamp
    return jittered_timestamps

# jitter timestamps
def single_jitter_run(timestamps, cluid, Nspk, jitter_ms, binwidth=0.0004, windowsize=0.05, alpha=0.001, dodeconv=False):
    ts_jitter = jitter_timestamps(timestamps, jitter_ms)
    # compute cross-correlogram 
    cch, cchbins = CCG(ts_jitter, cluid, binSize=binwidth, duration=windowsize, alpha=alpha)
    cch, ach1, nspks1, ach2, nspks2 = CCH3D_reshape(cch, Nspk)
    if dodeconv:
        # deconvolution
        cch = cchdeconv(cch, ach1, nspks1, ach2, nspks2)
    return cch

def run_jittered_ccg(timestamps, cluid, Nspk, jitter_ms=5, binwidth=0.0004, windowsize=0.05, alpha=0.001, num_jitters=1000, dodeconv=False):
    results = Parallel(n_jobs=-1)(delayed(single_jitter_run)(timestamps, cluid, Nspk, jitter_ms, binwidth, windowsize, alpha, dodeconv) for _ in tqdm(range(num_jitters)))
    return results

def cch_conv(CCH, W=5, WINTYPE='gaussian', HF=None, CALCP=True):
    CCH = np.array(CCH)
    if CCH.ndim == 1:
        CCH = CCH[:, None]
    
    nsamps, ncchs = CCH.shape
    WINTYPE = WINTYPE.lower()
    
    if HF is None:
        HF = {'gaussian': 0.6, 'median': 1}.get(WINTYPE, None)
        if HF is None:
            raise ValueError('unsupported window type')
    else:
        if HF < 0 or HF > 1:
            raise ValueError('HF not in range (0-1)')
    
    # 2. PREPARE THE CONVOLUTION WINDOW
    if WINTYPE in ['gauss', 'gaussian']:
        SDG = W/2
        if np.round(SDG)==SDG: # even W
            win = local_gausskernel(SDG, 6*SDG+1)
            cidx = int(SDG * 3 + 1)
        else:
            win = local_gausskernel(SDG, 6*SDG+2)
            cidx = int(SDG * 3 + 1.5)
    elif WINTYPE in ['median']:
        win = np.ones(W + 1 if W % 2 == 0 else W)
        cidx = W // 2 + 1
    else:
        raise ValueError('unsupported window type')
    
    win[cidx] = win[cidx] * (1 - HF)
    win = win / np.sum(win)
    
    if nsamps < (np.sqrt(2) * len(win)):
        raise ValueError('CCH-W mismatch (CCHs should be in columns; otherwise reduce W or elongate CCH)')
    
    # 3. COMPUTE A PREDICTOR BY CONVOLVING THE CCH WITH THE WINDOW
    if WINTYPE == 'median':
        pred = local_medfilt(CCH, W, HF)
    else:
        pred = local_firfilt(CCH, win)
    
    # 4. COMPUTE P-VALUE BASED ON A POISSON DISTRIBUTION WITH A CONTINUITY CORRECTION
    if CALCP:
        CCH = np.round(CCH).astype(int)
        pvals = 1 - poisson.cdf(CCH - 1, pred) - poisson.pmf(CCH, pred) * 0.5
    else:
        pvals = np.nan
    
    qvals = 1 - pvals
    
    return pvals, pred, qvals

def local_firfilt(x, W):
    C = len(W)
    D = int(np.ceil(np.ceil(C/2) - 1))
    x_padded = np.pad(x, ((D, D), (0, 0)), mode='reflect')
    Y = np.zeros((x.shape[0], x.shape[1]))
    for col in range(x.shape[1]):
        Y[:, col] = convolve(x_padded[:, col], W, mode='valid')
    return Y

def local_gausskernel(sigmaX, N):
    x = np.arange(-(N-1)/2, (N-1)/2 + 1)
    K = 1 / (2 * np.pi * sigmaX) * np.exp(-(x**2 / (2 * sigmaX**2)))
    return K

def local_medfilt(x, W, hf):
    hf = hf > 0.5
    W = (W // 2) * 2 + 1
    hw = W // 2
    x_padded = np.pad(x, ((hw, hw), (0, 0)), mode='reflect')
    y = np.zeros_like(x)
    
    for c in range(x.shape[1]):
        for r in range(x.shape[0]):
            window = x_padded[r:r+W, c]
            if hf:
                window = np.delete(window, hw)
            y[r, c] = np.median(window)
    
    return y

# calculate maximum segement with certain consecutive bins crossing threshold
def thresh_consect_slow(arr, threshold=0.001, consec_bins_min=4, consec_bins_max=15):
    max_width = 0
    max_start = 0
    max_end = 0
    
    current_start = None
    current_width = 0
    
    for i in range(len(arr)):
        if arr[i] < threshold:
            if current_start is None:
                current_start = i
            current_width += 1
        else:
            if current_width >= consec_bins_min and current_width < consec_bins_max:
                if current_width > max_width:
                    max_width = current_width
                    max_start = current_start
                    max_end = i
            current_start = None
            current_width = 0

    
    # Check if the last segment is below the threshold
    if current_width >= consec_bins_min and current_width < consec_bins_max:
        if current_width > max_width:
            max_width = current_width
            max_start = current_start
            max_end = len(arr)
    
    return max_width, (max_start, max_end)

# detect consecutive bins
def thresh_consec(arr, threshold=0.001, consec_bins_min=4, consec_bins_max=15):
    below_thresh = arr < threshold
    if np.nansum(below_thresh)== 0:
        return 0, (0, 0)  # Handle empty v case
    mask = np.convolve(below_thresh, np.ones(consec_bins_max), mode='valid') >= consec_bins_min
    if mask.any():
        start = np.argmax(mask)
        end = start + consec_bins_max
        while end < len(arr) and below_thresh[end]:
            end += 1
        return end - start, (start, end - 1)
    else:
        return 0, (0, 0)

# Calculates upper and lower bounds with Bonferroni correction.
# pred: Predicted values, alpha: Significance level, nBonf: Number of Bonferroni corrections.
def calculate_bounds(pred, alpha, nBonf):
    hiBound = poisson.ppf(1 - alpha / nBonf, pred)
    loBound = poisson.ppf(alpha / nBonf, pred)
    return hiBound, loBound

# check significance and return the strength
def check_sig(pvals, qvals, postbins, sig, sig_inh, bin_dur=0.0002, min_win_monosyn=0.005, 
              alpha=0.001, alpha2=0.01, sig_width_l=0.0008, sig_width_r=0.004):
    peak_bins = postbins[sig[postbins]]
    dip_bins = postbins[sig_inh[postbins]]
    
    # Check peak (or dip) p-value < 0.001            
    peak_pval_check = np.any(pvals[peak_bins] < alpha)
    dip_pval_check = np.any(qvals[dip_bins] < alpha)
    
    # Peak or dip width < 4 ms and alpha<0.01
    peak_width, _ = thresh_consec(pvals[peak_bins], alpha2, consec_bins_min=int(sig_width_l/bin_dur), consec_bins_max=int(sig_width_r/bin_dur))
    dip_width, _ = thresh_consec(qvals[peak_bins], alpha2, consec_bins_min=int(sig_width_l/bin_dur), consec_bins_max=int(sig_width_r/bin_dur))
    peak_width, dip_width = peak_width*bin_dur, dip_width*bin_dur
    width_check = ((peak_width >= sig_width_l) and (peak_width <= sig_width_r)) or ((dip_width >= sig_width_l) and (dip_width <= sig_width_r))
    
    # no overlap with 0-lag bin
    zero_lag_bin = len(sig)//2
    no_overlap = (zero_lag_bin not in peak_bins) and (zero_lag_bin not in dip_bins)

    return peak_pval_check, dip_pval_check, width_check, peak_width, dip_width, no_overlap


# all checks for excitatory connections
# Saxena et al., modified based on combination of multiple studies
def check_exc_conn(ccg, hiBound, loBound, pred, pvals, postbins, prebins, zero_lag_bins, 
                   bin_dur=0.0004, alpha=0.001, alpha2=0.01, sig_width=[0.0008,0.3], peaksd=3, widthsd=1):
    # check for significane for excitatory connections
    sig = ccg > hiBound
    baseline_ccg = ccg - pred
    std_baseline = np.nanstd(baseline_ccg)
    
    # 1) there should be a peak in 0.0008 ms - 0.0048 ms
    peak_post = sig[postbins]
    peak_pre = sig[prebins]

    # 2) peak p-value < alpha (0.001)
    peak_pval_check_post = pvals[postbins] < alpha
    peak_pval_check_pre = pvals[prebins] < alpha

    # 3) peak height > 2.5 sd
    peak_ht_post = np.nanmax(baseline_ccg[postbins])
    peak_ht_pre = np.nanmax(baseline_ccg[prebins])
    peakstd_check_post = peak_ht_post > peaksd*std_baseline
    peakstd_check_pre = peak_ht_pre > peaksd*std_baseline

    # 4) Peak’s width (contiguous bins with pvals < 0.01 or 0.5*peak_ht) in (0.8 - 4) ms
    sig_bins = ((pvals[postbins]<alpha2) & ((baseline_ccg[postbins]>0.5*peak_ht_post) | (baseline_ccg[postbins]>widthsd*std_baseline))).astype(int)
    contiguous_segments = np.split(sig_bins, np.where(np.diff(sig_bins) != 0)[0] + 1)
    widths_post = [len(seg)*bin_dur for seg in contiguous_segments if seg[0] == 1]
    width_check_post = any(sig_width[0] <= width <= sig_width[-1] for width in widths_post)
    
    sig_bins = ((pvals[prebins]<alpha2) & ((baseline_ccg[prebins]>0.5*peak_ht_pre) | (baseline_ccg[prebins]>widthsd*std_baseline))).astype(int)
    contiguous_segments = np.split(sig_bins, np.where(np.diff(sig_bins) != 0)[0] + 1)
    widths_pre = [len(seg)*bin_dur for seg in contiguous_segments if seg[0] == 1]
    width_check_pre = any(sig_width[0] <= width <= sig_width[-1] for width in widths_pre)
    
    # 5) Peak’s width did not overlap with zero-lag bin (indicative of common input)
    # modified to bring the range (-0.0004 - 0.0004)
    no_overlap = np.any(~sig[zero_lag_bins])
    if not no_overlap:
        zero_peak = np.nanmax(ccg[zero_lag_bins])
        if np.nanmax(ccg[postbins])>zero_peak or np.nanmax(ccg[prebins])>zero_peak:
            no_overlap = False 
    # else: 
    #     zero_peak = np.nanmax(ccg[zero_lag_bins])
    #     if zero_peak>np.nanmax(ccg[postbins]) and zero_peak>np.nanmax(ccg[prebins]):
    #         no_overlap = False
    
    # check all condition
    post_exc_check = np.any(peak_post) and np.any(peak_pval_check_post) and width_check_post and no_overlap and peakstd_check_post
    pre_exc_check = np.any(peak_pre) and np.any(peak_pval_check_pre) and width_check_pre and no_overlap and peakstd_check_pre

    return post_exc_check, pre_exc_check

# all checks for excitatory connections
# Saxena et al., modified based on combination of multiple studies
def check_inh_conn(ccg, hiBound, loBound, pred, qvals, postbins, prebins, zero_lag_bins, 
                   bin_dur=0.0004, alpha=0.001, alpha2=0.01, sig_width=[0.0008,0.3], dipsd=3, widthsd=1):
    # check for significane for excitatory connections
    sig = ccg > hiBound
    baseline_ccg = ccg - pred
    std_baseline = np.nanstd(baseline_ccg)
    sig = ccg < loBound

    # 1) there should be a dip in 0.0008 ms - 0.0048 ms
    dip_post = sig[postbins]
    dip_pre = sig[prebins]
    
    # 2) dip p-value < alpha (0.001)
    dip_pval_check_post = qvals[postbins] < alpha
    dip_pval_check_pre = qvals[prebins] < alpha
    
    # 3) dip amp < 1.5 sd
    dip_ht_post = np.nanmin(baseline_ccg[postbins])
    dip_ht_pre = np.nanmin(baseline_ccg[prebins])
    dipstd_check_post = dip_ht_post < dipsd*std_baseline
    dipstd_check_pre = dip_ht_pre < dipsd*std_baseline

    # 4) Dip’s width (contiguous bins with pvals < 0.01 or 0.5*peak_ht) in (0.8 - 4) ms
    sig_bins = ((qvals[postbins]<alpha2) & ((baseline_ccg[postbins]<0.5*dip_ht_post) | (baseline_ccg[postbins]<-widthsd*std_baseline))).astype(int)
    contiguous_segments = np.split(sig_bins, np.where(np.diff(sig_bins) != 0)[0] + 1)
    widths_post = [len(seg)*bin_dur for seg in contiguous_segments if seg[0] == 1]
    width_check_post = any(sig_width[0] <= width <= sig_width[-1] for width in widths_post)
    
    sig_bins = ((qvals[prebins]<alpha2) & ((baseline_ccg[prebins]<0.5*dip_ht_pre) | (baseline_ccg[prebins]<-widthsd*std_baseline))).astype(int)
    contiguous_segments = np.split(sig_bins, np.where(np.diff(sig_bins) != 0)[0] + 1)
    widths_pre = [len(seg)*bin_dur for seg in contiguous_segments if seg[0] == 1]
    width_check_pre = any(sig_width[0] <= width <= sig_width[-1] for width in widths_pre)
    
    # 5) Peak’s width did not overlap with zero-lag bin (indicative of common input)
    # modified to bring the range (-0.0004 - 0.0004)
    no_overlap = np.any(~sig[zero_lag_bins])
    if not no_overlap:
        zero_dip = np.nanmin(ccg[zero_lag_bins])
        if (np.nanmin(ccg[postbins]) < zero_dip) or (np.nanmin(ccg[prebins]) > zero_dip):
            no_overlap = False 
    else:
        zero_dip = np.nanmin(ccg[zero_lag_bins])
        if zero_dip<np.nanmin(ccg[postbins]) and zero_dip<np.nanmin(ccg[prebins]):
            no_overlap = False
    
    # check all condition
    post_inh_check = np.any(dip_post) and np.any(dip_pval_check_post) and no_overlap and dipstd_check_post and width_check_post
    pre_inh_check = np.any(dip_pre) and np.any(dip_pval_check_pre) and no_overlap and dipstd_check_pre and width_check_pre

    return post_inh_check, pre_inh_check



# compute significance with continuity correction English et al., 2017 
# with checks on max amp, width, common inputs, etc.
def ccg_sig_continuity(cch, t, nCells, WINTYPE='median', bin_dur=0.0002, min_win_monosyn=0.005, 
                       alpha=0.001, alpha2=0.01, sig_width_l=0.0008, sig_width_r=0.004):
    W = int(0.01/bin_dur) # 10 ms window for convolution

    # Initialize arrays for ouput
    Pval = np.full([len(t), nCells, nCells], np.nan)
    Pred = np.zeros([len(t), nCells, nCells])
    Bounds = np.zeros([cch.shape[0], nCells, nCells, 2])
    Pred[:] = np.nan
    Bounds[:] = np.nan
    sig_con = []
    sig_con_inh = []
    Pcausal = np.full((nCells, nCells), np.nan)
    Pcausal_inh = np.full((nCells, nCells), np.nan)
    syn_strength = np.full((nCells, nCells), np.nan)
    syn_ratio = np.full((nCells, nCells), np.nan)
    syn_strength_inh = np.full((nCells, nCells), np.nan)
    syn_ratio_inh = np.full((nCells, nCells), np.nan)

    # go through each pair 
    for refcid in range(nCells):
        for targetcid in range(refcid+1, nCells):
            # Median filter convolution based prediction
            ccg = cch[:, refcid, targetcid]

            # calculate null-distribution using convolution based on Stark & Abeles
            pvals, pred, qvals = cch_conv(ccg, W, WINTYPE=WINTYPE)
            pvals, pred, qvals = np.ravel(pvals), np.ravel(pred), np.ravel(qvals)
            
            # Store predicted values and pvalues for subsequent plotting
            Pred[:, refcid, targetcid] = pred
            Pval[:, refcid, targetcid] = pvals
            Pred[:, targetcid, refcid] = np.flipud(pred)
            Pval[:, targetcid, refcid] = np.flipud(pvals)
        
            # Calculate upper and lower limits with Bonferroni correction
            hiBound, loBound = calculate_bounds(pred, alpha, np.ceil(min_win_monosyn/bin_dur)*2)
            Bounds[:, refcid, targetcid, 0] = hiBound
            Bounds[:, refcid, targetcid, 1] = loBound
            Bounds[:, targetcid, refcid, 0] = np.flipud(hiBound)
            Bounds[:, targetcid, refcid, 1] = np.flipud(loBound)
        
            ############## EXCITATORY connections ##############
            # Poisson continuity correction
            # looking for whether postbins peak > prebins and reverse
            prebins = np.arange(np.round(len(ccg)/2 - 0.004/bin_dur) - 1, np.ceil(len(ccg)/2)).astype(int)
            postbins = np.arange(np.round(len(ccg)/2 + 0.0008/bin_dur) - 1, np.round(len(ccg)/2 + 0.005/bin_dur)).astype(int)
            sig = ccg > hiBound
            cchud = np.flipud(ccg)
            sigud = np.flipud(sig)
            sigpost = np.nanmax(ccg[postbins]) > poisson.ppf(1 - alpha, np.nanmax(ccg[prebins]))
            sigpre = np.nanmax(cchud[postbins]) > poisson.ppf(1 - alpha, np.nanmax(cchud[prebins]))
            
            if np.any(sigud[postbins]) and sigpre:
                sig_con.append([targetcid, refcid])
            if np.any(sig[postbins]) and sigpost:
                sig_con.append([refcid, targetcid])

            # check for significance for peaks
            pvals_causal = 1 - poisson.cdf(np.nanmax(ccg[postbins]) - 1, np.nanmax(ccg[prebins])) - 0.5*poisson.pmf(np.nanmax(ccg[postbins]), np.nanmax(ccg[prebins]))
            pvals_causalud = 1 - poisson.cdf(np.nanmax(cchud[postbins]) - 1, np.nanmax(cchud[prebins])) - 0.5*poisson.pmf(np.nanmax(cchud[postbins]), np.nanmax(cchud[prebins]))
            # can go negative for very small p-val - beyond comp. sig. dig
            if pvals_causalud < 0:
                pvals_causalud = 0
            if pvals_causal < 0:
                pvals_causal = 0
            Pcausal[refcid, targetcid] = pvals_causal
            Pcausal[targetcid, refcid] = pvals_causalud
        
            ########### INHIBITORY connections #########################
            sig_inh = ccg < loBound
            sigud_inh = np.flipud(sig_inh)
            sigpost_inh = np.nanmin(ccg[postbins]) < poisson.ppf(alpha, np.nanmin(ccg[prebins]))
            sigpre_inh = np.nanmin(cchud[postbins]) < poisson.ppf(alpha, np.nanmin(cchud[prebins]))
            
            if np.any(sigud_inh[postbins]) and sigpre_inh:
                sig_con_inh.append([targetcid, refcid])
            if np.any(sig_inh[postbins]) and sigpost_inh:
                sig_con_inh.append([refcid, targetcid])

            # check for significance in inhibitory dips
            pvals_causal = poisson.cdf(np.nanmin(ccg[postbins]), np.nanmin(ccg[prebins])) + 0.5*poisson.pmf(np.nanmin(ccg[postbins]), np.nanmin(ccg[prebins]))
            pvals_causalud = poisson.cdf(np.nanmin(cchud[postbins]), np.nanmin(cchud[prebins])) + 0.5*poisson.pmf(np.nanmin(cchud[postbins]), np.nanmin(cchud[prebins]))
            # can go negative for very small p-val - beyond comp. sig. dig
            if pvals_causalud < 0:
                pvals_causalud = 0
            if pvals_causal < 0:
                pvals_causal = 0
            Pcausal_inh[refcid, targetcid] = pvals_causal
            Pcausal_inh[targetcid, refcid] = pvals_causalud

            ## Additional check, for both sides refid -> targetid, targetid -> refid
            if sigpost or sigpre or sigpost_inh or sigpre_inh:
                # check significance in ref -> target
                peak_pval_check, dip_pval_check, width_check, peak_width, dip_width, no_overlap = check_sig( pvals, qvals, postbins, sig, sig_inh, bin_dur, 
                                                                                                            min_win_monosyn, alpha, alpha2, sig_width_l, sig_width_r)
                # Spike transmission probability and max peak/dip ratio
                excess_transmission, max_ratio = np.nan, np.nan
                excess_inhibition, min_ratio = np.nan, np.nan
                if peak_pval_check and width_check and no_overlap:
                    excess_transmission = np.nansum((ccg[postbins] - pred[postbins]) / np.nansum(ccg[postbins]))
                    max_ratio = np.nanmax(ccg[postbins]) / np.nanmean(pred[postbins])
                    syn_strength[refcid, targetcid] = excess_transmission
                    syn_ratio[refcid, targetcid] = max_ratio
                    
                if dip_pval_check and width_check and no_overlap:
                    excess_inhibition = np.nansum((ccg[postbins] - pred[postbins]) / np.nansum(ccg[postbins]))
                    min_ratio = np.nanmin(ccg[postbins]) / np.nanmean(pred[postbins])
                    syn_strength_inh[refcid, targetcid] = excess_inhibition
                    syn_ratio_inh[refcid, targetcid] = min_ratio
                

                # check significance reverse direction targetid -> refid
                pvalsud, predud, qvalsud = np.flipud(pvals), np.flip(pred), np.flipud(qvals)
                peak_pval_check, dip_pval_check, width_check, peak_width, dip_width, no_overlap = check_sig(pvalsud, qvalsud, postbins, sigud, sigud_inh, bin_dur, 
                                                                                                            min_win_monosyn, alpha, alpha2, sig_width_l, sig_width_r)
                # Spike transmission probability and max peak/dip ratio
                if peak_pval_check and width_check and no_overlap:
                    excess_transmission = np.nansum((cchud[postbins] - predud[postbins]) / np.nansum(cchud[postbins]))
                    max_ratio = np.nanmax(cchud[postbins]) / np.nanmean(predud[postbins])
                    syn_strength[targetcid, refcid] = excess_transmission
                    syn_ratio[targetcid, refcid] = max_ratio
                
                if dip_pval_check and width_check and no_overlap:
                    excess_inhibition = np.nansum((cchud[postbins] - predud[postbins]) / np.nansum(cchud[postbins]))
                    min_ratio = np.nanmin(cchud[postbins]) / np.nanmean(predud[postbins])
                    syn_strength_inh[targetcid, refcid] = excess_inhibition
                    syn_ratio_inh[targetcid, refcid] = min_ratio
    
    return Pval, Pred, Bounds, sig_con, sig_con_inh, Pcausal, Pcausal_inh, syn_strength, syn_ratio, syn_strength_inh, syn_ratio_inh
