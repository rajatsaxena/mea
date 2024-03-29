# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 16:38:51 2024

@author: jshobe
"""

#! /bin/env python
#
# Michael Gibson 17 July 2015
# Modified Adrian Foy Sep 2018

import glob
import numpy as np
import os, time
import cupy as cp # added by Rajat
import scipy.signal as spsig
from natsort import natsorted
from intanutil.read_header import read_header
from intanutil.get_bytes_per_data_block import get_bytes_per_data_block
from intanutil.read_one_data_block import read_one_data_block
from intanutil.notch_filter import notch_filter
from intanutil.data_to_result import data_to_result
#sys.stdout = open('intanRunLog.txt', 'w')

def read_data(filename):
    """Reads Intan Technologies RHD2000 data file generated by evaluation board GUI.
    
    Data are returned in a dictionary, for future extensibility.
    """

    tic = time.time()
    fid = open(filename, 'rb')
    filesize = os.path.getsize(filename)

    header = read_header(fid)

    print('Found {} amplifier channel{}.'.format(header['num_amplifier_channels'], plural(header['num_amplifier_channels'])))
    print('Found {} auxiliary input channel{}.'.format(header['num_aux_input_channels'], plural(header['num_aux_input_channels'])))
    print('Found {} supply voltage channel{}.'.format(header['num_supply_voltage_channels'], plural(header['num_supply_voltage_channels'])))
    print('Found {} board ADC channel{}.'.format(header['num_board_adc_channels'], plural(header['num_board_adc_channels'])))
    print('Found {} board digital input channel{}.'.format(header['num_board_dig_in_channels'], plural(header['num_board_dig_in_channels'])))
    print('Found {} board digital output channel{}.'.format(header['num_board_dig_out_channels'], plural(header['num_board_dig_out_channels'])))
    print('Found {} temperature sensors channel{}.'.format(header['num_temp_sensor_channels'], plural(header['num_temp_sensor_channels'])))
    print('')

    # Determine how many samples the data file contains.
    bytes_per_block = get_bytes_per_data_block(header)

    # How many data blocks remain in this file?
    data_present = False
    bytes_remaining = filesize - fid.tell()
    if bytes_remaining > 0:
        data_present = True

    if bytes_remaining % bytes_per_block != 0:
        raise Exception('Something is wrong with file size : should have a whole number of data blocks')

    num_data_blocks = int(bytes_remaining / bytes_per_block)

    num_amplifier_samples = header['num_samples_per_data_block'] * num_data_blocks
    num_aux_input_samples = int((header['num_samples_per_data_block'] / 4) * num_data_blocks)
    num_supply_voltage_samples = 1 * num_data_blocks
    num_board_adc_samples = header['num_samples_per_data_block'] * num_data_blocks
    num_board_dig_in_samples = header['num_samples_per_data_block'] * num_data_blocks
    num_board_dig_out_samples = header['num_samples_per_data_block'] * num_data_blocks

    record_time = num_amplifier_samples / header['sample_rate']

    if data_present:
        print('File contains {:0.3f} seconds of data.  Amplifiers were sampled at {:0.2f} kS/s.'.format(record_time, header['sample_rate'] / 1000))
    else:
        print('Header file contains no data.  Amplifiers were sampled at {:0.2f} kS/s.'.format(header['sample_rate'] / 1000))

    if data_present:
        # Pre-allocate memory for data.
        print('')
        print('Allocating memory for data...')

        data = {}
        if (header['version']['major'] == 1 and header['version']['minor'] >= 2) or (header['version']['major'] > 1):
            data['t_amplifier'] = np.zeros(num_amplifier_samples, dtype=np.int_)
        else:
            data['t_amplifier'] = np.zeros(num_amplifier_samples, dtype=np.uint)

        data['amplifier_data'] = np.zeros([header['num_amplifier_channels'], num_amplifier_samples], dtype=np.uint)
        data['aux_input_data'] = np.zeros([header['num_aux_input_channels'], num_aux_input_samples], dtype=np.uint)
        data['supply_voltage_data'] = np.zeros([header['num_supply_voltage_channels'], num_supply_voltage_samples], dtype=np.uint)
        data['temp_sensor_data'] = np.zeros([header['num_temp_sensor_channels'], num_supply_voltage_samples], dtype=np.uint)
        data['board_adc_data'] = np.zeros([header['num_board_adc_channels'], num_board_adc_samples], dtype=np.uint)
        
        # by default, this script interprets digital events (digital inputs and outputs) as booleans
        # if unsigned int values are preferred(0 for False, 1 for True), replace the 'dtype=np.bool' argument with 'dtype=np.uint' as shown
        # the commented line below illustrates this for digital input data; the same can be done for digital out
        
#        data['board_dig_in_data'] = np.zeros([header['num_board_dig_in_channels'], num_board_dig_in_samples], dtype=np.uint)
        data['board_dig_in_data'] = np.zeros([header['num_board_dig_in_channels'], num_board_dig_in_samples], dtype=np.bool_)
        data['board_dig_in_raw'] = np.zeros(num_board_dig_in_samples, dtype=np.uint)
        
        data['board_dig_out_data'] = np.zeros([header['num_board_dig_out_channels'], num_board_dig_out_samples], dtype=np.bool_)
        data['board_dig_out_raw'] = np.zeros(num_board_dig_out_samples, dtype=np.uint)

        # Read sampled data from file.
        print('Reading data from file...')

        # Initialize indices used in looping
        indices = {}
        indices['amplifier'] = 0
        indices['aux_input'] = 0
        indices['supply_voltage'] = 0
        indices['board_adc'] = 0
        indices['board_dig_in'] = 0
        indices['board_dig_out'] = 0

        print_increment = 10
        percent_done = print_increment
        for i in range(num_data_blocks):
            read_one_data_block(data, header, indices, fid)

            # Increment indices
            indices['amplifier'] += header['num_samples_per_data_block']
            indices['aux_input'] += int(header['num_samples_per_data_block'] / 4)
            indices['supply_voltage'] += 1
            indices['board_adc'] += header['num_samples_per_data_block']
            indices['board_dig_in'] += header['num_samples_per_data_block']
            indices['board_dig_out'] += header['num_samples_per_data_block']            

            fraction_done = 100 * (1.0 * i / num_data_blocks)
            if fraction_done >= percent_done:
                print('{}% done...'.format(percent_done))
                percent_done = percent_done + print_increment

        # Make sure we have read exactly the right amount of data.
        bytes_remaining = filesize - fid.tell()
        if bytes_remaining != 0: raise Exception('Error: End of file not reached.')



    # Close data file.
    fid.close()

    if (data_present):
        print('Parsing data...')

        # Extract digital input channels to separate variables.
        for i in range(header['num_board_dig_in_channels']):
            data['board_dig_in_data'][i, :] = np.not_equal(np.bitwise_and(data['board_dig_in_raw'], (1 << header['board_dig_in_channels'][i]['native_order'])), 0)

        # Extract digital output channels to separate variables.
        for i in range(header['num_board_dig_out_channels']):
            data['board_dig_out_data'][i, :] = np.not_equal(np.bitwise_and(data['board_dig_out_raw'], (1 << header['board_dig_out_channels'][i]['native_order'])), 0)

        # Scale voltage levels appropriately.
#        data['amplifier_data'] = np.multiply(0.195, (data['amplifier_data'].astype(np.int32) - 32768))      # units = microvolts
        data['aux_input_data'] = np.multiply(37.4e-6, data['aux_input_data'])               # units = volts
        data['supply_voltage_data'] = np.multiply(74.8e-6, data['supply_voltage_data'])     # units = volts
        if header['eval_board_mode'] == 1:
            data['board_adc_data'] = np.multiply(152.59e-6, (data['board_adc_data'].astype(np.int32) - 32768)) # units = volts
        elif header['eval_board_mode'] == 13:
            data['board_adc_data'] = np.multiply(312.5e-6, (data['board_adc_data'].astype(np.int32) - 32768)) # units = volts
        else:
            data['board_adc_data'] = np.multiply(50.354e-6, data['board_adc_data'])           # units = volts
        data['temp_sensor_data'] = np.multiply(0.01, data['temp_sensor_data'])               # units = deg C

        # Check for gaps in timestamps.
        num_gaps = np.sum(np.not_equal(data['t_amplifier'][1:]-data['t_amplifier'][:-1], 1))
        if num_gaps == 0:
            print('No missing timestamps in data.')
        else:
            print(filename, num_gaps)
            sfname = filename[:-3]+'npy'
            np.save(sfname,num_gaps)
            print('Warning: {0} gaps in timestamp data found.  Time scale will not be uniform!'.format(num_gaps))

        # Scale time steps (units = seconds).
        data['t_amplifier'] = data['t_amplifier'] / header['sample_rate']
        data['t_aux_input'] = data['t_amplifier'][range(0, len(data['t_amplifier']), 4)]
        data['t_supply_voltage'] = data['t_amplifier'][range(0, len(data['t_amplifier']), header['num_samples_per_data_block'])]
        data['t_board_adc'] = data['t_amplifier']
        data['t_dig'] = data['t_amplifier']
        data['t_temp_sensor'] = data['t_supply_voltage']

        # If the software notch filter was selected during the recording, apply the
        # same notch filter to amplifier data here.
        if False: #header['notch_filter_frequency'] > 0:
            print('Applying notch filter...')

            print_increment = 10
            percent_done = print_increment
            for i in range(header['num_amplifier_channels']):
                data['amplifier_data'][i,:] = notch_filter(data['amplifier_data'][i,:], header['sample_rate'], header['notch_filter_frequency'], 10)

                fraction_done = 100 * (i / header['num_amplifier_channels'])
                if fraction_done >= percent_done:
                    print('{}% done...'.format(percent_done))
                    percent_done += print_increment
    else:
        data = [];

    # Move variables to result struct.
    result = data_to_result(header, data, data_present)

    print('Done!  Elapsed time: {0:0.1f} seconds'.format(time.time() - tic))
    return result['t_amplifier'], result['amplifier_data'], data['board_dig_in_data'], data['board_adc_data'], result['frequency_parameters']['amplifier_sample_rate']


def plural(n):
    """Utility function to optionally pluralize words based on the value of n.
    """
    if n == 1:
        return ''
    else:
        return 's'


def decimateSig(arr):
    """function to decimate signal"""
    return spsig.decimate(arr, 5)


def decimateSig2(arr):
    """function to decimate signal"""
    return spsig.decimate(arr, 6)

# IBL shift code
# https://github.com/int-brain-lab/ibl-neuropixel/blob/fa5488d2a96a9f14fc3b86fbe4ef4091329213f6/src/neurodsp/fourier.py
def channel_shift(data, sample_shifts):
    """
    GPU Shifts channel signals via a Fourier transform to correct for different sampling times
    :param data: cupy array with shape (n_channels, n_times)
    :param sample_shifts: channel shifts, cupy array with shape (n_channels)
    :return: Aligned data, cupy array with shape (n_channels, n_times)
    """
    data = cp.array(data)
    sample_shifts = cp.array(sample_shifts)
    
    n_channels, n_times = data.shape

    dephas = cp.tile(- 2 * np.pi / n_times * cp.arange(n_times), (n_channels, 1))
    dephas += 2 * np.pi * (dephas < - np.pi)  # angles in the range (-pi,pi)
    dephas = cp.exp(1j * dephas * sample_shifts[:, cp.newaxis])

    data_shifted = cp.real(cp.fft.ifft(cp.fft.fft(data) * dephas))

    return cp.asnumpy(data_shifted)

# variable to calculate shift
# CHANGE this if needed
shift = np.tile(np.linspace(-1,0,32),16)

# this need to be changed for each animal
subsamplingfactor = 30
dirname = r'X:\Winny\Data\TR15_512ch_10-10-2023'
rawfname = 'TR15_512ch_10-10-2023_231010_163051'
opdirname = r'X:\Winny\Data\test'
aname = 'TR15'
saveLFP = True
saveAnalog = True

#####
lfp_filename = os.path.join(dirname,aname+'-lfp.npy')
digIn_filename = os.path.join(dirname, aname+'-digIn.npy')
analogIn_filename = os.path.join(dirname, aname+'-analogIn.npy')
analog_in = None
dig_in = None
amp_data_mmap = None
files = natsorted(glob.glob(os.path.join(dirname,rawfname,'*.rhd')))
for i, filename in enumerate(files):
    filename = os.path.basename(filename)
    if i==0:
        print("\n ***** Loading: " + filename)
        ts, amp_data, dig_in, analog_in, fs = read_data(os.path.join(dirname,rawfname,filename))
        amp_data_n  = []
        for c in range(amp_data.shape[0]):
            amp_data_n.append(np.array(channel_shift(np.array([amp_data[c]]), np.array([shift[c]]))[0] - 32768, dtype=np.int16))
        del amp_data
        amp_data_n = np.array(amp_data_n)
        arr = np.memmap(os.path.join(opdirname, filename[:-4]+'_shifted.bin'), dtype='int16', mode='w+', shape=amp_data_n[:256,:].T.shape)
        arr[:] = amp_data_n[:256,:].T
        del arr
        if saveLFP:
            # convert microvolts for lfp conversion
            digIn_ts_mmap = ts
            amp_data_n = np.multiply(0.195,  amp_data_n, dtype=np.float32)
            print("REAL FS = " + str(1./np.nanmedian(np.diff(ts))))
            size = amp_data_n.shape[1]
            ind = np.arange(0,size,subsamplingfactor)
            ts = ts[ind]
            fs = fs/float(subsamplingfactor)
            amp_ts_mmap = ts
            starts = amp_ts_mmap[-1]+1./fs
            amp_data_n = np.apply_along_axis(decimateSig,1,amp_data_n)
            amp_data_n = np.apply_along_axis(decimateSig2,1,amp_data_n)
            amp_data_mmap = amp_data_n
            amp_ts_mmap = ts
            del amp_data_n
    else:
        print("\n ***** Loading: " + filename)
        ts, amp_data, digIN, analogIN, fs = read_data(os.path.join(dirname,rawfname,filename))    
        amp_data_n  = []
        for c in range(amp_data.shape[0]):
            amp_data_n.append(np.array(channel_shift(np.array([amp_data[c]]), np.array([shift[c]]))[0] - 32768, dtype=np.int16))
        del amp_data
        amp_data_n = np.array(amp_data_n)
        arr = np.memmap(os.path.join(opdirname, filename[:-4]+'_shifted.bin'), dtype='int16', mode='w+', shape=amp_data_n[:256,:].T.shape)
        arr[:] = amp_data_n[:256,:].T
        del arr
        if saveLFP:
            # convert microvolts for lfp conversion
            amp_data_n = np.multiply(0.195,  amp_data_n, dtype=np.float32)
            print("REAL FS = " + str(1./np.nanmedian(np.diff(ts))))
            size = amp_data_n.shape[1]
            startind = np.where(ts>=starts)[0][0]
            ind = np.arange(startind,size,subsamplingfactor)
            ts = ts[ind]
            fs = fs/float(subsamplingfactor)
            starts = ts[-1]+1./fs
            amp_data_n = np.apply_along_axis(decimateSig,1,amp_data_n[:,startind:])
            amp_data_n = np.apply_along_axis(decimateSig2,1,amp_data_n)
            amp_data_mmap = np.concatenate((amp_data_mmap, amp_data_n), 1)
            dig_in = np.array(np.concatenate((dig_in, digIN), 1), dtype='uint8')
            if saveAnalog:
                analog_in = np.array(np.concatenate((analog_in, analogIN), 1), dtype=np.float32)
if saveLFP:
    np.save(lfp_filename, amp_data_mmap)
    np.save(digIn_filename, dig_in.T)
    if saveAnalog:
        np.save(analogIn_filename, analog_in.T)