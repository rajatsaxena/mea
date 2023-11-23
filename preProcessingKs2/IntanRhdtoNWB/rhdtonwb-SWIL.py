#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 01:02:50 2023

@author: rajat
"""
import pynwb
from pynwb import NWBFile
from ConvertIntanToNWB import *

# create data module for each animal

# swil8
swil8 = pynwb.file.Subject(
  age="P6M3D",
  subject_id="swil8",
  sex="M",
  species="Mus musculus",
  description = "batch 2",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound2/SWIL8/RawData/SWIL8_210715_191215.rhd',
                   nwb_filename='swil8merged.nwb',
                   session_description='fully interleaved + blocked novel',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil8,
                   manual_start_time=None)

# SWIL10
swil10 = pynwb.file.Subject(
  age="P6M8D",
  subject_id="swil10",
  sex="M",
  species="Mus musculus",
  description = "batch 2",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound2/SWIL10/RawData/SWIL10_210717_170731.rhd',
                   nwb_filename='swil10merged.nwb',
                   session_description='fully interleaved + blocked novel',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil10,
                   manual_start_time=None)

# SWIL105
swil105 = pynwb.file.Subject(
  age="P6M14D",
  subject_id="swil105",
  sex="M",
  species="Mus musculus",
  description = "batch 3",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound3/SWIL105/RawData/SWIL10_220415_160508.rhd',
                   nwb_filename='swil105merged.nwb',
                   session_description='interleaved blocks',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil105,
                   manual_start_time=None)

# SWIL11
swil11 = pynwb.file.Subject(
  age="P6M12D",
  subject_id="swil11",
  sex="M",
  species="Mus musculus",
  description = "batch 3",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound3/SWIL11/RawData/SWIL11_220413_171108.rhd',
                   nwb_filename='swil11merged.nwb',
                   session_description='interleaved blocks',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil11,
                   manual_start_time=None)

# SWIL11r2
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound3/SWIL11/SWIL11_220413_175914/SWIL11_220413_175914.rhd',
                   nwb_filename='swil11merged_v2.nwb',
                   session_description='interleaved blocks',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil11,
                   manual_start_time=None)

# SWIL12
swil12 = pynwb.file.Subject(
  age="P7M2D",
  subject_id="swil12",
  sex="M",
  species="Mus musculus",
  description = "batch 3",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound3/SWIL12/RawData/SWIL12_220602_183004.rhd',
                   nwb_filename='swil12merged.nwb',
                   session_description='interleaved blocks',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil12,
                   manual_start_time=None)

# convert SWIL12r2
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound3/SWIL12/SWIL12_220602_214406/SWIL12_220602_214406.rhd',
                   nwb_filename='swil12merged_v2.nwb',
                   session_description='interleaved blocks',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil12,
                   manual_start_time=None)

# SWIL13
swil13 = pynwb.file.Subject(
  age="P7M11D",
  subject_id="swil13",
  sex="M",
  species="Mus musculus",
  description = "batch 3",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound3/SWIL13/RawData/SWIL13_220714_171703.rhd',
                   nwb_filename='swil13merged.nwb',
                   session_description='interleaved blocks',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil13,
                   manual_start_time=None)

# SWIL15
swil15 = pynwb.file.Subject(
  age="P7M24D",
  subject_id="swil15",
  sex="M",
  species="Mus musculus",
  description = "batch 4",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound4/SWIL15/RawData/SWIL15_220827_152048.rhd',
                   nwb_filename='swil15merged.nwb',
                   session_description='interleaved blocks',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil15,
                   manual_start_time=None)

# SWIL18
swil18 = pynwb.file.Subject(
  age="P8M2D",
  subject_id="swil18",
  sex="M",
  species="Mus musculus",
  description = "batch 4",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound4/SWIL18/RawData/SWIL18_220908_163938.rhd',
                   nwb_filename='swil18merged.nwb',
                   session_description='interleaved blocks',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil18,
                   manual_start_time=None)

# SWIL19
swil19 = pynwb.file.Subject(
  age="P7M27D",
  subject_id="swil19",
  sex="M",
  species="Mus musculus",
  description = "batch 4",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound4/SWIL19/RawData/SWIL19_220901_160502.rhd',
                   nwb_filename='swil19merged.nwb',
                   session_description='interleaved blocks',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil19,
                   manual_start_time=None)

# SWIL20
swil20 = pynwb.file.Subject(
  age="P8M13D",
  subject_id="swil20",
  sex="M",
  species="Mus musculus",
  description = "batch 4",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound4/SWIL20/RawData/SWIL20_220922_163606.rhd',
                   nwb_filename='swil20merged.nwb',
                   session_description='interleaved blocks',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil20,
                   manual_start_time=None)


# SWIL22
swil22 = pynwb.file.Subject(
  age="P7M1D",
  subject_id="swil22",
  sex="M",
  species="Mus musculus",
  description = "batch 5",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound5/SWIL22/RawData/SWIL22_230112_155201.rhd',
                   nwb_filename='swil22merged.nwb',
                   session_description='interleaved blocks',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil22,
                   manual_start_time=None)

# SWIL23
swil23 = pynwb.file.Subject(
  age="P7M14D",
  subject_id="swil23",
  sex="M",
  species="Mus musculus",
  description = "batch 5",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound5/SWIL23/RawData/SWIL23_230125_155701.rhd',
                   nwb_filename='swil23merged.nwb',
                   session_description='interleaved blocks',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil23,
                   manual_start_time=None)

# SWIL24
swil24 = pynwb.file.Subject(
  age="P7M8D",
  subject_id="swil24",
  sex="M",
  species="Mus musculus",
  description = "batch 5",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound5/SWIL24/RawData/SWIL24_230119_160102.rhd',
                   nwb_filename='swil24merged.nwb',
                   session_description='interleaved blocks',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil24,
                   manual_start_time=None)

# SWIL25
swil25 = pynwb.file.Subject(
  age="P7M3D",
  subject_id="swil25",
  sex="M",
  species="Mus musculus",
  description = "batch 5",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound5/SWIL25/RawData/SWIL25_230114_162806.rhd',
                   nwb_filename='swil25merged.nwb',
                   session_description='interleaved blocks',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil25,
                   manual_start_time=None)

# SWIL26
swil26 = pynwb.file.Subject(
  age="P6M26D",
  subject_id="swil26",
  sex="M",
  species="Mus musculus",
  description = "batch 5",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound5/SWIL26/RawData/SWIL26_230113_155201.rhd',
                   nwb_filename='swil26merged.nwb',
                   session_description='interleaved blocks',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil26,
                   manual_start_time=None)

# SWIL3
swil3 = pynwb.file.Subject(
  age="P6M19D",
  subject_id="swil3",
  sex="M",
  species="Mus musculus",
  description = "batch 1",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound1/SWIL3_Remi/RawData/swil3_200521_181723.rhd',
                   nwb_filename='swil3merged.nwb',
                   session_description='interleaved blocks',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil3,
                   manual_start_time=None)

# SWIL4
swil4 = pynwb.file.Subject(
  age="P6M28D",
  subject_id="swil4",
  sex="M",
  species="Mus musculus",
  description = "batch 2",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound2/SWIL4-TD1/RawData/TD1_210420_183050.rhd',
                   nwb_filename='swil4merged.nwb',
                   session_description='fully interleaved',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil4,
                   manual_start_time=None)

# SWIL5
swil5 = pynwb.file.Subject(
  age="P7M1D",
  subject_id="swil5",
  sex="M",
  species="Mus musculus",
  description = "batch 2",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound2/SWIL5-TD2/RawData/TD2_210423_182211.rhd',
                   nwb_filename='swil5merged.nwb',
                   session_description='fully interleaved',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil5,
                   manual_start_time=None)

# swil6
swil6 = pynwb.file.Subject(
  age="P7M3D",
  subject_id="swil6",
  sex="M",
  species="Mus musculus",
  description = "batch 2",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound2/SWIL6-TD3/RawData/TD3_210430_170216.rhd',
                   nwb_filename='swil6merged.nwb',
                   session_description='fully interleaved',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil6,
                   manual_start_time=None)

# swil7
swil7 = pynwb.file.Subject(
  age="P6M29D",
  subject_id="swil7",
  sex="M",
  species="Mus musculus",
  description = "batch 2",
  )
convert_to_nwb(settings_filename=None,
                   intan_filename='/media/rajat/mcnlab_store2/Research/SPrecordings/Rajat_Data/Data-SWIL/SWILRound2/SWIL7-TD4/RawData/TD4_210427_173206.rhd',
                   nwb_filename='swil7merged.nwb',
                   session_description='fully interleaved',
                   blocks_per_chunk=2000,
                   use_compression=True,
                   compression_level=4,
                   lowpass_description=None,
                   highpass_description=None,
                   merge_files=True,
                   subject=swil7,
                   manual_start_time=None)
