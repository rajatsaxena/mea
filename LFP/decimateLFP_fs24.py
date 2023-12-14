#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 09:13:42 2023

@author: rajat
"""

import os
from load_intan_rhd import saveLFP

parent_dir = 'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound5'
animals = ['SWIL22', 'SWIL23', 'SWIL24', 'SWIL25', 'SWIL26']

parent_dir = 'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound4'
animals = ['SWIL18', 'SWIL19', 'SWIL20', 'SWIL15']

parent_dir = r'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound3'
animals = ['SWIL11', 'SWIL12', 'SWIL13', 'SWIL105']

parent_dir = r'Y:\Research\SPrecordings\Rajat_Data\Data-SWIL\SWILRound2'
animals = ['SWIL8', 'SWIL10'] #, 'SWIL4-TD1', 'SWIL5-TD2', 'SWIL6-TD3', 'SWIL7-TD4'] 

for i in range(len(animals)):
    print("Procesing: " + str(animals[i]))
    saveLFP(animals[i], os.path.join(parent_dir,animals[i],'RawData'), os.path.join(parent_dir,animals[i]))