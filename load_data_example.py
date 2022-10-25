#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run this script to automatically load the Picarro raw data inside the directories of interest
"""

from import_picarro_raw_data import import_picarro_raw_data


# Folders of interest, example
paths = [
        #'../Picarro Data/2022/10/10/',
        '../Picarro Data/2022/10/11/',
        '../Picarro Data/2022/10/12/',
        '../Picarro Data/2022/10/12/',
        '../Picarro Data/2022/10/13/'
         ]# Directory for daily data

#% List files and import
for file_path in paths:
    if 'Picarro_data' in locals():
        new_Picarro_data = import_picarro_raw_data(file_path, '17O') # Enable 17O mode reading
        Picarro_data = pd.concat([Picarro_data, new_Picarro_data])
        del(new_Picarro_data)
    else:
        Picarro_data = import_picarro_raw_data(file_path, '17O') # Enable 17O mode reading

