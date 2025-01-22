#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 17:08:10 2022

Import Picarro raw data file.

Function adapted for Picarro L2140i - HKDS2092

@author: daniele
"""

def import_picarro_raw_data(Picarro_data_path, op_mode, verbose = False):
    """
    Import raw data from Picarro file. Function willl load all available files 
    in Picarro_data_path. 
    col_nums variable must be checked against Picarro raw data file.
    
    Parameters
    ----------
    Picarro_data_path : string
        Path of the folder containing picarro data.  
    op_mode : TYPE
        'normal' should work in normal mode with L2140-i
        '17O' should work in O-17 mode with L2140-i
    verbose : boolean, default is false
        Returns information during loading of the data if set to True

    Returns
    -------
    Picarro data as a pandas dataframe.

    """
    import pandas as pd
    import numpy as np
    import os
    # Mode select
    col_names = ['DATE', 'TIME', 'H2O', 'Delta_18_16', 'Delta_D_H', 'ValveMask']
    col_names_extra = ['CavityPressure', 'CavityTemp', 'WarmBoxTemp', 'DasTemp', 'EtalonTemp', 'MPVPosition', 'OutletValve']
    #col_names = col_names + col_names_extra
    if op_mode  == 'normal':
        # Variable position in raw data file
        #col_nums = (0,1,2,8,15,16,17,21)
        pass
    elif op_mode == '17O':
        # Variable position in raw data file
        #col_nums = (0,1,2,8,9,10,11,12,14,15,16,17,18)
        col_names = col_names + ['Delta_17_16']
    # List files in directory
    data_filenames = os.listdir(Picarro_data_path)
    Picarro_data = pd.DataFrame()     
    # Import files
    for file in data_filenames:
        if file[-4:] == '.dat': # skip non-dat files
            if verbose:
                print('Reading: ', Picarro_data_path + file)
            df_dummy = pd.read_csv(Picarro_data_path + file,
                                   delim_whitespace = True,
                                   index_col = (0),
                                   parse_dates=[[0, 1]],
                                   #infer_datetime_format=True ,
                                   na_values=['NAN'],
                                   usecols=col_names,
                                   dtype={'FRAC_DAYS_SINCE_JAN1': np.float64})
            Picarro_data = pd.concat([Picarro_data, df_dummy], axis = 0)
    # Sort observation by index
    Picarro_data.index.name='Date'
    Picarro_data.sort_values(axis = 0, by = 'Date', inplace=True) # This now sorts in date order 
    #Picarro_data.index = Picarro_data.index.tz_localize('UTC') # Set timezone to UTC
    del(df_dummy)
    return Picarro_data