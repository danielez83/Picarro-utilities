#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Set of functions
Created on Fri Jul 19 12:19:30 2024

@author: daniele
"""

#%% Add path to Python
# import os
# import sys
# #cur_dir = os.path.dirname(__file__)
# #sys.path.append(os.path.join(cur_dir, "/functions"))

# # Imports
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import matplotlib.dates as mdates
# plt.style.use('seaborn-v0_8-whitegrid')


# General settings
#date_OI = '01.06.2021'
#correct_humidity = True

#%% import standard function
def import_std(my_paths, my_start_date, my_stop_date):
    from import_picarro_raw_data import import_picarro_raw_data
    import numpy as np
    import pandas as pd
    #List files and import
    for file_path in my_paths:
        if 'Picarro_data' in locals():
            new_Picarro_data = import_picarro_raw_data(file_path, 'normal')
            Picarro_data = pd.concat([Picarro_data, new_Picarro_data])
            del(new_Picarro_data)
        else:
            Picarro_data = import_picarro_raw_data(file_path, 'normal')
    
    # Subset data 
    start_winOI = pd.to_datetime(my_start_date, dayfirst = True)
    stop_winOI  = pd.to_datetime(my_stop_date, dayfirst = True)
    Picarro_date_mask = (Picarro_data.index > start_winOI) & (Picarro_data.index < stop_winOI)
    return Picarro_data[Picarro_date_mask]

#%% Humidity correction function
def humidity_correction(H2O_level_, uncorr_d18O, uncorr_dD, H2O_tie_point, correction_code):
    import numpy as np
    import pandas as pd
    
    # Correction functions
    if correction_code == 'Ramp_23.07.2024_bermuda':
        # off + a*exp(b*x)
        correction_d18O = lambda H2O_: 0.653065*np.exp(-0.000272*H2O_)
        correction_dD   = lambda H2O_: 34.338931*np.exp(-0.002426*H2O_)
    elif correction_code == 'Ramp_23.07.2024_m40':
        # off + a*exp(b*x)
        correction_d18O = lambda H2O_: 5.866788*np.exp(-0.001207*H2O_)
        correction_dD   = lambda H2O_: 34.092308*np.exp(-0.001219*H2O_)
    elif correction_code == 'Ramp_25.07.2024_m30':
        # off + (a/x) + (b*x)
        correction_d18O = lambda H2O_: (2584.716976/H2O_) + (0.000017*H2O_)
        correction_dD   = lambda H2O_: (25198.082048/H2O_) + (0.000408*H2O_)
    elif correction_code == 'Ramp_26.07.2024_bermuda':
        # off + a*exp(b*x)
        correction_d18O = lambda H2O_: 0.948753*np.exp(-0.000542*H2O_)
        correction_dD   = lambda H2O_: 11.205972*np.exp(-0.000770*H2O_)
    elif correction_code == 'Ramp_26.07.2024_m40':
        # off + a*exp(b*x)
        correction_d18O = lambda H2O_: 3.028881*np.exp(-0.000797*H2O_)
        correction_dD   = lambda H2O_: 20.739064*np.exp(-0.000829*H2O_)
    elif correction_code == 'Ramp_28.07.2024_bermuda':
        # off + a*exp(b*x)
        correction_d18O = lambda H2O_: 1.306959*np.exp(-0.000461*H2O_)
        correction_dD   = lambda H2O_: 16.223616*np.exp(-0.000987*H2O_)
    elif correction_code == 'Ramp_28.07.2024_B-Slap':
        # off + a*exp(b*x)
        correction_d18O = lambda H2O_: 2.384120*np.exp(-0.000531*H2O_)
        correction_dD   = lambda H2O_: 53.347717*np.exp(-0.001394*H2O_)
        
        
    # Estimate offset at  H2O level
    offset_d18O = correction_d18O(H2O_tie_point)
    offset_dD = correction_dD(H2O_tie_point)
    d18O_corrected = uncorr_d18O+correction_d18O(H2O_level_) - offset_d18O
    dD_corrected = uncorr_dD+correction_dD(H2O_level_) - offset_dD
    return (d18O_corrected, dD_corrected)

#%% Raw Picarro data calibration function
def calibrate_data(d18O_Raw, dD_Raw, date_OI, correct_humidity, correct_code):
    import numpy as np
    import pandas as pd
    from scipy.stats import linregress   
    from WIFVOS_standards import standards
    # Configuration -----------------------------------------------------
    basepath = '../Picarro Data/HBDS2212/'
    calibration_dates_file  = 'CalibrationAndRamps_20240723.xlsx'
    standards_d18O  = {}
    standards_dD  = {}
    for key in standards:
        standards_d18O[key] = standards[key]['d18O']
        standards_dD[key] = standards[key]['dD']
    # -------------------------------------------------------------------
    
    # Load Calibration spreadsheet
    df_calib = pd.read_excel(calibration_dates_file, sheet_name="Calibrations")
    df_calib['Date'] = pd.to_datetime(df_calib['Date'], dayfirst=True)
    df_calib['Start'] = pd.to_datetime(df_calib['Start'], dayfirst=True)
    df_calib['Stop'] = pd.to_datetime(df_calib['Stop'], dayfirst=True)

    # Subset the data for calibration
    mask_date = df_calib['Date'] == pd.to_datetime(date_OI, dayfirst=True)
    
    # Preallocate arrays
    d18O_meas = np.array([])
    d18O_true = np.array([])
    dD_meas = np.array([])
    dD_true = np.array([])
    H2O = np.array([])
    
    # Extract calibration data
    for std in df_calib[mask_date]['Standard'].unique():
        mask_std =  df_calib['Standard'] == std
        year = df_calib[mask_date & mask_std]['Start'].dt.year
        month = df_calib[mask_date & mask_std]['Start'].dt.month
        day = df_calib[mask_date & mask_std]['Start'].dt.day
        # Generate path
        strbuff = '%s%04d/%02d/%02d/' % (basepath, year.iloc[0], month.iloc[0], day.iloc[0])
        paths = [strbuff] # Directory for daily data
        start = df_calib[mask_date & mask_std]['Start'].iloc[0].strftime('%d-%m-%Y %H:%M')
        stop = df_calib[mask_date & mask_std]['Stop'].iloc[0].strftime('%d-%m-%Y %H:%M')
        curr_df = import_std(paths, start, stop)
        H2O_level = curr_df['H2O'].mean()
        # Here goes the humidity corection
        if correct_humidity:
                data_corrected = humidity_correction(curr_df['H2O'], 
                                                     curr_df['Delta_18_16'], 
                                                     curr_df['Delta_D_H'], 
                                                     H2O_level,
                                                     correct_code)
                d18O_corrected = data_corrected[0]
                dD_corrected = data_corrected[1]
        else:
            d18O_corrected = curr_df['Delta_18_16']
            dD_corrected = curr_df['Delta_D_H']
        
       
        
        d18O_meas = np.append(d18O_meas, d18O_corrected.mean())
        d18O_true = np.append(d18O_true, standards_d18O[std])
        dD_meas = np.append(dD_meas, dD_corrected.mean())
        dD_true = np.append(dD_true, standards_dD[std])
        H2O = np.append(H2O, curr_df['H2O'].mean())
    
    # Compute slopes and intercepts
    mod_d18O = linregress(d18O_meas, d18O_true)
    mod_dD = linregress(dD_meas, dD_true)
    model_structure = {'d18O' : mod_d18O, 'dD' : mod_dD, 'H2O' : H2O}
    return (d18O_Raw*model_structure['d18O'].slope + model_structure['d18O'].intercept,
            dD_Raw*model_structure['dD'].slope + model_structure['dD'].intercept)


#%% Humidity calibration
def calibrate_H2O(H2O_in, code):
    if code == 'H2O_23.07.2024':
        #slope = 0.7794593128208636
        #intercept = -408.02183176813924
        #calibrated_H2O = H2O_in*slope + intercept
        a = 2.71384315e-06  
        b = 7.27941261e-01
        c = -2.83703296e+02
        calibrated_H2O = a*H2O_in**2 + b*H2O_in + c
        
    elif code == 'H2O_25.07.2024':
        #slope = 0.7827297486870664
        #intercept = -246.1363592028656
        #calibrated_H2O = H2O_in*slope + intercept
        a = 1.07099683e-05
        b = 6.13045971e-01
        c = 1.43080137e+02
        calibrated_H2O = a*H2O_in**2 + b*H2O_in + c
    elif code == 'H2O_26.07.2024':
        #slope = 
        #intercept = 
        #calibrated_H2O = H2O_in*slope + intercept
        a = 5.45632805e-06  
        b = 6.91138071e-01
        c = -1.45274672e+01
        calibrated_H2O = a*H2O_in**2 + b*H2O_in + c
    elif code == 'H2O_all_DEPRECATED': # DON'T USE!
        a = 3.89852245e-06
        b = 7.11282411e-01
        c = -1.29986822e+02
        calibrated_H2O = a*H2O_in**2 + b*H2O_in + c
    elif code == 'H2O_30.07.2024': # GOOD
        slope = 0.8495326607880471
        intercept = 33.57303203475749
        calibrated_H2O = slope*H2O_in + intercept
    elif code == 'H2O_31.07.2024': # GOOD
        slope = 0.8407433060015952
        intercept = 52.31582324241481
        calibrated_H2O = slope*H2O_in + intercept
    elif code == 'H2O_All_HM40': # Excellent
        slope = 0.84423595009
        intercept = 51.28712608197
        calibrated_H2O = slope*H2O_in + intercept
    else: # No correction
        slope = 1
        intercept = 0
        calibrated_H2O = H2O_in*slope + intercept
        
    return calibrated_H2O

#%% Mass import of Picarro raw data
def mass_import_Picarro_data(paths,
                             start_date_str, 
                             stop_date_str,
                             start_win_rem = [],
                             stop_win_rem = [],
                             mode = 'normal'):
    from import_picarro_raw_data import import_picarro_raw_data
    import pandas as pd
    import numpy as np
    
    for file_path in paths:
        Picarro_data = import_picarro_raw_data(file_path, mode)
        
    # Subset data 
    start_winOI = pd.to_datetime(start_date_str, dayfirst = True)
    stop_winOI  = pd.to_datetime(stop_date_str, dayfirst = True)
    Picarro_date_mask = (Picarro_data.index > start_winOI) & (Picarro_data.index < stop_winOI)
    
    # Trim data
    if len(start_win_rem) > 0:
        for cur_win_start, cur_win_stop in zip(start_win_rem, stop_win_rem):
            start_winOI = pd.to_datetime(cur_win_start, dayfirst = True)
            stop_winOI  = pd.to_datetime(cur_win_stop, dayfirst = True)
            Picarro_date_mask_temp = (Picarro_data.index > start_winOI) & (Picarro_data.index < stop_winOI)
            Picarro_data[Picarro_date_mask_temp] = np.nan
            #Picarro_date_mask = Picarro_date_mask & np.logical_not(Picarro_date_mask_temp)
    
    #Picarro_data_subset = Picarro_data[Picarro_date_mask]
    return Picarro_data[Picarro_date_mask]


