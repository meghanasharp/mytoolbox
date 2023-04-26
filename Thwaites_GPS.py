#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 15:20:33 2023

@author: msharp
"""

# add OSU Glaciology GPS Functions to my path
import sys
sys.path.append("/Users/msharp/GitHub/Toolbox")

# import OS module
import os
# importing the zipfile module
from zipfile import ZipFile
#home made function from OSU-Glaciology Toolbox
# this require to have the function in the same folder
import GPS_functions as gf

#%%

#define my paths!
directory = "/Volumes/Extreme_SSD/PhD_MS/TARSAN/GPS/PPP/"
foldernames = ["PPP_16242_TAR2_netrs_N09_rinex_daily_2223/",
              "PPP_23174_TAR3_netrs_N07_rinex_daily_2223_part1/",
              "PPP_23174_TAR3_netrs_N07_rinex_daiy_2223_part2/",
              "PPP_25474_TAR1_netrs_N15_rinex_daily_2223/",
              "PPP_25479_TAR4_netrs_N12_rinex_daily_2223_part1/",
              "PPP_25479_TAR4_netrs_N12_rinex_daily_2223_part2/"]

for fn in foldernames:

    path_to_folder = directory + fn
    path_to_save = directory + "Concatenated/"


    output_filename = fn[:-1] + "_concat"
    
    #note: figure out what the timezone for the GPS data is!
    local_utc = 0
    
    df_base = gf.concat_gps_csv_files(path_to_folder = path_to_folder, 
                      path_to_save = path_to_save, 
                      output_filename = output_filename, 
                      local_utc = local_utc)
