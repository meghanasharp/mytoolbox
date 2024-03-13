#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 07:57:52 2024

@author: msharp
"""

import requests
import os

# %% define a function to save the files

def download_file(url, save_path):
    response = requests.get(url)
    with open(save_path, 'wb') as f:
        f.write(response.content)

# %% get the urls of each velocity pair
base_url = "https://nsidc.org/apps/itslive-search/velocities/urls/"

# Thwaites
params = {
  "bbox": "-107.562,-74.9689,-103.483,-75.1418",
  "start": "2022-01-01",
  "end": "2022-12-31",
  "percent_valid_pixels": 20, # percent of valid glacier pixels
  "min_interval": 7, # days of separation between image pairs
  "max_interval": 30, # days of separation between image pairs
  "version": 2, # version 1 requires an EDL access token header in the request
}
# This will return a list of NetCDF files in AWS S3 that can be accessed
# in the cloud or externally
velocity_pairs = requests.get(base_url, params=params)

print('Found ' + str(len(velocity_pairs.url)) + ' velocity pairs')

# %% 

# folder to save the files
nc_savefolder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs_TWITandTEIS/from_requests/nc/2022/'
png_savefolder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs_TWITandTEIS/from_requests/png/2022/'

# Check if the request was successful
if velocity_pairs.status_code == 200:
    
    
    # Iterate over the URLs in velocity_pairs
    for url in velocity_pairs.json():
        
        # Extract URLs from the pair
        nc_url = url['url']
        png_url = url['browse_image']
        
        # Extract filenames from URLs
        nc_filename = nc_url.split('/')[-1]
        png_filename = png_url.split('/')[-1]
        
        # Decide where to save the files
        nc_save_path = os.path.join(nc_savefolder, nc_filename)
        png_save_path = os.path.join(png_savefolder, png_filename)
        
        # Download the files
        download_file(nc_url, nc_save_path)
        download_file(png_url, png_save_path)
        print(f"Downloaded {nc_filename} and {png_filename}")
else:
    print("Failed to fetch velocity pairs:", velocity_pairs.status_code)