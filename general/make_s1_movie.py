#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 20 21:08:26 2025

@author: msharp
"""
# %%
import os
import matplotlib.pyplot as plt
import numpy as np

import matplotlib
#matplotlib.use('Agg') # Non-interactive backend for file saving
#matplotlib.use("inline")

import imageio
from PIL import Image

# %%
make_pngs = 0
make_mov = 1

# %%

img_folder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/Clipped_Data/'

files = [fn for fn in os.listdir(img_folder) if fn.endswith('.tiff') and not fn.startswith('._')]

files = sorted(files, reverse=True)

# %%

if make_pngs:
        
    fig, ax = plt.subplots()
    
    for img_fn in files:
        print(f"Making {img_fn}...")
        
        img = plt.imread(img_folder + img_fn)
        
        img_date = int(img_fn[14:22]) 
        
        ax.imshow(img, vmin= 0, vmax = 400, cmap = 'gray')
        
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(img_date)
                     
        #plt.show()
    
        out_folder = img_folder + '/PNGs/'
        plt.savefig(out_folder + img_fn[:-4] + 'png', dpi = 300)

# %% make GIF

if make_mov:
    
    # List of your PNG files in order
    png_folder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/Clipped_Data/PNGs/'
    fn_png = [fn for fn in os.listdir(png_folder) if fn.endswith('.png') and not fn.startswith('._')]
    
    fn_png = sorted(fn_png, reverse=False)

    out_gif = png_folder + '20141010_20260329.mp4'
    
    # Read and save as GIF
    with imageio.get_writer(out_gif, fps=20) as writer:
        for fn in fn_png:
            print(f"Writing {fn}...")
            
            image = imageio.imread(png_folder + fn)
            
            
            #im = Image.open(png_folder + fn)

            #crop_box = (2000, 2000, 7000, 7000)  # example values
            #im_cropped = im.crop(crop_box)
            
            # Convert to numpy for imageio
            #image = np.array(im_cropped)
            
            writer.append_data(image)
    
# %%
