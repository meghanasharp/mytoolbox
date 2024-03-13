#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 13:44:16 2024

@author: msharp
"""

import geopandas as gpd
from shapely.geometry import Point
from shapely.geometry import Polygon
import xarray as xr
import hvplot.pandas
import cartopy
import matplotlib as plt
from shapely import geometry


# %% Get ITS_LIVE catalog

itslive_catalog = gpd.read_file('https://its-live-data.s3.amazonaws.com/datacubes/catalog_v02.json')

# %% Search granules by a single point

def find_granule_by_point(input_point):
    '''returns url for the granule (zarr datacube) containing a specified point. point must be passed in epsg:4326
    '''
    catalog = gpd.read_file('https://its-live-data.s3.amazonaws.com/datacubes/catalog_v02.json')

    #make shapely point of input point
    p = gpd.GeoSeries([Point(input_point[0], input_point[1])],crs='EPSG:4326')
    #make gdf of point
    gdf = gdf = gpd.GeoDataFrame({'label': 'point', 
                                  'geometry':p})
    #find row of granule 
    granule = catalog.sjoin(gdf, how='inner')

    url = granule['zarr_url'].values[0]
    return url

# search for a granule. Should return a point.
#url = find_granule_by_point([95.180191, 30.645973])

point = [-105.614384, -74.930176]
url = find_granule_by_point(point)


# %% Open the datacube as an xarray.Dataset

def read_in_s3(http_url, chunks = 'auto'):
    #s3_url = http_url.replace('http','s3')   <-- as of Fall 2023, can pass http urls directly to xr.open_dataset()
    #s3_url = s3_url.replace('.s3.amazonaws.com','')

    datacube = xr.open_dataset(http_url, engine = 'zarr',
                                #storage_options={'anon':True},
                                chunks = chunks)

    return datacube

# open the datacube
dc = read_in_s3(url)

#inpath = "/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/"
#filename = "test_arrows.nc"
#dc = xr.open_dataset(inpath+filename)

# %% Visualize the footprint

def get_bounds_polygon(input_xr):
    
    xmin = input_xr.coords['x'].data.min()
    xmax = input_xr.coords['x'].data.max()

    ymin = input_xr.coords['y'].data.min()
    ymax = input_xr.coords['y'].data.max()

    pts_ls = [(xmin, ymin), (xmax, ymin),(xmax, ymax), (xmin, ymax), (xmin, ymin)]

    crs = f"epsg:{input_xr.mapping.spatial_epsg}"

    polygon_geom = Polygon(pts_ls)
    polygon = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[polygon_geom]) 
    
    return polygon

# get the bounds
bbox_dc = get_bounds_polygon(dc)

# visualize
# Note: this didn't work... I tried the examples online too and it seems to be a problem with those
# as well... i even tried downgrade hvplot
#poly =  bbox_dc.to_crs('EPSG:4326').hvplot(geo=True,legend=True,alpha=0.3, tiles='ESRI', color='red')

# %%

# To make it easier for now, we will only plot a few hundred time steps instead of the full time series
dc_sub = dc.isel(mid_date = slice(0,50))

# set a maximum value for the plot
dc_sub.v.mean(dim='mid_date').plot(vmax=300)


# %% Quiver plot

dc_test = dc.isel(mid_date=2)

# Figure size
size = (7, 8)

# Defining the figure
fig = plt.figure(figsize=size, facecolor='w', 
                 edgecolor='k')

#dc_test.plot.quiver(x = 'x', y= 'y', u = 'vx', v = 'vy')



