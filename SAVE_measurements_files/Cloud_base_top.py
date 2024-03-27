# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os as os
import xarray as xr
import numpy as np
import pandas as pd
import datetime
import scipy.io
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


path_to_cloud_base = "/home/delbeke/Documents/LIMA/SAVE_insitu/KIT_ceilometer/07/"
path_to_cloud_top = "/home/delbeke/Documents/LIMA/SAVE_insitu/KIT_ceilometer/UPS_cloud_top_radar/"
path_to_MESONH = "/home/delbeke/Documents/LIMA/DACCIWA_XXL_LIMA_BACK_RESULT/NPZ_FILES/"

in_situ_file_cloud_base = "Save_KIT_CM_20160703.nc"
in_situ_file_cloud_top = "cloud_top.csv"
in_situ_file_MESONH_MRC3D = "BCBackground_MRC3D_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_MRC3D)
mesonh_mrc3d = data["name"][1:]
in_situ_file_MESONH_MRC = "BCBackground_XXL_MRC_aerooff.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_MRC)
mesonh_mrc = data["name"][1:]
in_situ_file_MESONH_alt = "zhat.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_alt)
mesonh_alt = data["name"][:]

in_situ_products_cloud_base = xr.open_dataset(path_to_cloud_base + in_situ_file_cloud_base)
cloud_base = in_situ_products_cloud_base.cbh.values[:,0]
time_cloud_base = in_situ_products_cloud_base.time.values[:]
list_time_cloud_base = []
for iter_surf in range(len(time_cloud_base)):
    timestamp = pd.Timestamp(time_cloud_base[iter_surf])
    #print(timestamp)
    hour = timestamp.hour
    minute = timestamp.minute
    minute_decimal = minute/60
    date_surf = hour +  minute_decimal
    list_time_cloud_base.append(date_surf)
    #print(hour,minute,minute/60)
    
df = pd.read_csv (path_to_cloud_top + in_situ_file_cloud_top)
cloud_top = df['altitude (m)'].tolist()
list_time_cloud_top = df['hour'].tolist()

npixels_x = 240
npixels_y = 240
npixels_by_layer = npixels_x * npixels_y
list_time_occupied_layer = []
for iter_time_meso in range(np.shape(mesonh_mrc3d)[0]):
    list_altitude_occupied_layer = []
    for iter_altitude_meso in range(np.shape(mesonh_mrc3d)[1]):
        n_cloud_pixels = 0
        for iter_pixel_x in range(npixels_x):
            for iter_pixel_y in range(npixels_y):
                mrc_pixel = mesonh_mrc3d[iter_time_meso,iter_altitude_meso,iter_pixel_x,iter_pixel_y]
                if mrc_pixel >= 0.05:
                    n_cloud_pixels += 1
                    #print(mesonh_mrc3d[iter_time_meso,iter_altitude_meso,iter_pixel_x,iter_pixel_y])
        fraction_occupied_layer = n_cloud_pixels/npixels_by_layer
        pourcentage_occupied_layer = fraction_occupied_layer*100
        #print(pourcentage_occupied_layer)
        list_altitude_occupied_layer.append(pourcentage_occupied_layer)
        for iter_p in range(len(list_altitude_occupied_layer)):
            if list_altitude_occupied_layer[iter_p] == 0.0:
                list_altitude_occupied_layer[iter_p] = np.nan
    list_time_occupied_layer.append(list_altitude_occupied_layer)

list_cloud_time = []
for iter_time in range(24):
    list_cloud_time.append(iter_time)
    

list_base_cloud = []
for iter_time_bc_cloud in range(24):
    mean_mrc = mesonh_mrc[iter_time_bc_cloud,:]
    mean_mrc[mean_mrc>=0.05]
    #print("aaaaaaaaaaaaaaaaaaaaaaaaaaaa")
    #print(mean_mrc)
    base_cloud = 0 
    count = 0
    for iter_alt_bc_cloud in range(142):
        if count == 0 and mean_mrc[iter_alt_bc_cloud]>=0.05:
            base_cloud = mesonh_alt[iter_alt_bc_cloud]
            count += 1
        else:
            continue
    #print(base_cloud)
    list_base_cloud.append(base_cloud)



list_top_cloud = []
for iter_time_bc_cloud in range(24):
    mean_mrc = mesonh_mrc[iter_time_bc_cloud,:]
    mean_mrc[mean_mrc>=0.05]
    top_cloud = 0 
    for iter_alt_bc_cloud in range(142):
        count = 0
        if count == 0 and mean_mrc[iter_alt_bc_cloud]>=0.05:
            top_cloud = mesonh_alt[iter_alt_bc_cloud]
            count += 1
        else:
            continue
    list_top_cloud.append(top_cloud)


    
fig, ax = plt.subplots(figsize=(20,12))
plt.plot(list_cloud_time[2:-7],list_base_cloud[2:-7],color='blue',linewidth=4.0, label = 'Simulated CBH')
plt.plot(list_cloud_time[2:-7],list_top_cloud[2:-7],color='red',linewidth=4.0, label = 'Simulated CTH')
plt.scatter(list_time_cloud_base,cloud_base, color = 'k', label = 'Observed CBH')
plt.scatter(list_time_cloud_top,cloud_top, color = 'grey', label = 'Observed CTH')
cmap = plt.get_cmap('turbo')
for iter_time_plot in range(np.shape(mesonh_mrc3d)[0]):
    if iter_time_plot<2:
        continue
    else:
        list_time_plot = []
        for iter_alt in range(142):
            list_time_plot.append(iter_time_plot)
        cs = plt.scatter(list_time_plot,mesonh_alt,cmap=cmap, c=list_time_occupied_layer[iter_time_plot])
cbar = plt.colorbar(cs,pad = 0.01)
cbar.set_label('Cloud Fraction (%)', rotation=270,fontsize=20,labelpad=20)
cbar.ax.tick_params(labelsize=20)
plt.xlabel("Time UTC (h)", fontsize = 20)
plt.xticks(fontsize=20)
plt.xlim(0,18)
plt.ylabel("Altitude (m)", fontsize = 20 )
plt.yticks(fontsize=20)
plt.ylim(0,1400)
plt.grid(which='major', linestyle='-', linewidth='1', color='k')
plt.grid(which='minor', linestyle='-', linewidth='0.25', color='k')
plt.minorticks_on()
ax.annotate("", (6,1380), (6,1420),size = 100, arrowprops = dict(facecolor ='k', width=5, headlength = 0.1, headwidth=0.1))
ax.annotate("", (10,1380), (10,1420),size = 100, arrowprops = dict(facecolor ='k', width=5, headlength = 0.1, headwidth=0.1))
ax.annotate("", (17,1380), (17,1420),size = 100, arrowprops = dict(facecolor ='k', width=5, headlength = 0.1, headwidth=0.1))
ax.annotate("", (0.1,1400), (5.9,1400), size = 10, arrowprops = dict(facecolor ='r', width=5, headlength = 0.1, headwidth=0.1))
ax.annotate("", (6.1,1400), (9.9,1400), size = 10, arrowprops = dict(facecolor ='r', width=5, headlength = 0.1, headwidth=0.1))
ax.annotate("", (10.1,1400), (16.9,1400), size = 10, arrowprops = dict(facecolor ='r', width=5, headlength = 0.1, headwidth=0.1))
plt.text(2.0, 1420, 'jet phase', dict(size=20))
plt.text(6.75, 1420, 'stratus phase', dict(size=20))
plt.text(11.75, 1420, 'convective phase', dict(size=20))
plt.legend(fontsize=20)

