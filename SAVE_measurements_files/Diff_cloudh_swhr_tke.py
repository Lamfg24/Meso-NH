#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 19:40:19 2024

@author: delbeke
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

path_to_MESONH = "/home/delbeke/Documents/LIMA/DACCIWA_XXL_LIMA_BACK_RESULT/NPZ_FILES/"

in_situ_file_MESONH_MRC_aeroon = "BCBackground_XXL_MRC_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_MRC_aeroon)
mesonh_mrc_aeroon = data["name"][1:]

in_situ_file_MESONH_MRC_aerooff = "BCBackground_XXL_MRC_aerooff.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_MRC_aerooff)
mesonh_mrc_aerooff = data["name"][1:]

in_situ_file_MESONH_alt = "zhat.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_alt)
mesonh_alt = data["name"][:]

list_cloud_time = []
for iter_time in range(24):
    list_cloud_time.append(iter_time)
    

list_base_cloud_on = []
for iter_time_bc_cloud in range(24):
    mean_mrc_on = mesonh_mrc_aeroon[iter_time_bc_cloud,:]
    mean_mrc_on[mean_mrc_on>=0.05]
    #print("aaaaaaaaaaaaaaaaaaaaaaaaaaaa")
    #print(mean_mrc)
    base_cloud = 0 
    count = 0
    for iter_alt_bc_cloud in range(142):
        if count == 0 and mean_mrc_on[iter_alt_bc_cloud]>=0.05:
            base_cloud = mesonh_alt[iter_alt_bc_cloud]
            count += 1
        else:
            continue
    #print(base_cloud)
    list_base_cloud_on.append(base_cloud)



list_top_cloud_on = []
for iter_time_bc_cloud in range(24):
    mean_mrc_on = mesonh_mrc_aeroon[iter_time_bc_cloud,:]
    mean_mrc_on[mean_mrc_on>=0.05]
    top_cloud = 0 
    for iter_alt_bc_cloud in range(142):
        count = 0
        if count == 0 and mean_mrc_on[iter_alt_bc_cloud]>=0.05:
            top_cloud = mesonh_alt[iter_alt_bc_cloud]
            count += 1
        else:
            continue
    list_top_cloud_on.append(top_cloud)
    
    

list_base_cloud_off = []
for iter_time_bc_cloud in range(24):
    mean_mrc_off = mesonh_mrc_aerooff[iter_time_bc_cloud,:]
    mean_mrc_off[mean_mrc_off>=0.05]
    #print("aaaaaaaaaaaaaaaaaaaaaaaaaaaa")
    #print(mean_mrc)
    base_cloud = 0 
    count = 0
    for iter_alt_bc_cloud in range(142):
        if count == 0 and mean_mrc_off[iter_alt_bc_cloud]>=0.05:
            base_cloud = mesonh_alt[iter_alt_bc_cloud]
            count += 1
        else:
            continue
    #print(base_cloud)
    list_base_cloud_off.append(base_cloud)



list_top_cloud_off = []
for iter_time_bc_cloud in range(24):
    mean_mrc_off = mesonh_mrc_aerooff[iter_time_bc_cloud,:]
    mean_mrc_off[mean_mrc_off>=0.05]
    top_cloud = 0 
    for iter_alt_bc_cloud in range(142):
        count = 0
        if count == 0 and mean_mrc_off[iter_alt_bc_cloud]>=0.05:
            top_cloud = mesonh_alt[iter_alt_bc_cloud]
            count += 1
        else:
            continue
    list_top_cloud_off.append(top_cloud)
    
list_diff_base_cloud = []
for iter_diff in range(24):
    diff_base = list_base_cloud_on[iter_diff] - list_base_cloud_off[iter_diff]
    list_diff_base_cloud.append(diff_base)
    
list_diff_top_cloud = []
for iter_diff in range(24):
    diff_top = list_top_cloud_on[iter_diff] - list_top_cloud_off[iter_diff]
    list_diff_top_cloud.append(diff_top)
    
fig, ax = plt.subplots(figsize=(20,12))
plt.plot(list_cloud_time[2:-7],list_diff_base_cloud[2:-7],color='blue',linewidth=4.0,label = 'CBH')
plt.plot(list_cloud_time[2:-7],list_diff_top_cloud[2:-7],color='red',linewidth=4.0, linestyle = "--" ,label = 'CTH')
plt.xlabel("Time UTC (h)", fontsize = 20)
plt.xticks(fontsize=20)
plt.xlim(1,17)
plt.ylabel("Altitude Difference (m)", fontsize = 20 )
plt.yticks(fontsize=20)
plt.ylim(-25,15)
plt.grid(which='major', linestyle='-', linewidth='1', color='k')
plt.grid(which='minor', linestyle='-', linewidth='0.25', color='k')
plt.minorticks_on()
plt.legend(fontsize=20)


in_situ_file_MESONH_swhr_on = "BCBackground_XXL_DTRAD_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_swhr_on)
mesonh_swhr_aeroon = data["name"][1:]

in_situ_file_MESONH_swhr_off = "BCBackground_XXL_DTRAD_aerooff.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_swhr_off)
mesonh_swhr_aerooff = data["name"][1:]

list_diff_swhr = []
for iter_time in range(24):
    list_diff_swhr_alt = []
    for iter_alt in range(142):
        diff_tke_alt = mesonh_swhr_aeroon[iter_time,iter_alt] - mesonh_swhr_aerooff[iter_time,iter_alt]
        list_diff_swhr_alt.append(diff_tke_alt)
    list_diff_swhr.append(list_diff_swhr_alt)

fig, axs = plt.subplots(1,1,figsize=(15, 10))
cmap = plt.get_cmap('seismic')
for iter_time in range(6,17,2):
    list_plot_vertical = []
    for iter_time_plot in range(142):
        list_plot_vertical.append(iter_time)
    cs = plt.scatter(list_plot_vertical, mesonh_alt, c=list_diff_swhr[iter_time], cmap=cmap, linewidth = 5,vmin=-2,vmax=2.0)
cbar = plt.colorbar(cs,pad = 0.01)
cbar.set_label('SWHR Difference (K.$kd^{-1}$)', rotation=270,fontsize=20,labelpad=20)
cbar.ax.tick_params(labelsize=20)
plt.xlabel("Time UTC (h)", fontsize = 20)
plt.xticks(fontsize=20)
plt.xlim(5,17)
plt.ylabel("Altitude (m)", fontsize = 20 )
plt.yticks(fontsize=20)
plt.ylim(0,1400)
plt.annotate("", (6,1380), (6,1420),size = 100, arrowprops = dict(facecolor ='k', width=5, headlength = 0.1, headwidth=0.1))
plt.annotate("", (10,1380), (10,1420),size = 100, arrowprops = dict(facecolor ='k', width=5, headlength = 0.1, headwidth=0.1))
plt.annotate("", (17,1380), (17,1420),size = 100, arrowprops = dict(facecolor ='k', width=5, headlength = 0.1, headwidth=0.1))
plt.annotate("", (0.1,1400), (5.9,1400), size = 10, arrowprops = dict(facecolor ='r', width=5, headlength = 0.1, headwidth=0.1))
plt.annotate("", (6.1,1400), (9.9,1400), size = 10, arrowprops = dict(facecolor ='r', width=5, headlength = 0.1, headwidth=0.1))
plt.annotate("", (10.1,1400), (16.9,1400), size = 10, arrowprops = dict(facecolor ='r', width=5, headlength = 0.1, headwidth=0.1))
plt.text(6.75, 1420, 'stratus phase', dict(size=20))
plt.text(11.75, 1420, 'convective phase', dict(size=20))



in_situ_file_MESONH_tke_on = "BCBackground_XXL_TKET_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_tke_on)
mesonh_tke_aeroon = data["name"][1:]

in_situ_file_MESONH_tke_off = "BCBackground_XXL_TKET_aerooff.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_tke_off)
mesonh_tke_aerooff = data["name"][1:]

list_diff_tke = []
for iter_time in range(24):
    list_diff_tke_alt = []
    for iter_alt in range(142):
        diff_tke_alt = mesonh_tke_aeroon[iter_time,iter_alt] - mesonh_tke_aerooff[iter_time,iter_alt]
        list_diff_tke_alt.append(diff_tke_alt)
    list_diff_tke.append(list_diff_tke_alt)

fig, axs = plt.subplots(1,1,figsize=(15, 10))
cmap = plt.get_cmap('Spectral')
for iter_time in range(6,17,2):
    list_plot_vertical = []
    for iter_time_plot in range(142):
        list_plot_vertical.append(iter_time)
    cs = plt.scatter(list_plot_vertical, mesonh_alt, c=list_diff_tke[iter_time], cmap=cmap, linewidth = 5)
cbar = plt.colorbar(cs,pad = 0.01)
cbar.set_label('TKE Difference ($m^{2}.s^{-2}$)', rotation=270,fontsize=20,labelpad=20)
cbar.ax.tick_params(labelsize=20)
plt.xlabel("Time UTC (h)", fontsize = 20)
plt.xticks(fontsize=20)
plt.xlim(5,17)
plt.ylabel("Altitude (m)", fontsize = 20 )
plt.yticks(fontsize=20)
plt.ylim(0,1400)
plt.annotate("", (6,1380), (6,1420),size = 100, arrowprops = dict(facecolor ='k', width=5, headlength = 0.1, headwidth=0.1))
plt.annotate("", (10,1380), (10,1420),size = 100, arrowprops = dict(facecolor ='k', width=5, headlength = 0.1, headwidth=0.1))
plt.annotate("", (17,1380), (17,1420),size = 100, arrowprops = dict(facecolor ='k', width=5, headlength = 0.1, headwidth=0.1))
plt.annotate("", (0.1,1400), (5.9,1400), size = 10, arrowprops = dict(facecolor ='r', width=5, headlength = 0.1, headwidth=0.1))
plt.annotate("", (6.1,1400), (9.9,1400), size = 10, arrowprops = dict(facecolor ='r', width=5, headlength = 0.1, headwidth=0.1))
plt.annotate("", (10.1,1400), (16.9,1400), size = 10, arrowprops = dict(facecolor ='r', width=5, headlength = 0.1, headwidth=0.1))
plt.text(6.75, 1420, 'stratus phase', dict(size=20))
plt.text(11.75, 1420, 'convective phase', dict(size=20))


in_situ_file_MESONH_RADSURF_on = "BCBackground_XXL_RADSURF_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_RADSURF_on)
mesonh_RADSURF_aeroon = data["name"][1:]

in_situ_file_MESONH_RADSURF_off = "BCBackground_XXL_RADSURF_aerooff.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_RADSURF_off)
mesonh_RADSURF_aerooff = data["name"][1:]

list_diff_RADSURF = []
for iter_time in range(24):
    diff_RADSURF = mesonh_RADSURF_aeroon[iter_time] - mesonh_RADSURF_aerooff[iter_time]
    list_diff_RADSURF.append(diff_RADSURF)

fig, ax = plt.subplots(figsize=(20,12))
plt.plot(list_cloud_time,list_diff_RADSURF,color='blue',linewidth=4.0)
plt.xlabel("Time UTC (h)", fontsize = 20)
plt.xticks(fontsize=20)
plt.xlim(-1,24)
plt.ylabel("SWRADSURF Difference ($W.m^{-2}$)", fontsize = 20 )
plt.yticks(fontsize=20)
plt.ylim(-15,15)
plt.grid(which='major', linestyle='-', linewidth='1', color='k')
plt.grid(which='minor', linestyle='-', linewidth='0.25', color='k')
plt.minorticks_on()





