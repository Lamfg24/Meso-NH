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
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
import statsmodels.api as sm

in_situ_time = ["0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23"]

path_to_NRJ = "/home/delbeke/Documents/LIMA/SAVE_insitu/KIT_NRJ_FLUX/"
path_to_radiation = "/home/delbeke/Documents/LIMA/SAVE_insitu/KIT_radiation_1min/V1/"

in_situ_file_NRJ = "tower10Hz_gobe_inrab_20160703_000000_l3_v01_dacciwa.nc"
in_situ_file_radiation = "Save_KIT_EB1min_radiation_20160703.nc"

in_situ_products_NRJ = xr.open_dataset(path_to_NRJ + in_situ_file_NRJ)
in_situ_products_radiation = xr.open_dataset(path_to_radiation + in_situ_file_radiation)

h_surf = in_situ_products_NRJ.H.values
l_surf = in_situ_products_NRJ.LvE.values
time_surf = in_situ_products_NRJ.time.values
list_date_surf = []
for iter_surf in range(len(time_surf)):
    timestamp = pd.Timestamp(time_surf[iter_surf])
    #print(timestamp)
    hour = timestamp.hour
    minute = timestamp.minute
    minute_decimal = minute/60
    date_surf = hour +  minute_decimal
    list_date_surf.append(date_surf)
    #print(hour,minute,minute/60)
smoothed_h = sm.nonparametric.lowess(exog=list_date_surf, endog=h_surf, frac=0.13)
smoothed_l = sm.nonparametric.lowess(exog=list_date_surf, endog=l_surf, frac=0.13)


rad_short_down = in_situ_products_radiation.rsds.values
time_rad = in_situ_products_radiation.time.values
list_date_rad = []
for iter_rad in range(len(time_rad)):
    timestamp = pd.Timestamp(time_rad[iter_rad])
    #print(timestamp)
    hour = timestamp.hour
    minute = timestamp.minute
    minute_decimal = minute/60
    date_surf = hour +  minute_decimal
    list_date_rad.append(date_surf)
smoothed_rad = sm.nonparametric.lowess(exog=list_date_rad, endog=rad_short_down, frac=0.13)

fig, axs = plt.subplots(3,1,figsize=(30, 30))
axs[0].scatter(list_date_rad, rad_short_down, linewidth=2.0, label = "Observations", color = 'k' )
axs[0].plot(list_date_rad[0:-1], smoothed_rad[0:-1,1], linewidth=4.0, label = "Observations", color = 'k' )
axs[0].tick_params(axis="x", labelsize=25) 
axs[0].tick_params(axis="y", labelsize=25) 
# plt.xticks(fontsize=20)
axs[0].set_xlim(0, 24)
axs[0].set_xlabel('Time UTC (h)', fontsize = 30, y=0.1)
axs[0].set_ylabel('SWRADSURF (W.$m^{-1}$)', fontsize = 30)
#axs[0].xaxis.set_label_coords(0.5, -0.025)
# plt.yticks(fontsize=20)
axs[0].set_ylim(0,800)
axs[0].grid(which='major', linestyle='-', linewidth='1', color='k')
axs[0].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
axs[0].minorticks_on()


axs[1].scatter(list_date_surf, h_surf, linewidth=2.0, label = "Observations", color = 'k' )
axs[1].plot(list_date_surf[0:-1], smoothed_h[0:-1,1], linewidth=4.0, label = "Observations", color = 'k' )
axs[1].tick_params(axis="x", labelsize=25) 
axs[1].tick_params(axis="y", labelsize=25) 
# plt.xticks(fontsize=20)
axs[1].set_xlim(0, 24)
axs[1].set_xlabel('Time UTC (h)', fontsize = 30, y=0.1)
axs[1].set_ylabel('H (W.$m^{-1}$)', fontsize = 30)
#axs[0].xaxis.set_label_coords(0.5, -0.025)
# plt.yticks(fontsize=20)
axs[1].set_ylim(0,120)
axs[1].grid(which='major', linestyle='-', linewidth='1', color='k')
axs[1].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
axs[1].minorticks_on()

axs[2].scatter(list_date_surf, l_surf, linewidth=2.0, label = "Observations", color = 'k' )
axs[2].plot(list_date_surf[0:-1], smoothed_l[0:-1,1], linewidth=4.0, label = "Observations", color = 'k' )
axs[2].tick_params(axis="x", labelsize=25) 
axs[2].tick_params(axis="y", labelsize=25) 
# plt.xticks(fontsize=20)
axs[2].set_xlim(0, 24)
axs[2].set_xlabel('Time UTC (h)', fontsize = 30, y=0.1)
axs[2].set_ylabel('LE (W.$m^{-1}$)', fontsize = 30)
#axs[0].xaxis.set_label_coords(0.5, -0.025)
# plt.yticks(fontsize=20)
axs[2].set_ylim(0,350)
axs[2].grid(which='major', linestyle='-', linewidth='1', color='k')
axs[2].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
axs[2].minorticks_on()

    
    






