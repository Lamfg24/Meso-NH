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
int_time = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]

path_to_NRJ = "/home/delbeke/Documents/LIMA/SAVE_insitu/KIT_NRJ_FLUX/"
path_to_radiation = "/home/delbeke/Documents/LIMA/SAVE_insitu/KIT_radiation_1min/V1/"
path_to_MESONH = "/home/delbeke/Documents/LIMA/DACCIWA_XXL_LIMA_BACK_RESULT/NPZ_FILES/"

in_situ_file_NRJ = "tower10Hz_gobe_inrab_20160703_000000_l3_v01_dacciwa.nc"
in_situ_file_radiation = "Save_KIT_EB1min_radiation_20160703.nc"
in_situ_file_MESONH_radsurf = "BCBackground_XXL_RADSURF_aeroon.npz"
in_situ_file_MESONH_radsurf_std = "BCBackground_XXL_RADSURFstd_aeroon.npz"
in_situ_file_MESONH_H = "BCBackground_XXL_H_aeroon.npz"
in_situ_file_MESONH_H_std = "BCBackground_XXL_Hstd_aeroon.npz"
in_situ_file_MESONH_LE = "BCBackground_XXL_LE_aeroon.npz"
in_situ_file_MESONH_LE_std = "BCBackground_XXL_LEstd_aeroon.npz"

in_situ_products_NRJ = xr.open_dataset(path_to_NRJ + in_situ_file_NRJ)
in_situ_products_radiation = xr.open_dataset(path_to_radiation + in_situ_file_radiation)
data = np.load(path_to_MESONH + in_situ_file_MESONH_radsurf)
mesonh_radsurf = data["name"][1:]
data = np.load(path_to_MESONH + in_situ_file_MESONH_radsurf_std)
mesonh_radsurf_std = data["name"][1:]
data = np.load(path_to_MESONH + in_situ_file_MESONH_H)
mesonh_H = data["name"][1:]
data = np.load(path_to_MESONH + in_situ_file_MESONH_H_std)
mesonh_H_std = data["name"][1:]
data = np.load(path_to_MESONH + in_situ_file_MESONH_LE)
mesonh_LE = data["name"][1:]
data = np.load(path_to_MESONH + in_situ_file_MESONH_LE_std)
mesonh_LE_std = data["name"][1:]


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



fig, axs = plt.subplots(3,1,figsize=(20, 30))
axs[0].scatter(list_date_rad, rad_short_down, linewidth=2.0, label = "Observations", color = 'k' )
axs[0].plot(list_date_rad[0:-1], smoothed_rad[0:-1,1], linewidth=4.0, label = "Observations Fit", color = 'k' )
axs[0].plot(int_time, mesonh_radsurf, linewidth=4.0, label = "Model mean values", color = 'r' )
axs[0].fill_between(int_time, mesonh_radsurf + mesonh_radsurf_std, alpha=0.2, label = "Model standard deviation", color = 'r' )
axs[0].tick_params(axis="x", labelsize=25) 
axs[0].tick_params(axis="y", labelsize=25) 
# plt.xticks(fontsize=20)
axs[0].set_xlim(0, 24)
axs[0].set_xlabel('Time UTC (h)', fontsize = 30, y=0.1)
axs[0].set_ylabel('SWRADSURF (W.$m^{-2}$)', fontsize = 30)
#axs[0].xaxis.set_label_coords(0.5, -0.025)
# plt.yticks(fontsize=20)
axs[0].set_ylim(0,800)
axs[0].grid(which='major', linestyle='-', linewidth='1', color='k')
axs[0].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
axs[0].legend(fontsize=20)
axs[0].minorticks_on()
axs[0].text(-2.0, 825.0, '(a)', fontsize=25, fontweight = 'bold')


axs[1].scatter(list_date_surf, h_surf, linewidth=2.0, label = "Observations", color = 'k' )
axs[1].plot(list_date_surf[0:-1], smoothed_h[0:-1,1], linewidth=4.0, label = "Observations Fit", color = 'k' )
axs[1].plot(int_time, mesonh_H, linewidth=4.0, label = "Model mean values", color = 'r' )
axs[1].fill_between(int_time, mesonh_H + mesonh_H_std, alpha=0.2, label = "Model standard deviation", color = 'r' )
axs[1].tick_params(axis="x", labelsize=25) 
axs[1].tick_params(axis="y", labelsize=25) 
# plt.xticks(fontsize=20)
axs[1].set_xlim(0, 24)
axs[1].set_xlabel('Time UTC (h)', fontsize = 30, y=0.1)
axs[1].set_ylabel('H (W.$m^{-2}$)', fontsize = 30)
#axs[0].xaxis.set_label_coords(0.5, -0.025)
# plt.yticks(fontsize=20)
axs[1].set_ylim(0,120)
axs[1].grid(which='major', linestyle='-', linewidth='1', color='k')
axs[1].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
axs[1].legend(fontsize=20)
axs[1].minorticks_on()
axs[1].text(-2.0, 125.0, '(b)', fontsize=25, fontweight = 'bold')

axs[2].scatter(list_date_surf, l_surf, linewidth=2.0, label = "Observations", color = 'k' )
axs[2].plot(list_date_surf[0:-1], smoothed_l[0:-1,1], linewidth=4.0, label = "Observations Fit", color = 'k' )
axs[2].plot(int_time, mesonh_LE, linewidth=4.0, label = "Model mean values", color = 'r' )
axs[2].fill_between(int_time, mesonh_LE + mesonh_LE_std, alpha=0.2, label = "Model standard deviation", color = 'r' )
axs[2].tick_params(axis="x", labelsize=25) 
axs[2].tick_params(axis="y", labelsize=25) 
# plt.xticks(fontsize=20)
axs[2].set_xlim(0, 24)
axs[2].set_xlabel('Time UTC (h)', fontsize = 30, y=0.1)
axs[2].set_ylabel('LE (W.$m^{-2}$)', fontsize = 30)
#axs[0].xaxis.set_label_coords(0.5, -0.025)
# plt.yticks(fontsize=20)
axs[2].set_ylim(0,350)
axs[2].grid(which='major', linestyle='-', linewidth='1', color='k')
axs[2].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
axs[2].legend(fontsize=20)
axs[2].minorticks_on()
axs[2].text(-2.0, 365.0, '(c)', fontsize=25, fontweight = 'bold')
    
    






