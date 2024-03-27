#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 17:57:28 2024

@author: delbeke
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
import matplotlib.colors as mcolors

path_to_MESONH = "/home/delbeke/Documents/LIMA/DACCIWA_XXL_LIMA_BACK_RESULT/NPZ_FILES/"

in_situ_file_MESONH_MRC3D = "BCBackground_MRC3D_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_MRC3D)
mesonh_mrc3d = data["name"][1:]

list_integrated_mrc = []
for iter_mrc in range(24):
    one_mesonh_mrc3d = mesonh_mrc3d[iter_mrc]
    integrated_mrc_along_altitude = np.trapz(one_mesonh_mrc3d, axis=0)
    list_integrated_mrc.append(integrated_mrc_along_altitude)
    
length_x = 240 * 40
length_y = 240 * 40
length_x_array = np.linspace(0, length_x, 240)
length_y_array = np.linspace(0, length_y, 240)


fig, axs = plt.subplots(3,1,figsize=(10, 30))
cmap = plt.get_cmap('Greys_r')

cs = axs[0].pcolormesh(length_x_array, length_y_array, list_integrated_mrc[6], cmap=cmap, vmin=0, vmax=50)
cbar = plt.colorbar(cs,pad = 0.01)
cbar.set_label('LWP (g.$kg^{-1}$.m)', rotation=270,fontsize=20,labelpad=20)
cbar.ax.tick_params(labelsize=20)
axs[0].tick_params(axis="x", labelsize=25) 
axs[0].tick_params(axis="y", labelsize=25) 
axs[0].set_xlabel('x-axis length (m)', fontsize = 30, y=0.1)
axs[0].set_ylabel('y-axis length (m)', fontsize = 30)
axs[0].text(5.0, 9100.0, '(a)', fontsize=25, fontweight = 'bold', color = "bisque")

cs = axs[1].pcolormesh(length_x_array, length_y_array, list_integrated_mrc[12], cmap=cmap, vmin=0, vmax=50)
cbar = plt.colorbar(cs,pad = 0.01)
cbar.set_label('LWP (g.$kg^{-1}$.m)', rotation=270,fontsize=20,labelpad=20)
cbar.ax.tick_params(labelsize=20)
axs[1].tick_params(axis="x", labelsize=25) 
axs[1].tick_params(axis="y", labelsize=25) 
axs[1].set_xlabel('x-axis length (m)', fontsize = 30, y=0.1)
axs[1].set_ylabel('y-axis length (m)', fontsize = 30)
axs[1].text(5.0, 9100.0, '(c)', fontsize=25, fontweight = 'bold', color = "bisque")

cs = axs[2].pcolormesh(length_x_array, length_y_array, list_integrated_mrc[16], cmap=cmap, vmin=0, vmax=50)
cbar = plt.colorbar(cs,pad = 0.01)
cbar.set_label('LWP (g.$kg^{-1}$.m)', rotation=270,fontsize=20,labelpad=20)
cbar.ax.tick_params(labelsize=20)
axs[2].tick_params(axis="x", labelsize=25) 
axs[2].tick_params(axis="y", labelsize=25) 
axs[2].set_xlabel('x-axis length (m)', fontsize = 30, y=0.1)
axs[2].set_ylabel('y-axis length (m)', fontsize = 30)
axs[2].text(5.0, 9100.0, '(e)', fontsize=25, fontweight = 'bold', color = "bisque")
# plt.xticks(fontsize=20)
# plt.yticks(np.arange(0, max(zHAT), 50.0))
# plt.yticks(fontsize=20)
# plt.ylim(0,zlim)
# plt.title('Cloud_fraction_Day: ' + str(dtcur) ,fontsize=20, y=1.01)
# plt.savefig(MESO_file1+'_CLDFR.png',dpi=200)