# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os as os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Patch

zlim = 2000
in_situ_time = ["23","00","01","02","03","04","05","06","07","08","09","10","11"]

path_to_data = "/home/delbeke/Documents/LIMA/SAVE_insitu/UPS_Frequent Radiosoundings_Save_03_07_2016/"
#list_dir = os.listdir(path_to_data)
in_situ_file_23h = "LA201607022258.nc"
in_situ_file_00h = "LA201607022358.nc"
in_situ_file_01h = "LA201607030100.nc"
in_situ_file_02h = "LA201607030159.nc"
in_situ_file_03h = "LA201607030300.nc"
in_situ_file_04h = "LA201607030400.nc"
in_situ_file_05h = "LA201607030459.nc"
in_situ_file_06h = "LA201607030602.nc"
in_situ_file_07h = "LA201607030658.nc"
in_situ_file_08h = "LA201607030759.nc"
in_situ_file_09h = "LA201607030858.nc"
in_situ_file_10h = "LA201607031003.nc"
in_situ_file_11h = "LA201607031102.nc"

list_in_situ_product = []
in_situ_products_23h = xr.open_dataset(path_to_data + in_situ_file_23h)
list_in_situ_product.append(in_situ_products_23h)
in_situ_products_00h = xr.open_dataset(path_to_data + in_situ_file_00h)
list_in_situ_product.append(in_situ_products_00h)
in_situ_products_01h = xr.open_dataset(path_to_data + in_situ_file_01h)
list_in_situ_product.append(in_situ_products_01h)
in_situ_products_02h = xr.open_dataset(path_to_data + in_situ_file_02h)
list_in_situ_product.append(in_situ_products_02h)
in_situ_products_03h = xr.open_dataset(path_to_data + in_situ_file_03h)
list_in_situ_product.append(in_situ_products_03h)
in_situ_products_04h = xr.open_dataset(path_to_data + in_situ_file_04h)
list_in_situ_product.append(in_situ_products_04h)
in_situ_products_05h = xr.open_dataset(path_to_data + in_situ_file_05h)
list_in_situ_product.append(in_situ_products_05h)
in_situ_products_06h = xr.open_dataset(path_to_data + in_situ_file_06h)
list_in_situ_product.append(in_situ_products_06h)
in_situ_products_07h = xr.open_dataset(path_to_data + in_situ_file_07h)
list_in_situ_product.append(in_situ_products_07h)
in_situ_products_08h = xr.open_dataset(path_to_data + in_situ_file_08h)
list_in_situ_product.append(in_situ_products_08h)
in_situ_products_09h = xr.open_dataset(path_to_data + in_situ_file_09h)
list_in_situ_product.append(in_situ_products_09h)
in_situ_products_10h = xr.open_dataset(path_to_data + in_situ_file_10h)
list_in_situ_product.append(in_situ_products_10h)
in_situ_products_11h = xr.open_dataset(path_to_data + in_situ_file_11h)
list_in_situ_product.append(in_situ_products_11h)



###################" Jet phase ####################
color_plot = ["blue", "green", "red"]
c=0
fig, axs = plt.subplots(1,3,figsize=(30, 30))
for iter_list in range(1,7,2):
    alt = list_in_situ_product[iter_list].ALT.values
    rh = list_in_situ_product[iter_list].RH.values
    axs[0].plot(rh, alt, linewidth=4.0, label = in_situ_time[iter_list] + "h", color = color_plot[c] )
    axs[0].tick_params(axis="x", labelsize=25) 
    axs[0].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[0].set_xlim(60, 100)
    axs[0].set_xlabel('Relative Humidity (%)', fontsize = 30)
    axs[0].xaxis.set_label_coords(0.5, -0.025)
    # plt.yticks(fontsize=20)
    axs[0].set_ylabel('Altitude asl (m)', fontsize = 30)
    axs[0].set_ylim(0,zlim)
    axs[0].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[0].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[0].minorticks_on()
    #axs[0].legend(fontsize=20)
    #axs[0].set_title('Relative Humidity (%)',fontsize=20, y=1.01)
    axs[0].text(62.0, 1950.0, '(a)', fontsize=40, fontweight = 'bold')
    
    temp = list_in_situ_product[iter_list].TEMP.values
    axs[1].plot(temp, alt, linewidth=4.0, label = in_situ_time[iter_list] + "h", color = color_plot[c] )
    axs[1].tick_params(axis="x", labelsize=25) 
    axs[1].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[1].set_xlim(17, 25)
    axs[1].set_xlabel('Air Temperature (°C)', fontsize = 30)
    axs[1].xaxis.set_label_coords(0.5, -0.025)
    # plt.yticks(fontsize=20)
    axs[1].set_ylim(0,zlim)
    axs[1].axes.yaxis.set_ticklabels([])
    axs[1].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[1].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[1].minorticks_on()
    #axs[1].legend(fontsize=20)
    #axs[1].set_title('Air Temperature (°C)',fontsize=20, y=1.01)
    axs[1].text(17.5, 1950.0, '(b)', fontsize=40, fontweight = 'bold')
    
    
    ws = list_in_situ_product[iter_list].WSPD.values
    axs[2].plot(ws, alt, linewidth=4.0, label = in_situ_time[iter_list] + "h", color = color_plot[c] )
    axs[2].tick_params(axis="x", labelsize=25) 
    axs[2].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[2].set_xlim(0, 10)
    axs[2].set_xlabel('Wind Speed (m.$s^{-1}$)', fontsize = 30)
    axs[2].xaxis.set_label_coords(0.5, -0.0215)
    # plt.yticks(fontsize=20)
    axs[2].set_ylim(0,zlim)
    axs[2].axes.yaxis.set_ticklabels([])
    axs[2].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[2].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[2].minorticks_on()
    #axs[2].legend(fontsize=20)
    #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    axs[2].text(0.5, 1950.0, '(c)', fontsize=40, fontweight = 'bold')
    
    c+=1
    
plt.legend(fontsize=30,bbox_to_anchor=(-2.57, 0.1, 0.5, 1.0),handles=[Patch(color='blue', label='0h'), Patch(color='green', label='2h'), Patch(color='red', label='4h')])
    
    
    
###################" Stratus phase ####################
color_plot = ["blue", "green", "red"]
c=0
fig, axs = plt.subplots(1,3,figsize=(30, 30))
for iter_list in range(7,12,2):
    alt = list_in_situ_product[iter_list].ALT.values
    rh = list_in_situ_product[iter_list].RH.values
    axs[0].plot(rh, alt, linewidth=4.0, label = in_situ_time[iter_list] + "h", color = color_plot[c] )
    axs[0].tick_params(axis="x", labelsize=25) 
    axs[0].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[0].set_xlim(60, 100)
    axs[0].set_xlabel('Relative Humidity (%)', fontsize = 30)
    axs[0].xaxis.set_label_coords(0.5, -0.025)
    # plt.yticks(fontsize=20)
    axs[0].set_ylabel('Altitude asl (m)', fontsize = 30)
    axs[0].set_ylim(0,zlim)
    axs[0].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[0].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[0].minorticks_on()
    # axs[0].legend(fontsize=20)
    #axs[0].set_title('Relative Humidity (%)',fontsize=20, y=1.01)
    axs[0].text(62.0, 1950.0, '(a)', fontsize=40, fontweight = 'bold')
    
    temp = list_in_situ_product[iter_list].TEMP.values
    axs[1].plot(temp, alt, linewidth=4.0, label = in_situ_time[iter_list] + "h", color = color_plot[c] )
    axs[1].tick_params(axis="x", labelsize=25) 
    axs[1].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[1].set_xlim(17, 25)
    axs[1].set_xlabel('Air Temperature (°C)', fontsize = 30)
    axs[1].xaxis.set_label_coords(0.5, -0.025)
    # plt.yticks(fontsize=20)
    axs[1].set_ylim(0,zlim)
    axs[1].axes.yaxis.set_ticklabels([])
    axs[1].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[1].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[1].minorticks_on()
    # axs[1].legend(fontsize=20)
    #axs[1].set_title('Air Temperature (°C)',fontsize=20, y=1.01)
    axs[1].text(17.5, 1950.0, '(b)', fontsize=40, fontweight = 'bold')
    
    ws = list_in_situ_product[iter_list].WSPD.values
    axs[2].plot(ws, alt, linewidth=4.0, label = in_situ_time[iter_list] + "h", color = color_plot[c] )
    axs[2].tick_params(axis="x", labelsize=25) 
    axs[2].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[2].set_xlim(0, 10)
    axs[2].set_xlabel('Wind Speed (m.$s^{-1}$)', fontsize = 30)
    axs[2].xaxis.set_label_coords(0.5, -0.0215)
    # plt.yticks(fontsize=20)
    axs[2].set_ylim(0,zlim)
    axs[2].axes.yaxis.set_ticklabels([])
    axs[2].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[2].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[2].minorticks_on()
    # axs[2].legend(fontsize=20)
    #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    axs[2].text(0.5, 1950.0, '(c)', fontsize=40, fontweight = 'bold')
    
    c+=1

plt.legend(fontsize=30,bbox_to_anchor=(-2.57, 0.1, 0.5, 1.0),handles=[Patch(color='blue', label='6h'), Patch(color='green', label='8h'), Patch(color='red', label='10h')])
   