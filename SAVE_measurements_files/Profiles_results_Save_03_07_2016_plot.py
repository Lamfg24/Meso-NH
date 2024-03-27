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

path_to_MESONH = "/home/delbeke/Documents/LIMA/DACCIWA_XXL_LIMA_BACK_RESULT/NPZ_FILES/"

in_situ_file_MESONH_MRC = "BCBackground_XXL_MRC_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_MRC)
mesonh_mrc = data["name"][1:]
in_situ_file_MESONH_temp = "BCBackground_XXL_TEMP_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_temp)
mesonh_temp = data["name"][1:]
in_situ_file_MESONH_rh = "BCBackground_XXL_REHU_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_rh)
mesonh_rh = data["name"][1:]
in_situ_file_MESONH_equi_pot_temp = "BCBackground_XXL_THETAE_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_equi_pot_temp)
mesonh_equi_pot_temp = data["name"][1:]
in_situ_file_MESONH_u = "BCBackground_XXL_UT_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_u)
mesonh_u = data["name"][1:]
in_situ_file_MESONH_v = "BCBackground_XXL_VT_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_v)
mesonh_v = data["name"][1:]
in_situ_file_MESONH_nc =  "BCBackground_XXL_NCT_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_nc)
mesonh_nc = data["name"][1:]
in_situ_file_MESONH_na = "BCBackground_XXL_NAERO_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_na)
mesonh_na = data["name"][1:]
in_situ_file_MESONH_rc = "BCBackground_XXL_CRV_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_rc)
mesonh_rc = data["name"][1:]
in_situ_file_MESONH_tke ="BCBackground_XXL_TKET_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_tke)
mesonh_tke = data["name"][1:]
in_situ_file_MESONH_swhr = "BCBackground_XXL_DTRAD_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_swhr)
mesonh_swhr = data["name"][1:]
in_situ_file_MESONH_lwhr = "BCBackground_XXL_DTRAD_LW_aeroon.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_lwhr)
mesonh_lwhr = data["name"][1:]
in_situ_file_MESONH_alt = "zhat.npz"
data = np.load(path_to_MESONH + in_situ_file_MESONH_alt)
mesonh_alt = data["name"][:]
 
 
mesonh_ws = np.sqrt(mesonh_u*mesonh_u + mesonh_v*mesonh_v)


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

##################" Jet phase ####################
color_plot = ["blue", "green", "red"]
c=0
fig, axs = plt.subplots(1,10,figsize=(50, 30))
for iter_list in range(0,5,2):
    if iter_list == 0:
        base_cloud_plot = list_base_cloud[0]
        top_cloud_plot = list_top_cloud[0]
    if iter_list == 2:
        base_cloud_plot = list_base_cloud[2]
        top_cloud_plot = list_top_cloud[2]
    if iter_list == 4:
        base_cloud_plot = list_base_cloud[4]
        top_cloud_plot = list_top_cloud[4]
        
    
    alt = mesonh_alt[:]
    temp = mesonh_temp[iter_list,:]
    axs[0].plot(temp, alt, linewidth=4.0, color = color_plot[c] )
    axs[0].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[0].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[0].tick_params(axis="x", labelsize=25) 
    axs[0].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[0].set_xlim(16, 25)
    axs[0].set_xlabel('T (°C)', fontsize = 30)
    axs[0].xaxis.set_label_coords(0.5, -0.025)
    # plt.yticks(fontsize=20)
    axs[0].set_ylabel('Altitude (m)', fontsize = 30)
    axs[0].set_ylim(0,zlim)
    axs[0].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[0].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[0].minorticks_on()
    #axs[0].legend(fontsize=20)
    #axs[0].set_title('Relative Humidity (%)',fontsize=20, y=1.01)
    plt.text(-740.0, 1950.0, '(a)', fontsize=40, fontweight = 'bold')
    
    rh = mesonh_rh[iter_list,:]
    axs[1].plot(rh, alt, linewidth=4.0, color = color_plot[c] )
    axs[1].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[1].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[1].tick_params(axis="x", labelsize=25) 
    axs[1].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[1].set_xlim(55, 100)
    axs[1].set_xlabel('RH (%)', fontsize = 30)
    axs[1].xaxis.set_label_coords(0.5, -0.025)
    # plt.yticks(fontsize=20)
    axs[1].set_ylim(0,zlim)
    axs[1].axes.yaxis.set_ticklabels([])
    axs[1].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[1].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[1].minorticks_on()
    #axs[1].legend(fontsize=20)
    #axs[1].set_title('Air Temperature (°C)',fontsize=20, y=1.01)
    plt.text(-663.0, 1950.0, '(b)', fontsize=40, fontweight = 'bold')
    
    equi_pot_temp = mesonh_equi_pot_temp[iter_list,:]
    axs[2].plot(equi_pot_temp, alt, linewidth=4.0, color = color_plot[c] )
    axs[2].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[2].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[2].tick_params(axis="x", labelsize=25) 
    axs[2].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[2].set_xlim(330, 360)
    axs[2].set_xlabel('$\Theta_{e}$ (K)', fontsize = 30)
    axs[2].xaxis.set_label_coords(0.5, -0.025)
    # plt.yticks(fontsize=20)
    axs[2].set_ylim(0,zlim)
    axs[2].axes.yaxis.set_ticklabels([])
    axs[2].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[2].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[2].minorticks_on()
    #axs[2].legend(fontsize=20)
    #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-587.0, 1950.0, '(c)', fontsize=40, fontweight = 'bold')
    
    ws = mesonh_ws[iter_list,:]
    axs[3].plot(ws, alt, linewidth=4.0, color = color_plot[c] )
    axs[3].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[3].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[3].tick_params(axis="x", labelsize=25) 
    axs[3].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[3].set_xlim(0, 8)
    axs[3].set_xlabel('$w_{s}$ (m.$s^{-1}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[3].set_ylim(0,zlim)
    axs[3].axes.yaxis.set_ticklabels([])
    axs[3].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[3].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[3].minorticks_on()
    # #axs[].legend(fontsize=20)
    # #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-512.0, 1950.0, '(d)', fontsize=40, fontweight = 'bold')
    
    nc = mesonh_nc[iter_list,:]
    axs[4].plot(nc, alt, linewidth=4.0, color = color_plot[c] )
    axs[4].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[4].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[4].tick_params(axis="x", labelsize=25)
    axs[4].tick_params(axis="y", labelsize=25)
    # plt.xticks(fontsize=20)
    axs[4].set_xlim(0, 1000)
    axs[4].set_xlabel('$N_{c}$ (no.$cm^{-3}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[4].set_ylim(0,zlim)
    axs[4].axes.yaxis.set_ticklabels([])
    axs[4].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[4].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[4].minorticks_on()
    # #axs[].legend(fontsize=20)
    # #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-435.0, 1950.0, '(e)', fontsize=40, fontweight = 'bold')
    
    na = mesonh_na[iter_list,:]
    axs[5].plot(na, alt, linewidth=4.0, color = color_plot[c] )
    axs[5].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[5].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[5].tick_params(axis="x", labelsize=25) 
    axs[5].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    # axs[2].set_xlim(0, 10)
    axs[5].set_xlabel('$N_{a}$ (no.$cm^{-3}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[5].set_ylim(0,zlim)
    axs[5].axes.yaxis.set_ticklabels([])
    axs[5].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[5].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[5].minorticks_on()
    # #axs[].legend(fontsize=20)
    # #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-360.0, 1950.0, '(f)', fontsize=40, fontweight = 'bold')
    
    rc = mesonh_rc[iter_list,:]
    axs[6].plot(rc, alt, linewidth=4.0, color = color_plot[c] )
    axs[6].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[6].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[6].tick_params(axis="x", labelsize=25) 
    axs[6].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[6].set_xlim(0, 6)
    axs[6].set_xlabel('$r_{c}$ ($\mu$m)', fontsize = 30)
    axs[6].xaxis.set_label_coords(0.5, -0.023)
    # plt.yticks(fontsize=20)
    axs[6].set_ylim(0,zlim)
    axs[6].axes.yaxis.set_ticklabels([])
    axs[6].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[6].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[6].minorticks_on()
    # #axs[].legend(fontsize=20)
    # #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-284.0, 1950.0, '(g)', fontsize=40, fontweight = 'bold')
    
    tke = mesonh_tke[iter_list,:]
    axs[7].plot(tke, alt, linewidth=4.0, color = color_plot[c] )
    axs[7].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[7].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[7].tick_params(axis="x", labelsize=25) 
    axs[7].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[7].set_xlim(0, 0.25)
    axs[7].set_xlabel('TKE ($m^{2}.s^{-2}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[7].set_ylim(0,zlim)
    axs[7].axes.yaxis.set_ticklabels([])
    axs[7].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[7].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[7].minorticks_on()
    # #axs[].legend(fontsize=20)
    # #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-208.0, 1950.0, '(h)', fontsize=40, fontweight = 'bold')
    
    swhr = mesonh_swhr[iter_list,:]
    axs[8].plot(swhr, alt, linewidth=4.0, color = color_plot[c] )
    axs[8].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[8].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[8].tick_params(axis="x", labelsize=25) 
    axs[8].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    # axs[2].set_xlim(0, 10)
    axs[8].set_xlabel('SWHR (K.$d^{-1}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[8].set_ylim(0,zlim)
    axs[8].axes.yaxis.set_ticklabels([])
    axs[8].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[8].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[8].minorticks_on()
    plt.text(-130.0, 1950.0, '(i)', fontsize=40, fontweight = 'bold')
  
    
    lwhr = mesonh_lwhr[iter_list,:]
    axs[9].plot(lwhr, alt, linewidth=4.0, color = color_plot[c] )
    axs[9].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[9].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[9].tick_params(axis="x", labelsize=25) 
    axs[9].tick_params(axis="y", labelsize=25) 
    # # plt.xticks(fontsize=20)
    axs[9].set_xlabel('LWHR (K.$d^{-1}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[9].set_ylim(0,zlim)
    axs[9].axes.yaxis.set_ticklabels([])
    axs[9].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[9].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[9].minorticks_on()
    plt.text(-57.0, 1950.0, '(j)', fontsize=40, fontweight = 'bold')
    
    c+=1
    
plt.legend(fontsize=30,bbox_to_anchor=(-10.62, 0.105, 0.5, 1.0),handles=[Patch(color='blue', label='0h'), Patch(color='green', label='2h'), Patch(color='red', label='4h')])
plt.show() 
    
    
##################" stratus phase ####################
color_plot = ["blue", "green", "red"]
c=0
fig, axs = plt.subplots(1,10,figsize=(50, 30))
for iter_list in range(6,11,2):
    if iter_list == 6:
        base_cloud_plot = list_base_cloud[6]
        top_cloud_plot = list_top_cloud[6]
    if iter_list == 8:
        base_cloud_plot = list_base_cloud[8]
        top_cloud_plot = list_top_cloud[8]
    if iter_list == 10:
        base_cloud_plot = list_base_cloud[10]
        top_cloud_plot = list_top_cloud[10]
    
    alt = mesonh_alt[:]
    temp = mesonh_temp[iter_list,:]
    axs[0].plot(temp, alt, linewidth=4.0, color = color_plot[c] )
    axs[0].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[0].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[0].tick_params(axis="x", labelsize=25) 
    axs[0].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[0].set_xlim(15, 25)
    axs[0].set_xlabel('T (°C)', fontsize = 30)
    axs[0].xaxis.set_label_coords(0.5, -0.025)
    # plt.yticks(fontsize=20)
    axs[0].set_ylabel('Altitude (m)', fontsize = 30)
    axs[0].set_ylim(0,zlim)
    axs[0].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[0].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[0].minorticks_on()
    #axs[0].legend(fontsize=20)
    #axs[0].set_title('Relative Humidity (%)',fontsize=20, y=1.01)
    plt.text(-1155.0, 1950.0, '(a)', fontsize=40, fontweight = 'bold')
    
    rh = mesonh_rh[iter_list,:]
    axs[1].plot(rh, alt, linewidth=4.0, color = color_plot[c] )
    axs[1].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[1].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[1].tick_params(axis="x", labelsize=25) 
    axs[1].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[1].set_xlim(40, 100)
    axs[1].set_xlabel('RH (%)', fontsize = 30)
    axs[1].xaxis.set_label_coords(0.5, -0.025)
    # plt.yticks(fontsize=20)
    axs[1].set_ylim(0,zlim)
    axs[1].axes.yaxis.set_ticklabels([])
    axs[1].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[1].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[1].minorticks_on()
    #axs[1].legend(fontsize=20)
    #axs[1].set_title('Air Temperature (°C)',fontsize=20, y=1.01)
    plt.text(-1035.0, 1950.0, '(b)', fontsize=40, fontweight = 'bold')
    
    equi_pot_temp = mesonh_equi_pot_temp[iter_list,:]
    axs[2].plot(equi_pot_temp, alt, linewidth=4.0, color = color_plot[c] )
    axs[2].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[2].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[2].tick_params(axis="x", labelsize=25) 
    axs[2].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[2].set_xlim(330, 355)
    axs[2].set_xlabel('$\Theta_{e}$ (K)', fontsize = 30)
    axs[2].xaxis.set_label_coords(0.5, -0.025)
    # plt.yticks(fontsize=20)
    axs[2].set_ylim(0,zlim)
    axs[2].axes.yaxis.set_ticklabels([])
    axs[2].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[2].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[2].minorticks_on()
    #axs[2].legend(fontsize=20)
    #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-910.0, 1950.0, '(c)', fontsize=40, fontweight = 'bold')
    
    ws = mesonh_ws[iter_list,:]
    axs[3].plot(ws, alt, linewidth=4.0, color = color_plot[c] )
    axs[3].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[3].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[3].tick_params(axis="x", labelsize=25) 
    axs[3].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[3].set_xlim(0, 10)
    axs[3].set_xlabel('$w_{s}$ (m.$s^{-1}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[3].set_ylim(0,zlim)
    axs[3].axes.yaxis.set_ticklabels([])
    axs[3].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[3].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[3].minorticks_on()
    # #axs[].legend(fontsize=20)
    # #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-792.0, 1950.0, '(d)', fontsize=40, fontweight = 'bold')
    
    nc = mesonh_nc[iter_list,:]
    axs[4].plot(nc, alt, linewidth=4.0, color = color_plot[c] )
    axs[4].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[4].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[4].tick_params(axis="x", labelsize=25) 
    axs[4].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[4].set_xlim(0, 1200)
    axs[4].set_xlabel('$N_{c}$ (no.$cm^{-3}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[4].set_ylim(0,zlim)
    axs[4].axes.yaxis.set_ticklabels([])
    axs[4].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[4].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[4].minorticks_on()
    # #axs[].legend(fontsize=20)
    # #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-672.0, 1950.0, '(e)', fontsize=40, fontweight = 'bold')
    
    na = mesonh_na[iter_list,:]
    axs[5].plot(na, alt, linewidth=4.0, color = color_plot[c] )
    axs[5].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[5].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[5].tick_params(axis="x", labelsize=25) 
    axs[5].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    # axs[2].set_xlim(0, 10)
    axs[5].set_xlabel('$N_{a}$ (no.$cm^{-3}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[5].set_ylim(0,zlim)
    axs[5].axes.yaxis.set_ticklabels([])
    axs[5].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[5].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[5].minorticks_on()
    # #axs[].legend(fontsize=20)
    # #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-553.0, 1950.0, '(f)', fontsize=40, fontweight = 'bold')
    
    rc = mesonh_rc[iter_list,:]
    axs[6].plot(rc, alt, linewidth=4.0, color = color_plot[c] )
    axs[6].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[6].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[6].tick_params(axis="x", labelsize=25) 
    axs[6].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[6].set_xlim(0, 8)
    axs[6].set_xlabel('$r_{c}$ ($\mu$m)', fontsize = 30)
    axs[6].xaxis.set_label_coords(0.5, -0.023)
    # plt.yticks(fontsize=20)
    axs[6].set_ylim(0,zlim)
    axs[6].axes.yaxis.set_ticklabels([])
    axs[6].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[6].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[6].minorticks_on()
    # # #axs[].legend(fontsize=20)
    # # #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-435.0, 1950.0, '(g)', fontsize=40, fontweight = 'bold')
    
    tke = mesonh_tke[iter_list,:]
    axs[7].plot(tke, alt, linewidth=4.0, color = color_plot[c] )
    axs[7].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[7].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[7].tick_params(axis="x", labelsize=25) 
    axs[7].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[7].set_xlim(0, 0.18)
    axs[7].set_xlabel('TKE ($m^{2}.s^{-2}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[7].set_ylim(0,zlim)
    axs[7].axes.yaxis.set_ticklabels([])
    axs[7].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[7].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[7].minorticks_on()
    # # #axs[].legend(fontsize=20)
    # # #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-310.0, 1950.0, '(h)', fontsize=40, fontweight = 'bold')
    
    swhr = mesonh_swhr[iter_list,:]
    axs[8].plot(swhr, alt, linewidth=4.0, color = color_plot[c] )
    axs[8].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[8].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[8].tick_params(axis="x", labelsize=25) 
    axs[8].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[8].set_xlim(0, 50)
    axs[8].set_xlabel('SWHR (K.$d^{-1}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[8].set_ylim(0,zlim)
    axs[8].axes.yaxis.set_ticklabels([])
    axs[8].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[8].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[8].minorticks_on()
    plt.text(-190.0, 1950.0, '(i)', fontsize=40, fontweight = 'bold')
  
    
    lwhr = mesonh_lwhr[iter_list,:]
    axs[9].plot(lwhr, alt, linewidth=4.0, color = color_plot[c] )
    axs[9].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[9].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[9].tick_params(axis="x", labelsize=25) 
    axs[9].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[9].set_xlim(-80, 20)
    axs[9].set_xlabel('LWHR (K.$d^{-1}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[9].set_ylim(0,zlim)
    axs[9].axes.yaxis.set_ticklabels([])
    axs[9].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[9].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[9].minorticks_on()
    plt.text(-75.0, 1950.0, '(j)', fontsize=40, fontweight = 'bold')
    
    c+=1
    
plt.legend(fontsize=30,bbox_to_anchor=(-10.62, 0.105, 0.5, 1.0),handles=[Patch(color='blue', label='6h'), Patch(color='green', label='8h'), Patch(color='red', label='10h')])
    
  
    
# ###################" convection phase ####################
color_plot = ["blue", "green", "red","cyan"]
c=0
fig, axs = plt.subplots(1,10,figsize=(50, 30))
for iter_list in [12,14,16,17]:
    if iter_list == 12:
        base_cloud_plot = list_base_cloud[12]
        top_cloud_plot = list_top_cloud[12]
    if iter_list == 14:
        base_cloud_plot = list_base_cloud[14]
        top_cloud_plot = list_top_cloud[14]
    if iter_list == 16:
        base_cloud_plot = list_base_cloud[16]
        top_cloud_plot = list_top_cloud[16]
    if iter_list == 17:
        base_cloud_plot = list_base_cloud[17]
        top_cloud_plot = list_top_cloud[17]
        
    alt = mesonh_alt[:]
    temp = mesonh_temp[iter_list,:]
    axs[0].plot(temp, alt, linewidth=4.0, color = color_plot[c] )
    axs[0].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[0].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[0].tick_params(axis="x", labelsize=25) 
    axs[0].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[0].set_xlim(15, 30)
    axs[0].set_xlabel('T (°C)', fontsize = 30)
    axs[0].xaxis.set_label_coords(0.5, -0.025)
    # plt.yticks(fontsize=20)
    axs[0].set_ylabel('Altitude (m)', fontsize = 30)
    axs[0].set_ylim(0,zlim)
    axs[0].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[0].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[0].minorticks_on()
    #axs[0].legend(fontsize=20)
    #axs[0].set_title('Relative Humidity (%)',fontsize=20, y=1.01)
    plt.text(-918.0, 1950.0, '(a)', fontsize=40, fontweight = 'bold')
    
    rh = mesonh_rh[iter_list,:]
    axs[1].plot(rh, alt, linewidth=4.0, color = color_plot[c] )
    axs[1].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[1].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[1].tick_params(axis="x", labelsize=25) 
    axs[1].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[1].set_xlim(0, 100)
    axs[1].set_xlabel('RH (%)', fontsize = 30)
    axs[1].xaxis.set_label_coords(0.5, -0.025)
    # plt.yticks(fontsize=20)
    axs[1].set_ylim(0,zlim)
    axs[1].axes.yaxis.set_ticklabels([])
    axs[1].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[1].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[1].minorticks_on()
    #axs[1].legend(fontsize=20)
    #axs[1].set_title('Air Temperature (°C)',fontsize=20, y=1.01)
    plt.text(-825.0, 1950.0, '(b)', fontsize=40, fontweight = 'bold')
    
    equi_pot_temp = mesonh_equi_pot_temp[iter_list,:]
    axs[2].plot(equi_pot_temp, alt, linewidth=4.0, color = color_plot[c] )
    axs[2].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[2].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[2].tick_params(axis="x", labelsize=25) 
    axs[2].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    # axs[2].set_xlim(315, 360)
    axs[2].set_xlabel('$\Theta_{e}$ (K)', fontsize = 30)
    axs[2].xaxis.set_label_coords(0.5, -0.025)
    # plt.yticks(fontsize=20)
    axs[2].set_ylim(0,zlim)
    axs[2].axes.yaxis.set_ticklabels([])
    axs[2].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[2].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[2].minorticks_on()
    #axs[2].legend(fontsize=20)
    #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-725.0, 1950.0, '(c)', fontsize=40, fontweight = 'bold')
    
    ws = mesonh_ws[iter_list,:]
    axs[3].plot(ws, alt, linewidth=4.0, color = color_plot[c] )
    axs[3].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[3].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[3].tick_params(axis="x", labelsize=25) 
    axs[3].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[3].set_xlim(0, 10)
    axs[3].set_xlabel('$w_{s}$ (m.$s^{-1}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[3].set_ylim(0,zlim)
    axs[3].axes.yaxis.set_ticklabels([])
    axs[3].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[3].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[3].minorticks_on()
    #axs[].legend(fontsize=20)
    #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-632.0, 1950.0, '(d)', fontsize=40, fontweight = 'bold')
    
    nc = mesonh_nc[iter_list,:]
    axs[4].plot(nc, alt, linewidth=4.0, color = color_plot[c] )
    axs[4].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[4].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[4].tick_params(axis="x", labelsize=25) 
    axs[4].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[4].set_xlim(0, 800)
    axs[4].set_xlabel('$N_{c}$ (no.$cm^{-3}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[4].set_ylim(0,zlim)
    axs[4].axes.yaxis.set_ticklabels([])
    axs[4].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[4].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[4].minorticks_on()
    #axs[].legend(fontsize=20)
    #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-535.0, 1950.0, '(e)', fontsize=40, fontweight = 'bold')
    
    na = mesonh_na[iter_list,:]
    axs[5].plot(na, alt, linewidth=4.0, color = color_plot[c] )
    axs[5].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[5].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[5].tick_params(axis="x", labelsize=25) 
    axs[5].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    # axs[2].set_xlim(0, 10)
    axs[5].set_xlabel('$N_{a}$ (no.$cm^{-3}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[5].set_ylim(0,zlim)
    axs[5].axes.yaxis.set_ticklabels([])
    axs[5].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[5].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[5].minorticks_on()
    #axs[].legend(fontsize=20)
    #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-440.0, 1950.0, '(f)', fontsize=40, fontweight = 'bold')
    
    rc = mesonh_rc[iter_list,:]
    axs[6].plot(rc, alt, linewidth=4.0, color = color_plot[c] )
    axs[6].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[6].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[6].tick_params(axis="x", labelsize=25) 
    axs[6].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[6].set_xlim(0, 6)
    axs[6].set_xlabel('$r_{c}$ ($\mu$m)', fontsize = 30)
    axs[6].xaxis.set_label_coords(0.5, -0.023)
    # plt.yticks(fontsize=20)
    axs[6].set_ylim(0,zlim)
    axs[6].axes.yaxis.set_ticklabels([])
    axs[6].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[6].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[6].minorticks_on()
    #axs[].legend(fontsize=20)
    #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-345.0, 1950.0, '(g)', fontsize=40, fontweight = 'bold')
    
    tke = mesonh_tke[iter_list,:]
    axs[7].plot(tke, alt, linewidth=4.0, color = color_plot[c] )
    axs[7].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[7].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[7].tick_params(axis="x", labelsize=25) 
    axs[7].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[7].set_xlim(0, 0.3)
    axs[7].set_xlabel('TKE ($m^{2}.s^{-2}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[7].set_ylim(0,zlim)
    axs[7].axes.yaxis.set_ticklabels([])
    axs[7].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[7].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[7].minorticks_on()
    #axs[].legend(fontsize=20)
    #axs[2].set_title('Wind Speed (m.$s^{-1}$)',fontsize=20, y=1.01)
    plt.text(-247.0, 1950.0, '(h)', fontsize=40, fontweight = 'bold')
    
    swhr = mesonh_swhr[iter_list,:]
    axs[8].plot(swhr, alt, linewidth=4.0, color = color_plot[c] )
    axs[8].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[8].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[8].tick_params(axis="x", labelsize=25) 
    axs[8].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[8].set_xlim(0, 50)
    axs[8].set_xlabel('SWHR (K.$d^{-1}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[8].set_ylim(0,zlim)
    axs[8].axes.yaxis.set_ticklabels([])
    axs[8].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[8].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[8].minorticks_on()
    plt.text(-150.0, 1950.0, '(i)', fontsize=40, fontweight = 'bold')
  
    
    lwhr = mesonh_lwhr[iter_list,:]
    axs[9].plot(lwhr, alt, linewidth=4.0, color = color_plot[c] )
    axs[9].axhline(y=base_cloud_plot, color=color_plot[c], linestyle='-.',linewidth=3.0)
    axs[9].axhline(y=top_cloud_plot, color=color_plot[c], linestyle='--',linewidth=3.0)
    axs[9].tick_params(axis="x", labelsize=25) 
    axs[9].tick_params(axis="y", labelsize=25) 
    # plt.xticks(fontsize=20)
    axs[9].set_xlim(-60, 20)
    axs[9].set_xlabel('LWHR (K.$d^{-1}$)', fontsize = 30)
    # plt.yticks(fontsize=20)
    axs[9].set_ylim(0,zlim)
    axs[9].axes.yaxis.set_ticklabels([])
    axs[9].grid(which='major', linestyle='-', linewidth='1', color='k')
    axs[9].grid(which='minor', linestyle='-', linewidth='0.25', color='k')
    axs[9].minorticks_on()
    plt.text(-57.0, 1950.0, '(j)', fontsize=40, fontweight = 'bold')
    
    c+=1
    
plt.legend(fontsize=30,bbox_to_anchor=(-10.5, 0.13, 0.5, 1.0),handles=[Patch(color='blue', label='12h'), Patch(color='green', label='14h'), Patch(color='red', label='16h'), Patch(color='cyan', label='17h')])
    