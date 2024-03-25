#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 13:18:12 2021

@author: dell
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt


MESO_path = "./"

SWF_NET=[]

for i in range(1,5):
    MESO_path = "./"
    MESO_path = MESO_path + "background_23_02_aerooff/"
    MESO_file = "BACK1.1.AEROF.00"+str(i)+"dgg.nc"
    print(MESO_path + MESO_file)
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    MRr = MESO_products.SWF_NET.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    mMRr = np.mean(np.mean(MRr,axis=1),axis=1)
    SWF_NET.append(mMRr)
    ZHAt = MESO_products.ZHAT.values
    np.savez("zhat",ZHAt)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "background_03_05_aerooff/"
    MESO_file = "BACK2.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    MRr = MESO_products.SWF_NET.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    mMRr = np.mean(np.mean(MRr,axis=1),axis=1)
    SWF_NET.append(mMRr)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "background_06_08_aerooff/"
    MESO_file = "BACK3.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    MRr = MESO_products.SWF_NET.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    mMRr = np.mean(np.mean(MRr,axis=1),axis=1)
    SWF_NET.append(mMRr)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "background_09_11_aerooff/"
    MESO_file = "BACK4.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    MRr = MESO_products.SWF_NET.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    mMRr = np.mean(np.mean(MRr,axis=1),axis=1)
    SWF_NET.append(mMRr)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "background_12_14_aerooff/"
    MESO_file = "BACK5.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    MRr = MESO_products.SWF_NET.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    mMRr = np.mean(np.mean(MRr,axis=1),axis=1)
    SWF_NET.append(mMRr)



for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "background_15_17_aerooff/"
    MESO_file = "BACK6.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    MRr = MESO_products.SWF_NET.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    mMRr = np.mean(np.mean(MRr,axis=1),axis=1)
    SWF_NET.append(mMRr)
   
for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "background_18_20_aerooff/"
    MESO_file = "BACK7.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    MRr = MESO_products.SWF_NET.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    mMRr = np.mean(np.mean(MRr,axis=1),axis=1)
    SWF_NET.append(mMRr)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "background_21_23_aerooff/"
    MESO_file = "BACK8.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    MRr = MESO_products.SWF_NET.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    mMRr = np.mean(np.mean(MRr,axis=1),axis=1)
    SWF_NET.append(mMRr)

	
np.savez("BCBackground_XXL_SWF_NET_aerooff",SWF_NET)





