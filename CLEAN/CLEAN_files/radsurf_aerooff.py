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

RADSURF=[]

for i in range(1,5):
    MESO_path = "./"
    MESO_path = MESO_path + "clean_23_02_aerooff/"
    MESO_file = "CLEA1.1.AEROF.00"+str(i)+"dgg.nc"
    print(MESO_path + MESO_file)
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    RADSWD_VIS = MESO_products.RADSWD_VIS.values[0,3:-3,3:-3] # net shortwave flux (W.m-2)
    RADSWD_NIR = MESO_products.RADSWD_NIR.values[0,3:-3,3:-3]
    RADLWD = MESO_products.RADLWD.values[0,3:-3,3:-3]
    RAD = RADSWD_VIS + RADSWD_NIR
    mRAD = np.mean(RAD)
    RADSURF.append(mRAD)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "clean_03_05_aerooff/"
    MESO_file = "CLEA2.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    RADSWD_VIS = MESO_products.RADSWD_VIS.values[0,3:-3,3:-3] # net shortwave flux (W.m-2)
    RADSWD_NIR = MESO_products.RADSWD_NIR.values[0,3:-3,3:-3]
    RADLWD = MESO_products.RADLWD.values[0,3:-3,3:-3]
    RAD = RADSWD_VIS + RADSWD_NIR
    mRAD = np.mean(RAD)
    RADSURF.append(mRAD)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "clean_06_08_aerooff/"
    MESO_file = "CLEA3.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    RADSWD_VIS = MESO_products.RADSWD_VIS.values[0,3:-3,3:-3] # net shortwave flux (W.m-2)
    RADSWD_NIR = MESO_products.RADSWD_NIR.values[0,3:-3,3:-3]
    RADLWD = MESO_products.RADLWD.values[0,3:-3,3:-3]
    RAD = RADSWD_VIS + RADSWD_NIR
    mRAD = np.mean(RAD)
    RADSURF.append(mRAD)



for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "clean_09_11_aerooff/"
    MESO_file = "CLEA4.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    RADSWD_VIS = MESO_products.RADSWD_VIS.values[0,3:-3,3:-3] # net shortwave flux (W.m-2)
    RADSWD_NIR = MESO_products.RADSWD_NIR.values[0,3:-3,3:-3]
    RADLWD = MESO_products.RADLWD.values[0,3:-3,3:-3]
    RAD = RADSWD_VIS + RADSWD_NIR
    mRAD = np.mean(RAD)
    RADSURF.append(mRAD)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "clean_12_14_aerooff/"
    MESO_file = "CLEA5.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    RADSWD_VIS = MESO_products.RADSWD_VIS.values[0,3:-3,3:-3] # net shortwave flux (W.m-2)
    RADSWD_NIR = MESO_products.RADSWD_NIR.values[0,3:-3,3:-3]
    RADLWD = MESO_products.RADLWD.values[0,3:-3,3:-3]
    RAD = RADSWD_VIS + RADSWD_NIR
    mRAD = np.mean(RAD)
    RADSURF.append(mRAD)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "clean_15_17_aerooff/"
    MESO_file = "CLEA6.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    RADSWD_VIS = MESO_products.RADSWD_VIS.values[0,3:-3,3:-3] # net shortwave flux (W.m-2)
    RADSWD_NIR = MESO_products.RADSWD_NIR.values[0,3:-3,3:-3]
    RADLWD = MESO_products.RADLWD.values[0,3:-3,3:-3]
    RAD = RADSWD_VIS + RADSWD_NIR
    mRAD = np.mean(RAD)
    RADSURF.append(mRAD)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "clean_18_20_aerooff/"
    MESO_file = "CLEA7.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    RADSWD_VIS = MESO_products.RADSWD_VIS.values[0,3:-3,3:-3] # net shortwave flux (W.m-2)
    RADSWD_NIR = MESO_products.RADSWD_NIR.values[0,3:-3,3:-3]
    RADLWD = MESO_products.RADLWD.values[0,3:-3,3:-3]
    RAD = RADSWD_VIS + RADSWD_NIR
    mRAD = np.mean(RAD)
    RADSURF.append(mRAD)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "clean_21_23_aerooff/"
    MESO_file = "CLEA8.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    RADSWD_VIS = MESO_products.RADSWD_VIS.values[0,3:-3,3:-3] # net shortwave flux (W.m-2)
    RADSWD_NIR = MESO_products.RADSWD_NIR.values[0,3:-3,3:-3]
    RADLWD = MESO_products.RADLWD.values[0,3:-3,3:-3]
    RAD = RADSWD_VIS + RADSWD_NIR
    mRAD = np.mean(RAD)
    RADSURF.append(mRAD)

	
np.savez("Clean_RADSURF_aerooff",RADSURF)



