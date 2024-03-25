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

REHU=[]

for i in range(1,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_23_02_aeroon/"
    MESO_file = "LOCA1.1.AERON.00"+str(i)+"dgg.nc"
    print(MESO_path + MESO_file)
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    MRc = MESO_products.REHU.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    REHU.append(MRc)
    print(len(REHU))
    ZHAt = MESO_products.ZHAT.values
    np.savez("zhat",ZHAt)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_03_05_aeroon/"
    MESO_file = "LOCA2.1.AERON.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    MRc = MESO_products.REHU.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    REHU.append(MRc)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_06_08_aeroon/"
    MESO_file = "LOCA3.1.AERON.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    MRc = MESO_products.REHU.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    REHU.append(MRc)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_09_11_aeroon/"
    MESO_file = "LOCA4.1.AERON.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    MRc = MESO_products.REHU.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    REHU.append(MRc)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_12_14_aeroon/"
    MESO_file = "LOCA5.1.AERON.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    MRc = MESO_products.REHU.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    REHU.append(MRc)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_15_17_aeroon/"
    MESO_file = "LOCA6.1.AERON.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    MRc = MESO_products.REHU.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    REHU.append(MRc)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_18_20_aeroon/"
    MESO_file = "LOCA7.1.AERON.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    MRc = MESO_products.REHU.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    REHU.append(MRc)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_21_23_aeroon/"
    MESO_file = "LOCA8.1.AERON.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    MRc = MESO_products.REHU.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    REHU.append(MRc)

   	
np.savez("Local_REHU3D_aeroon",REHU)




