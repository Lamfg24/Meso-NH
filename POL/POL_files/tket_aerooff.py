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

TKET_l = []

for i in range(1,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_23_02_aerooff/"
    MESO_file = "LOCA1.1.AEROF.00"+str(i)+"dgg.nc"
    print(MESO_path + MESO_file)
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_03_05_aerooff/"
    MESO_file = "LOCA2.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_06_08_aerooff/"
    MESO_file = "LOCA3.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_09_11_aerooff/"
    MESO_file = "LOCA4.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)



for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_12_14_aerooff/"
    MESO_file = "LOCA5.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)



for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_15_17_aerooff/"
    MESO_file = "LOCA6.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_18_20_aerooff/"
    MESO_file = "LOCA7.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_21_23_aerooff/"
    MESO_file = "LOCA8.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)

	
np.savez("Local_TKET_aerooff",TKET_l)




