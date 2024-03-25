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
    MESO_path = MESO_path + "background_23_02_aerooff/"
    MESO_file = "DACC1.1.AEROF.00"+str(i)+"dgg.nc"
    print(MESO_path + MESO_file)
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "background_03_05_aerooff/"
    MESO_file = "DACC2.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "background_06_08_aerooff/"
    MESO_file = "DACC3.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "background_09_11_aerooff/"
    MESO_file = "DACC4.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)



for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "background_12_14_aerooff/"
    MESO_file = "DACC5.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)



for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "background_15_17_aerooff/"
    MESO_file = "DACC6.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "background_18_20_aerooff/"
    MESO_file = "DACC7.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "background_21_23_aerooff/"
    MESO_file = "DACC8.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)

	
np.savez("Background_TKET_aerooff",TKET_l)




