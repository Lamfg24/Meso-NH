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
    MESO_path = MESO_path + "clean_23_02_aeroon/"
    MESO_file = "CLEA1.1.AERON.00"+str(i)+"dgg.nc"
    print(MESO_path + MESO_file)
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "clean_03_05_aeroon/"
    MESO_file = "CLEA2.1.AERON.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "clean_06_08_aeroon/"
    MESO_file = "CLEA3.1.AERON.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "clean_09_11_aeroon/"
    MESO_file = "CLEA4.1.AERON.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)



for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "clean_12_14_aeroon/"
    MESO_file = "CLEA5.1.AERON.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)



for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "clean_15_17_aeroon/"
    MESO_file = "CLEA6.1.AERON.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "clean_18_20_aeroon/"
    MESO_file = "CLEA7.1.AERON.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "clean_21_23_aeroon/"
    MESO_file = "CLEA8.1.AERON.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    TKET = MESO_products.TKET.values[0,:,3:-3,3:-3] # turbulent kinetic energy (m2/s2)
    mTKET = np.mean(np.mean(TKET,axis=1),axis=1)
    TKET_l.append(mTKET)


	
np.savez("Clean_TKET_aeroon",TKET_l)




