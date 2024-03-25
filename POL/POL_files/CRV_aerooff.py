#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 13:18:12 2021

@author: dell
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

water_density = 1000

MESO_path = "./"

NCT=[]
LWC=[]
CRV=[]

for i in range(1,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_23_02_aerooff/"
    MESO_file = "LOCA1.1.AEROF.00"+str(i)+"dgg.nc"
    print(MESO_path + MESO_file)
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    NCt = MESO_products.NCT.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    LWc = MESO_products.LWC.values[0,:,3:-3,3:-3]
    Rvcd = np.power((3*LWc/1000)/(4*np.pi*water_density*NCt*(1E6)),1/3)*1000000
    CRv = np.ma.masked_equal(Rvcd, 0)
    mCRv = np.nanmean(np.nanmean(CRv,axis=1),axis=1)
    CRV.append(mCRv)
    ZHAt = MESO_products.ZHAT.values
    np.savez("zhat",ZHAt)

for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_03_05_aerooff/"
    MESO_file = "LOCA2.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    NCt = MESO_products.NCT.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    LWc = MESO_products.LWC.values[0,:,3:-3,3:-3]
    Rvcd = np.power((3*LWc/1000)/(4*np.pi*water_density*NCt*(1E6)),1/3)*1000000
    CRv = np.ma.masked_equal(Rvcd, 0)
    mCRv = np.nanmean(np.nanmean(CRv,axis=1),axis=1)
    CRV.append(mCRv)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_06_08_aerooff/"
    MESO_file = "LOCA3.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    NCt = MESO_products.NCT.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    LWc = MESO_products.LWC.values[0,:,3:-3,3:-3]
    Rvcd = np.power((3*LWc/1000)/(4*np.pi*water_density*NCt*(1E6)),1/3)*1000000
    CRv = np.ma.masked_equal(Rvcd, 0)
    mCRv = np.nanmean(np.nanmean(CRv,axis=1),axis=1)
    CRV.append(mCRv)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_09_11_aerooff/"
    MESO_file = "LOCA4.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    NCt = MESO_products.NCT.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    LWc = MESO_products.LWC.values[0,:,3:-3,3:-3]
    Rvcd = np.power((3*LWc/1000)/(4*np.pi*water_density*NCt*(1E6)),1/3)*1000000
    CRv = np.ma.masked_equal(Rvcd, 0)
    mCRv = np.nanmean(np.nanmean(CRv,axis=1),axis=1)
    CRV.append(mCRv)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_12_14_aerooff/"
    MESO_file = "LOCA5.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    NCt = MESO_products.NCT.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    LWc = MESO_products.LWC.values[0,:,3:-3,3:-3]
    Rvcd = np.power((3*LWc/1000)/(4*np.pi*water_density*NCt*(1E6)),1/3)*1000000
    CRv = np.ma.masked_equal(Rvcd, 0)
    mCRv = np.nanmean(np.nanmean(CRv,axis=1),axis=1)
    CRV.append(mCRv)



for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_15_17_aerooff/"
    MESO_file = "LOCA6.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    NCt = MESO_products.NCT.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    LWc = MESO_products.LWC.values[0,:,3:-3,3:-3]
    Rvcd = np.power((3*LWc/1000)/(4*np.pi*water_density*NCt*(1E6)),1/3)*1000000
    CRv = np.ma.masked_equal(Rvcd, 0)
    mCRv = np.nanmean(np.nanmean(CRv,axis=1),axis=1)
    CRV.append(mCRv)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_18_20_aerooff/"
    MESO_file = "LOCA7.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    NCt = MESO_products.NCT.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    LWc = MESO_products.LWC.values[0,:,3:-3,3:-3]
    Rvcd = np.power((3*LWc/1000)/(4*np.pi*water_density*NCt*(1E6)),1/3)*1000000
    CRv = np.ma.masked_equal(Rvcd, 0)
    mCRv = np.nanmean(np.nanmean(CRv,axis=1),axis=1)
    CRV.append(mCRv)


for i in range(2,5):
    MESO_path = "./"
    MESO_path = MESO_path + "local_21_23_aerooff/"
    MESO_file = "LOCA8.1.AEROF.00"+str(i)+"dgg.nc"
    MESO_products = xr.open_dataset(MESO_path + MESO_file)
    NCt = MESO_products.NCT.values[0,:,3:-3,3:-3] # mixing ratio for cloud (g/kg)
    LWc = MESO_products.LWC.values[0,:,3:-3,3:-3]
    Rvcd = np.power((3*LWc/1000)/(4*np.pi*water_density*NCt*(1E6)),1/3)*1000000
    CRv = np.ma.masked_equal(Rvcd, 0)
    mCRv = np.nanmean(np.nanmean(CRv,axis=1),axis=1)
    CRV.append(mCRv)

	
np.savez("Loca_CRV_aerooff",CRV)





