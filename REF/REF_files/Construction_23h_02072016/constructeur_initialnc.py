#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 13:44:31 2021

@author: dell
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from numpy import loadtxt



#MESO_path = '../001_prep_ideal_Chem_23h_background/'
MESO_path = '../modif/'
MESO_path = "../background_formation/"

Profil_path = './'
file_M = 'initial_file.nc'
#file_M = '15h.nc'

lines_w = loadtxt(Profil_path+"profil_subsi")

for i in range(1,10):
    print(i,"i")
    # print(i-1,"h")
    print("-----")
    
    lines = loadtxt(Profil_path+"profil_th/profil_th_0"+str(i))
    lines_rh = loadtxt(Profil_path+"profil_rv/profil_rv_0"+str(i))
    
    
    fh = netCDF4.Dataset(MESO_path+file_M, mode='r+')
    fh.variables['TENDTHFRC00'+str(i)][:]= lines
    fh.variables['TENDRVFRC00'+str(i)][:]= lines_rh
    fh.variables['WFRC00'+str(i)][:]= lines_w
    fh.close()
#print(lines)
#print(lines_rh)
    
for i in range(10,49):
    print(i,"i")
    # print(i-1,"h")
    print("-----")
    
    lines = loadtxt(Profil_path+"profil_th/profil_th_"+str(i))
    lines_rh = loadtxt(Profil_path+"profil_rv/profil_rv_"+str(i))
    
    
    fh = netCDF4.Dataset(MESO_path+file_M, mode='r+')
    fh.variables['TENDTHFRC0'+str(i)][:]= lines
    fh.variables['TENDRVFRC0'+str(i)][:]= lines_rh
    fh.variables['WFRC0'+str(i)][:]= lines_w
    fh.close()
#print(lines)
#print(lines_rh)


for i in range(49,50):
    print(i,"i")
    # print(i-1,"h")
    print("-----")

    lines = loadtxt(Profil_path+"profil_th/profil_th_"+str(47))
    lines_rh = loadtxt(Profil_path+"profil_rv/profil_rv_"+str(47))


    fh = netCDF4.Dataset(MESO_path+file_M, mode='r+')
    fh.variables['TENDTHFRC0'+str(i)][:]= lines
    fh.variables['TENDRVFRC0'+str(i)][:]= lines_rh
    fh.variables['WFRC0'+str(i)][:]= lines_w
    fh.close()
#print(lines)
#print(lines_rh)


for i in range(50,58):
    print(i,"i")
    print(i-48,"i-48")
    print("-----")

    lines = loadtxt(Profil_path+"profil_th/profil_th_0"+str(i-48))
    lines_rh = loadtxt(Profil_path+"profil_rv/profil_rv_0"+str(i-48))


    fh = netCDF4.Dataset(MESO_path+file_M, mode='r+')
    fh.variables['TENDTHFRC0'+str(i)][:]= lines
    fh.variables['TENDRVFRC0'+str(i)][:]= lines_rh
    fh.variables['WFRC0'+str(i)][:]= lines_w
    fh.close()
#print(lines)
#print(lines_rh)

for i in range(58,97):
    print(i,"i")
    print(i-48,"i-48")
    print("-----")

    lines = loadtxt(Profil_path+"profil_th/profil_th_"+str(i-48))
    lines_rh = loadtxt(Profil_path+"profil_rv/profil_rv_"+str(i-48))


    fh = netCDF4.Dataset(MESO_path+file_M, mode='r+')
    fh.variables['TENDTHFRC0'+str(i)][:]= lines
    fh.variables['TENDRVFRC0'+str(i)][:]= lines_rh
    fh.variables['WFRC0'+str(i)][:]= lines_w
    fh.close()


