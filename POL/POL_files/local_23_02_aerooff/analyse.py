import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

zlim=2000

middle = 121
for i in range(1,5):
    MESO_path = "./"
    if i <10:
        MESO_file1 = "LOCA1.1.AEROF.00"+str(i)+".nc"
        #MESO_file1 = "03_0"+str(i)+"h.nc"
    else:
        MESO_file1 = "LOCA1.1.AEROF.0"+str(i)+".nc"
        #MESO_file1 = "03_"+str(i)+"h.nc"

    MESO_products = xr.open_dataset(MESO_path + MESO_file1)
    rVT = MESO_products.RVT.values[0,:,:,:]
    tHT = MESO_products.THT.values[0,:,:,:]
    zHAT = MESO_products.ZHAT.values
    dtcur = MESO_products.DTCUR.values
    lon = MESO_products.nj.values
    lat = MESO_products.ni.values
    lev = MESO_products.level.values
    wt = MESO_products.WT.values[0,:,:,:]

    cloud = MESO_products.CLDFR.values[0,:,:,:]
    fig = plt.figure(figsize=(20,12))
    cmap = plt.get_cmap('rainbow')
    cs = plt.pcolormesh(lon, zHAT, cloud[:,middle,:], cmap=cmap,vmin=0,vmax=1)
    cbar = plt.colorbar(cs,pad = 0.01)
    cbar.set_label('CLDFR', rotation=270,fontsize=20,labelpad=20)
    cbar.ax.tick_params(labelsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(np.arange(0, max(zHAT), 50.0))
    plt.yticks(fontsize=20)
    plt.ylim(0,zlim)
    plt.title('Cloud_fraction_Day: ' + str(dtcur) ,fontsize=20, y=1.01)
    plt.savefig(MESO_file1+'_CLDFR.png',dpi=200)

    fig = plt.figure(figsize=(20,12))
    cmap = plt.get_cmap('rainbow')
    cs = plt.pcolormesh(lon, zHAT, wt[:,middle,:], cmap=cmap)
    cbar = plt.colorbar(cs,pad = 0.01)
    cbar.set_label('w m.s-1', rotation=270,fontsize=20,labelpad=20)
    cbar.ax.tick_params(labelsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(np.arange(0, max(zHAT), 50.0))
    plt.yticks(fontsize=20)
    plt.ylim(0,zlim)
    plt.title('w: ' + str(dtcur) ,fontsize=20, y=1.01)
    plt.savefig(MESO_file1+'_WT.png',dpi=200)


    fig = plt.figure(figsize=(20,12))
    cmap = plt.get_cmap('gist_ncar')
    cs = plt.pcolormesh(lon, zHAT, tHT[:,middle,:], cmap=cmap,vmin=290,vmax=320)
    cbar = plt.colorbar(cs,pad = 0.01)
    cbar.set_label('PT K', rotation=270,fontsize=20,labelpad=20)
    cbar.ax.tick_params(labelsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(np.arange(0, max(zHAT), 50.0))
    plt.yticks(fontsize=20)
    plt.ylim(0,zlim)
    plt.title('PT: ' + str(dtcur) ,fontsize=20, y=1.01)
    plt.savefig(MESO_file1+'_THT.png',dpi=200)

    fig = plt.figure(figsize=(20,12))
    cmap = plt.get_cmap('gist_ncar')
    cs = plt.pcolormesh(lon, zHAT, rVT[:,middle,:], cmap=cmap,vmin=0.010,vmax=0.020)
    cbar = plt.colorbar(cs,pad = 0.01)
    cbar.set_label('rVT kg.kg', rotation=270,fontsize=20,labelpad=20)
    cbar.ax.tick_params(labelsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(np.arange(0, max(zHAT), 50.0))
    plt.yticks(fontsize=20)
    plt.ylim(0,zlim)
    plt.title('rVT: ' + str(dtcur) ,fontsize=20, y=1.01) 
    plt.savefig(MESO_file1+'_RVT.png',dpi=200)

