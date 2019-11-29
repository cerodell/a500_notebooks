#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 14:03:45 2019

@author: rodell
"""

from matplotlib import pyplot as plt
from cr500.utils import plt_set
import numpy as np
import xarray as xr
import context




filein = str(context.pro_data_dir)+str('/xr/var_ds.zarr')
var_ds = xr.open_zarr(filein)



#var_ds.sel(time = slice('2001-01-02','2001-01-03'))

tmep_10m = var_ds.TA[:,0].mean()






###############################################################################################
"""########################## Time avereage wind speed profile ############################"""
###############################################################################################
### Group by hour
avg_hour = var_ds.groupby('time.hour').mean(dim='time')

datetime = np.array(avg_hour.hour)
z_list = np.array(avg_hour.z)
z_str = [str(i) for i in z_list]

print(z_str[0][:-2])

fig, ax = plt.subplots(1,1, figsize=(12,10))
fig.suptitle('Average Wind Speed Profile \n 2001 - 2019 ', fontsize= plt_set.title_size, fontweight="bold")
#fig.autofmt_xdate()

for i in range(len(z_list)-1):
    wsp = np.array(avg_hour.F[:,i])
    ax.plot(datetime, wsp , label = z_str[i][:-2])
    ax.set_xlabel("Time average hour UTC \n (Local - 1)", fontsize= plt_set.label)
    ax.set_ylabel("Wind Speed $(ms^-1)$", fontsize= plt_set.label)
    ax.tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
    ax.xaxis.grid(color='gray', linestyle='dashed')
    ax.yaxis.grid(color='gray', linestyle='dashed')
    ax.set_facecolor('lightgrey')
#    ax.legend(loc='best')
    ax.legend(loc='upper right', bbox_to_anchor=(.80,1.0), shadow=True, ncol=6, title='Height AGL (meters)')

    
    
    



###############################################################################################
"""########################## Time avereage wind speed profile ############################"""
###############################################################################################




fig, ax = plt.subplots(3,1, figsize=(16,10))
fig.suptitle('Air Quality Variables', fontsize=16, fontweight="bold")
fig.autofmt_xdate()

#ax[0].set_title("Particulate Matter 1 $\u03BCg/m^3$", fontsize = label)
#ax[0].set_xlabel("Datetime (PDT)", fontsize = label)
ax[0].set_ylabel("PM 1 ($\u03BCg/m^3$)", fontsize = label)
for i in range(len(pm1_f)):
    ax[0].plot(Date_f[i],pm1_f[i], alpha = 0.8, label = "Sensor  " + str_aq[i])
ax[0].legend(loc ='upper left')
ax[0].tick_params(axis='both', which='major', labelsize=tick_size)
ax[0].get_legend().remove()
ax[0].set_xticklabels([])
ax[0].xaxis.grid(color='gray', linestyle='dashed')
ax[0].yaxis.grid(color='gray', linestyle='dashed')
ax[0].set_ylim(0,65)
ax[0].set_facecolor('lightgrey')
xfmt = DateFormatter('%m/%d %H:%M')
ax[0].xaxis.set_major_formatter(xfmt)
ax[0].legend(loc='upper center', bbox_to_anchor=(.50,1.22), shadow=True, ncol=4)

###############################################################################################
"""########################## Time avereage wind speed profile ############################"""
###############################################################################################



inversion = avg_hour.where(avg_hour.TA[:,0] > avg_hour.TA[:,6], drop=True)
inversion.F[:,6].plot(size=6, x='hour')


