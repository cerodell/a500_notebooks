#!/usr/bin/env python3
# %% [markdown]
#
# # This notebook......


# %%

# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 14:03:45 2019

@author: rodell
"""

from matplotlib import pyplot as plt
from cr500.utils import plt_set, constants
import numpy as np
import xarray as xr
import context


# %% [markdown]
#
# # Read out .zarr files for Doppler and Insitu measurements

# %%

filein = str(context.pro_data_dir)+str('/xr/var_ds.zarr')
var_ds = xr.open_zarr(filein)

filein = str(context.pro_data_dir)+str('/xr/doppler_ds.zarr')
doppler_ds = xr.open_zarr(filein)

save = str(context.pro_data_dir)+str('/Images/')

# %% [markdown]
#
# # Solve for Kinematic Sensible heat flux... add to var_ds

# %%

F_sfc = ((var_ds.H) / (((var_ds.P0 * 100) / (constants.R * var_ds.TA[:,-1])) * constants.Cp))
var_ds.update({'F_sfc':F_sfc})


# %% [markdown]
#
# # Time avereage wind speed profile
#
# ######## Plot the time avereage wind speed profile fromt he cabauw tower. Data for 2001 to 2019

# %%

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
fig.suptitle('Hourly Time Averaged Wind Speed Profile \n 2001 - 2019 ', fontsize= plt_set.title_size, fontweight="bold")
# Get array of Time Averaged Wind Speed
wsp = np.array(avg_hour.F)

for i in range(len(z_list)-1):
    ax.plot(datetime, wsp[:,i] , label = z_str[i][:-2])
    ax.set_xlabel("Hour UTC \n (Local - 1)", fontsize= plt_set.label)
    ax.set_ylabel("Wind Speed $(ms^-1)$", fontsize= plt_set.label)
    ax.tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
    ax.xaxis.grid(color='gray', linestyle='dashed')
    ax.yaxis.grid(color='gray', linestyle='dashed')
    ax.set_facecolor('lightgrey')
#    ax.legend(loc='best')
    ax.legend(loc='upper right', bbox_to_anchor=(.80,1.0), shadow=True, ncol=6, title='Height AGL (meters)')
   
fig.savefig(save + 'Avg_Wind_Speed_Profile')


# %% [markdown]
#
# # Wind speed profile at Hieght
#
# ######## Plot  wind speed profile at height for the cabauw tower. Data for 2001 to 2019

# %%
###############################################################################################
"""########################## Wind speed profile at Hieght ############################"""
###############################################################################################

### Group by hour
avg_hour = var_ds.groupby('time.hour').mean(dim='time')

datetime = np.array(avg_hour.hour)
z_list = np.array(avg_hour.z)
datestr = [str(i) for i in datetime]
datestr = datestr[0::4]

fig, ax = plt.subplots(1,1, figsize=(12,10))
fig.suptitle('Hourly Time Averaged Wind Speed Profile at Height \n 2001 - 2019 ', fontsize= plt_set.title_size, fontweight="bold")

for i in range(len(datestr)):
    ax.plot(wsp.T[:,i],z_list, label = datestr[i])
    ax.set_xlabel("Wind Speed $(ms^-1)$", fontsize= plt_set.label)
    ax.set_ylabel("Height Above Ground Level  \n (AGL)", fontsize= plt_set.label)
    ax.tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
    ax.xaxis.grid(color='gray', linestyle='dashed')
    ax.yaxis.grid(color='gray', linestyle='dashed')
    ax.set_facecolor('lightgrey')
    #    ax.legend(loc='best')
    ax.legend(loc='upper right', bbox_to_anchor=(.80,1.0), shadow=True, ncol=6, title='Height AGL (meters)')
       
fig.savefig(save + 'Time_Averaged_Wind_Speed_Profile_at_Height')

    


# %% [markdown]
#
# # Time avereage wind speed profile by season


# %%

###############################################################################################
"""##################### Time avereage wind speed profile by season #########################"""
###############################################################################################
#### Group by hour
avg_seasons = var_ds.groupby('time.season').groups

#z_list = np.array(avg_seasons.z)
#z_str = [str(i) for i in z_list]
#
#print(z_str[0][:-2])
#
#fig, ax = plt.subplots(2,2, figsize=(12,10))
#fig.suptitle('Hourly Time Averaged Wind Speed Profile \n 2001 - 2019 ', fontsize= plt_set.title_size, fontweight="bold")
##fig.autofmt_xdate()
#
##for i in range(len(z_list)-1):
#wsp = np.array(avg_seasons.F[0,4])
#ax[0,0].plot(datetime, wsp , label = z_str[i][:-2])
#ax[0,0].set_xlabel("Hour UTC \n (Local - 1)", fontsize= plt_set.label)
#ax[0,0].set_ylabel("Wind Speed $(ms^-1)$", fontsize= plt_set.label)
#ax[0,0].tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
#ax[0,0].xaxis.grid(color='gray', linestyle='dashed')
#ax[0,0].yaxis.grid(color='gray', linestyle='dashed')
#ax[0,0].set_facecolor('lightgrey')
##    ax.legend(loc='best')
#ax[0,0].legend(loc='upper right', bbox_to_anchor=(.80,1.0), shadow=True, ncol=6, title='Height AGL (meters)')
#   
##fig.savefig(save + 'Time_Averaged_Wind_Speed_Profile')



# %% [markdown]
#
# # Merge Doppler and and Surface/Tower Data
#
#

# %%
###############################################################################################
"""##################### Merge Doppler and and Surface/Tower Data  #########################"""
###############################################################################################
#
#var_2018 = var_ds.sel(time = slice('2018-05-01', '2018-09-01'))
#avg_hour_2018 = var_2018.groupby('time.hour').mean(dim='time')
#
#doppler_hour = doppler_ds.groupby('time.hour').mean(dim='time')
#
#xarray_list = [avg_hour_2018, doppler_hour]
#master_ds = xr.merge(xarray_list)
#master_ds = xr.combine_nested(xarray_list, 'z')

#z_list = np.array(doppler_hour.height_1st_interval)
#z_list = z_list[0,:]


# %% [markdown]
#
# # Apply conditional statements to var_ds to find a time of stable unstable and natural surface layer. 

# %%


temp_200m = np.array(var_ds.TA[:,0])
temp_2m = np.array(var_ds.TA[:,-1])

#inversion = var_ds.where(var_ds.P0 > 999, -1)


inversion = var_ds.where(var_ds.TA[:,0] > var_ds.TA[:,-1], drop=False)


temp_200m_inv = np.array(inversion.TA[:,0])
temp_2m_inv = np.array(inversion.TA[:,-1])


























