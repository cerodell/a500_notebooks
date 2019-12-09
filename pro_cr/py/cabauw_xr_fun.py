# %% [markdown]
#
# # This notebook......


# %%

# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 14:03:45 2019

@author: rodell
"""
from cr500.utils import plt_set, constants, wind_eq
from matplotlib import pyplot as plt
import scipy.stats as stats
import seaborn as sns
import xarray as xr
import numpy as np
import context

## np.where hates NaN values so use this to ignor some runtime error
np.warnings.filterwarnings('ignore')

# %% [markdown]
#
# # Read out .zarr files for Doppler and Insitu measurements

# %%

filein = str(context.pro_data_dir)+str('/xr/var_ds.zarr')
var_ds_i = xr.open_zarr(filein)

filein = str(context.pro_data_dir)+str('/xr/doppler_ds.zarr')
doppler_ds = xr.open_zarr(filein)


save = str(context.pro_data_dir)+str('/Images/')



# %% [markdown]
#
# # Solve for Virtual Temperature, Air Density, Kinematic Sensible heat flux, and the Obukhov length ... add all to var_ds

# #### Should add the derivation with Latex
# %%


##  Air Density
rho = var_ds_i.P0*1.e2/(constants.Rd*(var_ds_i.TA002 + 273.15))

## Fleagle and bussinger eq. 6.31
Eb = var_ds_i.H + 0.02*var_ds_i.LE

## Virtural temperature 
Tv = var_ds_i.TA002 + 273.15  + 0.61*var_ds_i.Q002*1.e-3

## add all to var_ds
var_ds_i.update({'rho':rho})
var_ds_i.update({'Eb':Eb})
var_ds_i.update({'Tv':Tv})

## the Obukhov length
L = - var_ds_i.Tv*constants.Cp*rho*var_ds_i.UST**3./(constants.k*constants.g*var_ds_i.Eb)
#good = np.abs(Eb) > 1

## add to var_ds
var_ds_i.update({'L':L})

var_ds = var_ds_i.where((var_ds_i.L > -150) & (var_ds_i.L < 150), drop=False)
min_L, max_L = round(np.nanmin(var_ds.L),3), round(np.nanmax(var_ds.L),3)
print(f"Min (L) {min_L}  & Max (L) {max_L}")

## Make an array of tower heights..will be use a lot for ploting later...
z_list = np.array(var_ds.z)
z_str = [str(i) for i in z_list]




# %% [markdown]
#
# ## Plot stability as a function of day  averaged of the 18 years period from 2001-2018
# #### Stability as definded by $$/phi = (z/L)$$ 
# #### L is the Obukhov length and z is the heihgt


# %%
###############################################################################################
"""############################## Plot the Obukhov length  ##################################"""
###############################################################################################
## Apply conditional statements
phi = (var_ds.z[-1]/var_ds.L)
var_ds.update({'phi':(('time'), phi)})
var_ds_i = var_ds.where((var_ds.phi > -150) & (var_ds.phi < 150), drop=False)
min_phi, max_phi, avg_phi = round(np.nanmin(var_ds.phi),3), round(np.nanmax(var_ds.phi),3), round(np.nanmean(np.array(var_ds.phi)),3)
print(f"Min (z/L) {min_phi} hPa & Max (z/L) {max_phi}, Avg (z/L) {avg_phi}")

avg_day = var_ds_i.groupby('time.dayofyear').mean(dim='time')


fig, ax = plt.subplots(1,1, figsize=(18,10))
fig.suptitle('Obukhov Length Scale anual variability \n an average of 2001-2018', fontsize= plt_set.title_size, fontweight="bold")

#ax.plot(avg_day.dayofyear, avg_day.phi ,color = 'k')
ax.plot(var_ds_i.time, var_ds_i.phi ,color = 'k')


ax.set_xlabel("Day of Year", fontsize= plt_set.label)
ax.set_ylabel("Stability (z/L)", fontsize= plt_set.label)
ax.tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
ax.xaxis.grid(color='gray', linestyle='dashed')
ax.yaxis.grid(color='gray', linestyle='dashed')
ax.set_facecolor('lightgrey')
   
fig.savefig(save + 'stability_z_L')


# %%
print(f"Min (z/L) {min_phi} hPa & Max (z/L) {max_phi}, Avg (z/L) {avg_phi}")

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


# %%

###############################################################################################
"""########################## Time avereage wind speed profile ############################"""
###############################################################################################
### Group by hour
#avg_hour = var_ds.groupby('time.hour').mean(dim='time')
#
#datetime = np.array(avg_hour.hour)
#z_list = np.array(avg_hour.z)
#z_str = [str(i) for i in z_list]
#
#print(z_str[0][:-2])
#
#fig, ax = plt.subplots(1,1, figsize=(12,10))
#fig.suptitle('Hourly Time Averaged Temperaure Profile \n 2001 - 2019 ', fontsize= plt_set.title_size, fontweight="bold")
## Get array of Time Averaged Wind Speed
#wsp = np.array(avg_hour.F)

#for i in range(len(z_list)-1):
#    ax.plot(datetime, avg_hour.TA[:,i] , label = z_str[i][:-2])
#    ax.set_xlabel("Hour UTC \n (Local - 1)", fontsize= plt_set.label)
#    ax.set_ylabel("Temperature $(K)$", fontsize= plt_set.label)
#    ax.tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
#    ax.xaxis.grid(color='gray', linestyle='dashed')
#    ax.yaxis.grid(color='gray', linestyle='dashed')
#    ax.set_facecolor('lightgrey')
##    ax.legend(loc='best')
#    ax.legend(loc='upper right', bbox_to_anchor=(.55,1.0), shadow=True, ncol=6, title='Height AGL (meters)')
#  
#dat_height = np.meshgrid(datetime,z_list)
#ax.contourf(dat_height[0], avg_hour.TA, z_list)
#ax.set_xlabel("Hour UTC \n (Local - 1)", fontsize= plt_set.label)
#ax.set_ylabel("Temperature $(K)$", fontsize= plt_set.label)
#ax.tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
#ax.xaxis.grid(color='gray', linestyle='dashed')
#ax.yaxis.grid(color='gray', linestyle='dashed')
#ax.set_facecolor('lightgrey')
##    ax.legend(loc='best')
#ax.legend(loc='upper right', bbox_to_anchor=(.55,1.0), shadow=True, ncol=6, title='Height AGL (meters)')
#
#fig.savefig(save + 'Avg_Hourly_Temp_Profile')

#avg_hour.TA.plot(x = 'time', y = 'z')



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
fig.suptitle('Hourly Time Averaged Wind Speed Profile at Height \n 2001 - 2018 ', fontsize= plt_set.title_size, fontweight="bold")

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
#fig.suptitle('Hourly Time Averaged Wind Speed Profile \n 2001 - 2018 ', fontsize= plt_set.title_size, fontweight="bold")
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

