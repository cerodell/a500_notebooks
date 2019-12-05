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

from cr500.utils import plt_set, constants, wind_eq
from matplotlib import pyplot as plt
import seaborn as sns
import xarray as xr
import numpy as np
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
# # Solve for Virtual Temperature, Air Density, Kinematic Sensible heat flux, and the Obukhov length ... add all to var_ds

# #### Should add the derivation with Latex
# %%
##  Air Density
rho = var_ds.P0*1.e2/(constants.Rd*(var_ds.TA002 + 273.15))

## Fleagle and bussinger eq. 6.31
Eb = var_ds.H + 0.02*var_ds.LE

## Virtural temperature 
Tv = var_ds.TA002 + 273.15  + 0.61*var_ds.Q002*1.e-3

## add all to var_ds
var_ds.update({'rho':rho})
var_ds.update({'Eb':Eb})
var_ds.update({'Tv':Tv})

## the Obukhov length
L = - var_ds.Tv*constants.Cp*rho*var_ds.UST**3./(constants.k*constants.g*var_ds.Eb)
#good = np.abs(Eb) > 1


## add to var_ds
var_ds.update({'L':L})

### MY OLD METHOD MAY DELETE!!!!!!!
#F_sfc = ((var_ds.H) / (((var_ds.P0 * 100) / (constants.R * var_ds.TA[:,-1])) * constants.Cp))
#var_ds.update({'F_sfc':F_sfc})

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
avg_hour = var_ds.groupby('time.hour').mean(dim='time')

datetime = np.array(avg_hour.hour)
z_list = np.array(avg_hour.z)
z_str = [str(i) for i in z_list]

print(z_str[0][:-2])

fig, ax = plt.subplots(1,1, figsize=(12,10))
fig.suptitle('Hourly Time Averaged Temperaure Profile \n 2001 - 2019 ', fontsize= plt_set.title_size, fontweight="bold")
# Get array of Time Averaged Wind Speed
wsp = np.array(avg_hour.F)

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
# # Apply conditional statements to var_ds to find a time of stable, unstable and nutral surface layer. 


# %%

###############################################################################################
"""##################### Apply conditional statements to var_ds #########################"""
###############################################################################################
## Solve for the dT/dz
grad = np.gradient(var_ds.TA, var_ds.z, axis = 1)
grad_mean = np.mean(grad, axis =1)

##  Add to the DataArray
var_ds.update({'dTdz_mean':(('time'), grad_mean)})
var_ds.update({'dTdz':(('time','z'), grad)})

## Drop any and all rain events
meso_ds = var_ds.where((var_ds.RAIN < 0.0001), drop=False)

### np.where hate NaN values so use this to ignor some runtime error
np.warnings.filterwarnings('ignore')

## Apply conditional statements
stable_i = var_ds.where((meso_ds.dTdz_mean < 0), drop=False)

unstable = var_ds.where(meso_ds.dTdz_mean > 0.01 , drop=False)

nutral = var_ds.where((meso_ds.dTdz_mean > 0) & (meso_ds.dTdz_mean < 0.006), drop=False)


#uts = np.array(stable.UST)
#wsp = np.array(stable.F)
#l_ob = np.array(stable.L)

# %% [markdown]

# # Solve for wind speed at height in stable conditions using Log Linear Wind Equation

# ##### $$M(z) = ({u_*} / k) * [ln(z/z_o) + 6 * (z/L)]$$ 

# %% 

###############################################################################################
"""################## Stable Conditions using Log Linear Wind Equation #######################"""
###############################################################################################
## Make play data for height of loglin model
#height = np.arange(10,201,1)

m_z = [] 
## I hate this loop its slow and stupid
for i in range(len(z_list)):
    m_z_i = wind_eq.loglin_stable(stable_i.UST,z_list[i],stable_i.L)
    m_z.append(m_z_i)
m_z = np.stack(m_z)

## Add to Stable DataArray
stable_i.update({'mz_stable':(('time','z'), m_z.T)})

## Remove any odd negative values...some occur and I dont know why 
stable = stable_i.where((stable_i.mz_stable > 0) & (stable_i.mz_stable < 40), drop=False)
#stable.where(stable.mz_stable > 0, drop=False)



# %% 

###############################################################################################
"""####################### Plot Stable Conditions mean and do Stats  #########################"""
###############################################################################################

fig, ax = plt.subplots(1,1, figsize=(12,10))
fig.suptitle('Log Linear Wind Profile in Stable Surface Layer \n Average 2001 - 2018', fontsize= plt_set.title_size, fontweight="bold")

ax.scatter(stable.F.mean(dim='time'),z_list, color = 'red', marker='+', label = 'in-situ')
ax.plot(stable.mz_stable.mean(dim='time'),z_list, color = 'k', label = 'model')

ax.set_xlabel("Wind Speed $(ms^-1)$", fontsize= plt_set.label)
ax.set_ylabel("Height Above Ground Level  \n (AGL)", fontsize= plt_set.label)
ax.tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
ax.xaxis.grid(color='gray', linestyle='dashed')
ax.yaxis.grid(color='gray', linestyle='dashed')
ax.set_facecolor('lightgrey')
ax.legend(loc='best')
   
fig.savefig(save + 'Log_Linear_Wind_Profile_Stable_Surface_Layer')


# %%
#fig, ax = plt.subplots(1,1, figsize=(12,10))
#fig.suptitle('Log Wind Profile in Nutrual Surface Layer', fontsize= plt_set.title_size, fontweight="bold")
sns.jointplot(x= stable.F, y= stable.mz_stable, kind='hex', bins=15)

# %% [markdown]

# # Solve for wind speed at height in nutral conditions using Log  Wind Equation

# ##### $$M(z) = ({u_*} / k) * ln(z/z_o)$$ 
# %% 

###############################################################################################
"""################## Nutrual Conditions using Linear Wind Equation #######################"""
###############################################################################################

## Make play data for height of loglin model
#height = np.arange(10,201,1)
m_z = [] 

## I hate this loop its slow and stupid
for i in range(len(z_list)):
    m_z_i = wind_eq.logwind_neutral(nutral.UST,z_list[i])
    m_z.append(m_z_i)
m_z = np.stack(m_z)

## Add to Stable DataArray
nutral.update({'mz_neutral':(('time','z'), m_z.T)})

## Remove any odd negative values...some occur and I dont know why 
nutral = nutral.where((nutral.mz_neutral > 0), drop=False)

# %%

###############################################################################################
"""####################### Plot Neutral Conditions mean and do Stats  #########################"""
###############################################################################################
fig, ax = plt.subplots(1,1, figsize=(12,10))
fig.suptitle('Log Wind Profile in Nutrual Surface Layer', fontsize= plt_set.title_size, fontweight="bold")


ax.scatter(nutral.F.mean(dim='time'),z_list, color = 'red', marker='+', label = 'in-situ')
ax.plot(nutral.mz_neutral.mean(dim='time'),z_list, color = 'k', label = 'model')

ax.set_xlabel("Wind Speed $(ms^-1)$", fontsize= plt_set.label)
ax.set_ylabel("Height Above Ground Level  \n (AGL)", fontsize= plt_set.label)
ax.tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
ax.xaxis.grid(color='gray', linestyle='dashed')
ax.yaxis.grid(color='gray', linestyle='dashed')
ax.set_facecolor('lightgrey')
ax.legend(loc='best')
   
fig.savefig(save + 'Log_Wind_Profile_Nutrual_Surface_Layer')

# %%
#fig, ax = plt.subplots(1,1, figsize=(12,10))
#fig.suptitle('Log Wind Profile in Nutrual Surface Layer', fontsize= plt_set.title_size, fontweight="bold")
sns.jointplot(x= nutral.F, y= nutral.mz_neutral, kind='hex', bins=15)


# %% [markdown]
#
# # Apply conditional statements to var_ds to find a time of stable unstable surface layer. 

# ##### ADD EQ!!!!!!!!!!!!!!


# %% 

###############################################################################################
"""################## Unstable Conditions using Radix Wind Equation #######################"""
###############################################################################################


"""#######  THIS SECTION NEEDS WORK 
## Make play data for height of loglin model
#height = np.arange(10,201,1)
m_z = [] 

## I hate this loop its slow and stupid
for i in range(len(z_list)):
    m_z_i = wind_eq.RxL()
    m_z.append(m_z_i)
m_z = np.stack(m_z)

## Add to Stable DataArray
unstable.update({'mz_unstable':(('time','z'), m_z.T)})

## Remove any odd negative values...some occur and I dont know why 
unstable = unstable.where((unstable.mz_unstable > 0), drop=False)



# %% [markdown]
#
# # Apply conditional statements to var_ds to find  natural surface layer. 


# %% 



# %% 

#sns.jointplot(x= stable2.F[:,-2], y= stable2.m_z_login[:,-2], kind='kde')


"""






