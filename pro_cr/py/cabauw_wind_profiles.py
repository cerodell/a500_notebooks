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

var_ds = var_ds_i.where((var_ds_i.L > -1500) & (var_ds_i.L < 1500), drop=False)
min_L, max_L = round(np.nanmin(var_ds.L),3), round(np.nanmax(var_ds.L),3)
print(f"Min (L) {min_L}  & Max (L) {max_L}")

## Make an array of tower heights..will be use a lot for ploting later...
z_list = np.array(var_ds.z)
z_str = [str(i) for i in z_list]



# %% [markdown]
#
# # Merge Doppler and and Surface/Tower Data
# # TRY AND GET DOPPLER DATA WORKING WITH OTHER DATA ITS KILLING ME 
#

# %%
###############################################################################################
"""##################### Merge Doppler and and Surface/Tower Data  #########################"""
###############################################################################################

#rounded = np.array(doppler_ds.time, dtype='datetime64[h]')
#doppler_ds.update({'time':rounded})

#doppler_ds.time.sort_index()

#dtest = doppler_ds.time.dt.floor('min')
#dtest.resample(time='6H'

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

## Solve for the dT/dz...will be used later for the stability classes 
grad_dtdz = np.gradient(var_ds.TA, var_ds.z, axis = 1)
grad_dtdz_mean = np.mean(grad_dtdz, axis =1)

##  Add to the DataArray
var_ds.update({'dTdz_mean':(('time'), grad_dtdz_mean)})
var_ds.update({'dTdz':(('time','z'), grad_dtdz)})

## Drop any and all rain events
#meso_ds_i = var_ds.where((var_ds.RAIN < 0.0001), drop=False)
min_u_str, max_u_str = np.nanmin(var_ds.UST), np.nanmax(var_ds.UST)
print(f"Min UST {min_u_str} ms^-1 & Max UST {max_u_str} ms^-1")

## Drop any and all rain events
meso_ds_i = var_ds.where((var_ds.RAIN < 0.0001), drop=False)
min_rain, max_rain = np.nanmin(meso_ds_i.RAIN), np.nanmax(meso_ds_i.RAIN)
print(f"Min Rain {min_rain} cm & Max Rain {max_rain} cm")

## Drop low press events....should likely expand on this of press gradient over time.....
meso_ds = meso_ds_i.where((meso_ds_i.P0 > 990. ), drop=False)
min_press, max_press = round(np.nanmin(meso_ds.P0),3), round(np.nanmax(meso_ds.P0),3)
print(f"Min Press {min_press} hPa & Max Press {max_press} hPa")


# %% [markdown]

# ### Filtter for dirrefnt stability classes by comapring dT/dz to the dry adiabtic laps rate

# %%

## Apply conditional statements
stable_i = var_ds.where((meso_ds.dTdz_mean < 0), drop=False)
#
nutral_i = var_ds.where((meso_ds.dTdz_mean > 0) & (meso_ds.dTdz_mean < 0.003), drop=False)
##nutral_i = var_ds.where(meso_ds.dTdz_mean == 0, drop=False)
#
unstable_i = var_ds.where(meso_ds.dTdz_mean > 0.01 , drop=False)

## Solve for the dF/dz...will be use to see where the top tow heigfht have
## nearly zero change....this is need for the radix equation 
grad_dfdz = np.gradient(unstable_i.F, unstable_i.z, axis = 1)
grad_dfdz_mean = np.mean(grad_dfdz, axis =1)

##  Add to the DataArray
unstable_i.update({'dFdz_mean':(('time'), grad_dfdz_mean)})
unstable_i.update({'dFdz':(('time','z'), grad_dfdz)})

unstable_ii = unstable_i.where((unstable_i.dFdz[:,0] > -0.001) & (unstable_i.dFdz[:,0] < 0.001))


"""
I have been tesing this method and it yeild little differncec as comapred 
to my other arropch to define stability

## Apply conditional statements
#phi = (meso_ds.z[-1]/meso_ds.L)
#meso_ds.update({'phi':(('time'), phi)})
#min_phi, max_phi = round(np.nanmin(meso_ds.L),3), round(np.nanmax(meso_ds.L),3)
#print(f"Min (z/L) {min_phi} hPa & Max (z/L) {max_phi}")
#
#stable_i = var_ds.where((meso_ds.phi > 0.), drop=False)
#
#nutral_i = var_ds.where(meso_ds.phi == 0, drop=False)
#
#unstable_i = var_ds.where(meso_ds.phi < 0. , drop=False)
"""


# %% [markdown]

# # Show the Distribution of dT/dz for each stability conditon


# %% 

f, axes = plt.subplots(1, 3, figsize=(16, 8), sharex=True)
sns.kdeplot(stable_i.dTdz_mean, shade=True, color="r", ax=axes[0], label = 'Stable')
sns.kdeplot(nutral_i.dTdz_mean, shade=True, color="b", ax=axes[1], label = 'Nutral')
sns.kdeplot(unstable_i.dTdz_mean, shade=True, color="g", ax=axes[2], label = 'Unstable')
sns.set_style("darkgrid", {"axes.facecolor": ".9"})



# %% [markdown]

# # Solve for wind speed at height in stable conditions using Log Linear Wind Equation

# ### $$M(z) = ({u_*} / k) * [ln(z/z_o) + 6 * (z/L)]$$ 

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

sum_me = np.array(stable.F[:,-4])
stable_non_nans = (~np.isnan(sum_me)).sum()
print(f'We have {stable_non_nans} number of times of a Stable BL')


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
## Hexbin plot or 2D histogram
j = sns.jointplot(x= stable.F, y= stable.mz_stable, kind='hex', bins=15)
j.annotate(stats.pearsonr)


# %% [markdown]

# # Solve for wind speed at height in nutral conditions using Log  Wind Equation

# ### $$M(z) = ({u_*} / k) * ln(z/z_o)$$ 
# %% 

###############################################################################################
"""################## Nutrual Conditions using Linear Wind Equation #######################"""
###############################################################################################

## Make play data for height of loglin model
#height = np.arange(10,201,1)
m_z = [] 

## I hate this loop its slow and stupid
for i in range(len(z_list)):
    m_z_i = wind_eq.logwind_neutral(nutral_i.UST,z_list[i])
    m_z.append(m_z_i)
m_z = np.stack(m_z)

## Add to Stable DataArray
nutral_i.update({'mz_neutral':(('time','z'), m_z.T)})

## Remove any odd negative values...some occur and I dont know why 
nutral = nutral_i.where((nutral_i.mz_neutral > 0), drop=False)

sum_me = np.array(nutral.F[:,-4])
nutral_non_nans = (~np.isnan(sum_me)).sum()
print(f'We have {nutral_non_nans} number of times of a Nutral BL')


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
## Hexbin plot or 2D histogram
j = sns.jointplot(x= nutral.F, y= nutral.mz_neutral, kind='hex', bins=10)
j.annotate(stats.pearsonr)


# %% [markdown]
#
# # Apply conditional statements to var_ds to find a time of stable unstable surface layer. 

# ### $$M(z) = M_B * (({\zeta}^D)^A)*exp[A*(1-\zeta ^D)]$$
# ### $$\zeta = (1/C)*(z/z_i)*(w_*/u_*)^B$$


# %% 

###############################################################################################
"""################## Unstable Conditions using Radix Wind Equation #######################"""
###############################################################################################


#######  THIS SECTION NEEDS WORK 
## Make dummy data for w_str casue i dont have any :(
## Data range was taken from stulls Minessota field campaign
lengeth = len(np.array(unstable_ii.P0))
w_str = np.random.uniform(low=1., high=3.5, size=(lengeth,))
z_i = np.random.uniform(low=180., high=250, size=(lengeth,))
unstable_ii.update({'w_str':(('time'), w_str)})
unstable_ii.update({'z_i':(('time'), z_i)})

dim = [] 

## I hate this loop its slow and stupid
for i in range(len(z_list)):
    dim_i = wind_eq.dimRxL(unstable_ii.w_str,unstable_ii.UST, z_list[i], unstable_ii.z_i)
    dim.append(dim_i)

dim_stack = np.stack(dim)

## Add to Stable DataArray
unstable_ii.update({'dim':(('time','z'), dim_stack.T)})

### Drop any and all dim events less than zero or greater than 1
unstable_iii = unstable_ii.where((unstable_ii.dim > 0) & (unstable_ii.dim < 100000), drop=False)

min_rain, max_rain = np.nanmin(unstable_iii.dim), np.nanmax(unstable_iii.dim)
print(f"Min Dim {min_rain} cm & Max Dim {max_rain}")
#
sum_me = np.array(unstable_iii.F[:,-4])
unstable_non_nans = (~np.isnan(sum_me)).sum()
print(f'We have {unstable_non_nans} number of times of a Unstable BL')

mz_unstable = wind_eq.RxL(unstable_iii.dim,unstable_iii.F[:,0])

### Add to Stable DataArray
unstable_iii.update({'mz_unstable':(('time','z'), mz_unstable)})
#unstable_iv  = unstable_iii.where(unstable_iii.dim.any() > 1 ,unstable_iii.mz_unstable.any(), unstable_iii.F[:,0].any())


### Remove any odd negative values...some occur and I dont know why 
unstable = unstable_iii.where((unstable_iii.mz_unstable > 0), drop=False)

sum_me = np.array(unstable.F[:,-4])
unstable_non_nans = (~np.isnan(sum_me)).sum()
print(f'We have {unstable_non_nans} number of times of a Unstable BL')


# %% [markdown]
#
# # Apply conditional statements to var_ds to find  natural surface layer. 


# %% 
## Drop any and all rain events
min_dim, max_dim, max_dim_ave = np.nanmin(unstable.dim), np.nanmax(unstable.dim), np.nanmean(unstable.dim)
print(f"Min Dim {min_rain} cm & Max Dim {max_rain} & Avg Dim {max_dim_ave}")

min_mz_unstable, max_mz_unstable = np.nanmin(unstable.mz_unstable), np.nanmax(unstable.mz_unstable)
print(f"Min mz {min_mz_unstable} cm & Max mz {max_mz_unstable}")

# %%

###############################################################################################
"""####################### Plot Unstable Conditions mean and do Stats  #########################"""
###############################################################################################
fig, ax = plt.subplots(1,1, figsize=(12,10))
fig.suptitle('Log Wind Profile in Unstable Radix Layer', fontsize= plt_set.title_size, fontweight="bold")


ax.scatter(unstable.F.mean(dim='time'),z_list, color = 'red', marker='+', label = 'in-situ')
ax.plot(unstable.mz_unstable.mean(dim='time'),z_list, color = 'k', label = 'model')

ax.set_xlabel("Wind Speed $(ms^-1)$", fontsize= plt_set.label)
ax.set_ylabel("Height Above Ground Level  \n (AGL)", fontsize= plt_set.label)
ax.tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
ax.xaxis.grid(color='gray', linestyle='dashed')
ax.yaxis.grid(color='gray', linestyle='dashed')
ax.set_facecolor('lightgrey')
ax.legend(loc='best')


fig.savefig(save + 'Radix_Layer_Wind_Profile')


# %%
## Hexbin plot or 2D histogram
j = sns.jointplot(x= unstable.F, y= unstable.mz_unstable, kind='hex', bins=10)
j.annotate(stats.pearsonr)








