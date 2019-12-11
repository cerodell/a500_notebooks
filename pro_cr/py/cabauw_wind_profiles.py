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
from scipy.optimize import curve_fit
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
#
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
# # TRY TO GET DOPPLER DATA WORKING BUT I AM LOOSING THE BATTLE
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
timedop = np.array(doppler_ds.time)
dtime = np.diff(timedop)
plt.plot(dtime)
#print(np.min(dtime))
#ztime = np.sum(np.diff(timedop))
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
meso_ds_i = var_ds.where((var_ds.RAIN < 0.0001), drop=False)
min_rain, max_rain = np.nanmin(meso_ds_i.RAIN), np.nanmax(meso_ds_i.RAIN)
print(f"Min Rain {min_rain} cm & Max Rain {max_rain} cm")

## Drop low press events....should likely expand on this of press gradient over time.....
meso_ds = meso_ds_i.where((meso_ds_i.P0 > 990. ), drop=False)
min_press, max_press = round(np.nanmin(meso_ds.P0),3), round(np.nanmax(meso_ds.P0),3)
print(f"Min Press {min_press} hPa & Max Press {max_press} hPa")


# %% [markdown]
#
# ### Filtter for dirrefnt stability classes by comapring dT/dz to the dry adiabtic laps rate

# %%

# Apply conditional statements
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
#I have been tesing this method and it yeild little differance as comapred 
#to my other arropch to define stability

# Apply conditional statements
phi = (meso_ds.z[-1]/meso_ds.L)
meso_ds.update({'phi':(('time'), phi)})
min_phi, max_phi = round(np.nanmin(meso_ds.L),3), round(np.nanmax(meso_ds.L),3)
print(f"Min (z/L) {min_phi} hPa & Max (z/L) {max_phi}")

stable_i = var_ds.where((meso_ds.phi > 0.), drop=False)

nutral_i = var_ds.where((meso_ds.phi > -0.009) & (meso_ds.phi < 0.009), drop=False)

unstable_ii = var_ds.where(meso_ds.phi < 0. , drop=False)

"""
sum_me = np.array(stable_i.F[:,-4])
nutral_non_nans = (~np.isnan(sum_me)).sum()
print(f'We have {nutral_non_nans} number of times of a Stable BL')

sum_me = np.array(nutral_i.F[:,-4])
nutral_non_nans = (~np.isnan(sum_me)).sum()
print(f'We have {nutral_non_nans} number of times of a Neutral BL')

sum_me = np.array(unstable_ii.F[:,-4])
nutral_non_nans = (~np.isnan(sum_me)).sum()
print(f'We have {nutral_non_nans} number of times of a Unstable BL')
# %% [markdown]
#
# # Show the Distribution of dT/dz for each stability conditon


# %%

f, axes = plt.subplots(1, 3, figsize=(16, 8), sharex=True)
sns.kdeplot(stable_i.dTdz_mean, shade=True, color="r", ax=axes[0], label = 'Stable')
sns.kdeplot(nutral_i.dTdz_mean, shade=True, color="b", ax=axes[1], label = 'Nutral')
sns.kdeplot(unstable_ii.dTdz_mean, shade=True, color="g", ax=axes[2], label = 'Unstable')
sns.set_style("darkgrid", {"axes.facecolor": ".9"})



# %% [markdown]
#
# # Solve for wind speed at height in stable conditions using Log Linear Wind Equation
#
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
fig.suptitle('Wind Profile in Stable Surface Layer \n Average 2001 - 2018', fontsize= plt_set.title_size, fontweight="bold")

ax.scatter(stable.F.mean(dim='time'),z_list, color = 'red', marker='+', label = 'in-situ')
ax.plot(stable.mz_stable.mean(dim='time'),z_list, color = 'k', label = 'model')

ax.set_xlabel("Wind Speed $(ms^-1)$", fontsize= plt_set.label)
ax.set_ylabel("Height Above Ground Level  \n (AGL)", fontsize= plt_set.label)
ax.tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
ax.xaxis.grid(color='gray', linestyle='dashed')
ax.yaxis.grid(color='gray', linestyle='dashed')
ax.set_facecolor('lightgrey')
ax.legend(loc='best')
   
fig.savefig(save + 'Stable_Surface_Layer')


# %%
## Hexbin plot or 2D histogram
j = sns.jointplot(x= stable.F, y= stable.mz_stable, kind='hex', bins=15)
j.annotate(stats.pearsonr)
j.fig.set_size_inches(8,8)
plt.subplots_adjust(top=0.9)
j.fig.suptitle('Stable Conditions \n  Log Linear Wind Eq (model) to Observation (in-situ)')
j.savefig(save + "Stats_Stable.png")








# %% [markdown]
#
# # Solve for wind speed at height in nutral conditions using Log  Wind Equation
#
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
fig.suptitle('Wind Profile in Nutrual Surface Layer', fontsize= plt_set.title_size, fontweight="bold")


ax.scatter(nutral.F.mean(dim='time'),z_list, color = 'red', marker='+', label = 'in-situ')
ax.plot(nutral.mz_neutral.mean(dim='time'),z_list, color = 'k', label = 'model')

ax.set_xlabel("Wind Speed $(ms^-1)$", fontsize= plt_set.label)
ax.set_ylabel("Height Above Ground Level  \n (AGL)", fontsize= plt_set.label)
ax.tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
ax.xaxis.grid(color='gray', linestyle='dashed')
ax.yaxis.grid(color='gray', linestyle='dashed')
ax.set_facecolor('lightgrey')
ax.legend(loc='best')
   
fig.savefig(save + 'Nutrual_Surface_Layer')

# %%
## Hexbin plot or 2D histogram
j = sns.jointplot(x= nutral.F, y= nutral.mz_neutral, kind='hex', bins=15)
j.annotate(stats.pearsonr)
j.fig.set_size_inches(8,8)
plt.subplots_adjust(top=0.9)
j.fig.suptitle('Nutrual Conditions \n Log Wind Eq (model) to Observation (in-situ)')
j.savefig(save + "Stats_Nutrual.png")


# %% [markdown]
#
# # Apply conditional statements to var_ds to find a time of stable unstable surface layer. 
#
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

print(np.nanmax(unstable_ii.F[:,0]))

## Make U_str not equal zero
unstable_iii = unstable_ii.where((unstable_ii.UST > 0.000001), drop=False)
min_u_str, max_u_str = np.nanmin(unstable_iii.UST), np.nanmax(unstable_iii.UST)
print(f"Min UST {min_u_str} ms^-1 & Max UST {max_u_str} ms^-1")

# %%




## Make play data for height of loglin model
#height = np.arange(10,201,1)
m_z, dim = [] , []

## I hate this loop its slow and stupid
for i in range(len(z_list)):
    mz_unstable, dim_i = wind_eq.RxLtest(unstable_iii.w_str, unstable_iii.UST, z_list[i],unstable_iii.z_i, unstable_iii.F[:,0])
    dim.append(dim_i)
    m_z.append(m_z_i)

dim = np.stack(dim)   
m_z = np.stack(m_z)


### Add to radix eq derived wsp to DataArray
unstable_iii.update({'dim':(('time','z'), dim.T)})
unstable_iii.update({'mz_unstable':(('time','z'), m_z.T)})

max_f = np.nanmax(unstable_iii.mz_unstable[:,0])
print(f"If NAN you are fucked --> {max_f}")

### Drop any and all dim events less than zero or greater than 1
unstable_iv = unstable_iii.where((unstable_iii.dim > 0), drop=False)
#unstable_v = unstable_iv.where(unstable_iv.dim < 1, unstable_iv.mz_unstable, unstable_iv.F[:,0].all())

unstable_v = unstable_iv
#min_dim, max_dim, avg_dim = np.nanmin(unstable_iv.dim), np.nanmax(unstable_iv.dim), np.nanmean(unstable_iv.dim)
min_dim, max_dim, avg_dim = np.nanmin(unstable_v.dim), np.nanmax(unstable_v.dim), np.nanmean(unstable_v.dim)
print(f"Min Dim {min_dim} & Max Dim {max_dim} & Avg Dim {avg_dim}")

max_f = np.nanmax(unstable_v.F[:,0])
print(f"If NAN you are fucked --> {max_f}")
#plt.plot(unstable_iv.mz_unstable)
#### Remove any odd negative values...some occur and I dont know why 
unstable = unstable_v.where((unstable_v.mz_unstable > 0), drop=False)

###Check to see how many cases we have left
sum_me = np.array(unstable.F[:,-4])
unstable_non_nans = (~np.isnan(sum_me)).sum()
print(f'We have {unstable_non_nans} number of times of a Unstable BL')

max_f = np.nanmax(unstable.F[:,0])
print(f"If NAN you are fucked --> {max_f}")


# %%
#
#dim = [] 
#
### I hate this loop its slow and stupid
#for i in range(len(z_list)):
#    dim_i = wind_eq.dimRxL(unstable_iii.w_str,unstable_iii.UST, z_list[i], unstable_iii.z_i)
#    dim.append(dim_i)
#
#dim_stack = np.stack(dim)
#
### Add dim to data array
#unstable_iii.update({'dim':(('time','z'), dim_stack.T)})
#
##min_dim, max_dim, avg_dim = np.nanmin(unstable_iii.dim[:,2]), np.nanmax(unstable_iii.dim[:,2]), np.nanmean(unstable_iii.dim[:,2])
#min_dim, max_dim, avg_dim = np.nanmin(unstable_iii.dim), np.nanmax(unstable_iii.dim), np.nanmean(unstable_iii.dim)
#print(f"Min Dim {min_dim} & Max Dim {max_dim} & Avg Dim {avg_dim}")
##unstable_iv = unstable_iii
#### Drop any and all dim events less than zero or greater than 1
#unstable_iv = unstable_iii.where((unstable_iii.dim > 0), drop=False)
#unstable_v = unstable_iv.where(unstable_iv.dim > 1, unstable_iv.

#min_dim, max_dim, avg_dim = np.nanmin(unstable_iv.dim[:,-1]), np.nanmax(unstable_iv.dim), np.nanmean(unstable_iv.dim)
#print(f"Min Dim {min_dim} & Max Dim {max_dim} & Avg Dim {avg_dim}")

#unstable_v = unstable_iv.where((unstable_iv.F[:,0] > 2), drop=False)
#
###Check to see how many cases we have left
#sum_me = np.array(unstable_v.F[:,-4])
#unstable_non_nans = (~np.isnan(sum_me)).sum()
#print(f'We have {unstable_non_nans} number of times of a Unstable BL')
#
#
#
### Solve for wsp using the radix wind eq!!!
#mz_unstable = wind_eq.RxL(unstable_v.dim,unstable_v.F[:,0])
#
#### Add to radix eq derived wsp to DataArray
#unstable_v.update({'mz_unstable':(('time','z'), mz_unstable)})
#
#print(np.nanmax(unstable_v.F[:,0]))
#print(np.nanmax(unstable_v.z[-4]))
#
#
##plt.plot(unstable_iv.mz_unstable)
##### Remove any odd negative values...some occur and I dont know why 
#unstable = unstable_v.where((unstable_v.mz_unstable > 0), drop=False)
#
####Check to see how many cases we have left
#sum_me = np.array(unstable.F[:,-4])
#unstable_non_nans = (~np.isnan(sum_me)).sum()
#print(f'We have {unstable_non_nans} number of times of a Unstable BL')


# %% [markdown]
#
# # Apply conditional statements to var_ds to find  natural surface layer. 


# %%
## Drop any and all rain events
min_dim, max_dim, max_dim_ave = round(np.nanmin(unstable.dim),3), round(np.nanmax(unstable.dim),3), round(np.nanmean(unstable.dim),3)
print(f"Min Dim {min_dim} & Max Dim {max_dim} & Avg Dim {max_dim_ave}")

min_mz_unstable, max_mz_unstable = round(np.nanmin(unstable.mz_unstable),3), round(np.nanmax(unstable.mz_unstable),3)
print(f"Min mz {min_mz_unstable} cm & Max mz {max_mz_unstable}")

# %%

###############################################################################################
"""####################### Plot Unstable Conditions mean and do Stats  #########################"""
###############################################################################################
fig, ax = plt.subplots(1,1, figsize=(12,10))
fig.suptitle('Wind Profile in Unstable Radix Layer', fontsize= plt_set.title_size, fontweight="bold")


ax.scatter(unstable.F.mean(dim='time'),z_list, color = 'red', marker='+', label = 'in-situ')
ax.plot(unstable.mz_unstable.mean(dim='time'),z_list, color = 'k', label = 'model')

ax.set_xlabel("Wind Speed $(ms^-1)$", fontsize= plt_set.label)
ax.set_ylabel("Height Above Ground Level  \n (AGL)", fontsize= plt_set.label)
ax.tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
ax.xaxis.grid(color='gray', linestyle='dashed')
ax.yaxis.grid(color='gray', linestyle='dashed')
ax.set_facecolor('lightgrey')
ax.legend(loc='best')


fig.savefig(save + 'Unstable_Wind_Profile')


# %%
## Hexbin plot or 2D histogram

sns.set(rc={'figure.figsize':(11.7,8.27)})
j = sns.jointplot(x= unstable.F, y= unstable.mz_unstable, kind='hex', bins=10)
j.annotate(stats.pearsonr)
j.fig.set_size_inches(8,8)
plt.subplots_adjust(top=0.9)
j.fig.suptitle('Unstable Conditions \n Radix Eq (model) to Observation (in-situ)')
j.savefig(save + "Stats_Unstable.png")



# %%
#
##def wind_func(z, a0, a1, a2, a3):
##    'nonlinear function in a and to fit to data'
##    fit = a0 + a1*z + a2*z**2. + a3*np.log(z)
##    return fit
#
#unstable_F = unstable.where(unstable.F[:,-4] == np.nan, drop=False)
#
#def RxL(dim,m_bl, A, D):  
#    '''
#    This Function is used for the wind profile in the Radix Layer
#
#    '''   
#    m_z = m_bl * ((dim**D)**A) * np.exp(A*(1-dim**D))
#    m_z = np.array(m_z)
#    return m_z
#
#    
#def rmnan(var):
#    a = np.array(var)
#    a = a[~np.isnan(a)]
#    a[a < 1E308]
#    return a
#
#wsp     = rmnan(unstable.F.mean(dim='time'))
#wsp_200 = rmnan(unstable.F[:,0])
#u_str   = rmnan(unstable.UST)
#w_str   = rmnan(unstable.w_str)
#dim_z   = rmnan(unstable.dim)
#
#mx_hiehgt = rmnan(unstable_iii.z_i)
#
#
#zz = unstable_iii.z[0:6]
#
##### Do the curve fit 
#popt, pcov = curve_fit(RxL,zz,wsp)
#print(popt.shape)
#
##### Create a z axis of the same shape as function and run funtion 
##z_int = np.arange(1,201,0.5)
#fit_f = RxL(dim_z, wsp_200, *popt)
#
## %%
#
#fig, ax = plt.subplots(1,1, figsize=(12,10))
#fig.suptitle('Wind Profile in Unstable Radix Layer', fontsize= plt_set.title_size, fontweight="bold")
#ax.plot(fit_f,z_list, color ='k', linewidth = 2)



# %% [markdown]
#
# # Look a the seasonality difference of the varied wind profile equations

# #### This shows mean wind speed profile of the observations versus the wind profile equations 

# %%


import pandas as pd

seasons_stable = stable.groupby('time.season').mean('time')
seasons_unstable = unstable.groupby('time.season').mean('time')
seasons_nutral = nutral.groupby('time.season').mean('time')


# Quick plot to show the results
notnull = pd.notnull(seasons_stable.F[0])

color_season = ['blue', 'green', 'red', 'orange']
fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(12,10))
for i, season in enumerate(('DJF', 'MAM', 'JJA', 'SON')):
    seasons_stable['F'].sel(season=season).where(notnull).plot(
        ax=axes[i, 0], y= 'z', color = color_season[i])
    
    seasons_stable['mz_stable'].sel(season=season).where(notnull).plot(
        ax=axes[i, 0], y= 'z', color = color_season[i], linestyle='dashed')

    seasons_unstable['F'].sel(season=season).where(notnull).plot(
        ax=axes[i, 1], y= 'z', color = color_season[i])
    
    seasons_unstable['mz_unstable'].sel(season=season).where(notnull).plot(
        ax=axes[i, 1], y= 'z', color = color_season[i], linestyle='dashed')

    seasons_nutral['F'].sel(season=season).where(notnull).plot(
        ax=axes[i, 2], y= 'z', color = color_season[i])
    
    seasons_nutral['mz_neutral'].sel(season=season).where(notnull).plot(
        ax=axes[i, 2], y= 'z', color = color_season[i], linestyle='dashed')

    axes[i, 0].set_ylabel(season)
    axes[i, 1].set_ylabel('')
    axes[i, 2].set_ylabel('')





for ax in axes.flat:
    ax.axes.get_xaxis().set_ticklabels([])
    ax.axes.get_yaxis().set_ticklabels([])
    ax.axes.axis('tight')
    ax.axes.set_xlim(2,12)
    ax.axes.set_xlabel('')
    ax.axes.set_title('')

axes[0, 0].set_title('Stable')
axes[0, 1].set_title('Unstable')
axes[0, 2].set_title('Neutral')

plt.tight_layout()

#fig.suptitle('Seasonal Wind Speed Profiles', fontsize=16, y=1.02)

# %% [markdown]
#
# # Look a the seasonality difference of the varied wind profile equations

# #### This shows the distribution of wind speed profiles of the observations versus the wind profile equations 


# %%

#
#import pandas as pd
#
seasons_stable_all = stable.groupby('time.season')
#seasons_unstable = unstable.groupby('time.season').mean('time')
#seasons_nutral = nutral.groupby('time.season').mean('time')
#
#
## Quick plot to show the results
#notnull = pd.notnull(seasons_stable.F[0])
#
#color_season = ['blue', 'green', 'red', 'orange']
#fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(12,10))
#for i, season in enumerate(('DJF', 'MAM', 'JJA', 'SON')):
#    seasons_stable['F'].sel(season=season).where(notnull).plot(
#        ax=axes[i, 0], y= 'z', color = color_season[i])
#    
#    seasons_stable['mz_stable'].sel(season=season).where(notnull).plot(
#        ax=axes[i, 0], y= 'z', color = color_season[i], linestyle='dashed')
#
#    seasons_unstable['F'].sel(season=season).where(notnull).plot(
#        ax=axes[i, 1], y= 'z', color = color_season[i])
#    
#    seasons_unstable['mz_unstable'].sel(season=season).where(notnull).plot(
#        ax=axes[i, 1], y= 'z', color = color_season[i], linestyle='dashed')
#
#    seasons_nutral['F'].sel(season=season).where(notnull).plot(
#        ax=axes[i, 2], y= 'z', color = color_season[i])
#    
#    seasons_nutral['mz_neutral'].sel(season=season).where(notnull).plot(
#        ax=axes[i, 2], y= 'z', color = color_season[i], linestyle='dashed')
#
#    axes[i, 0].set_ylabel(season)
#    axes[i, 1].set_ylabel('')
#    axes[i, 2].set_ylabel('')
#
#
#
#
#
#for ax in axes.flat:
#    ax.axes.get_xaxis().set_ticklabels([])
#    ax.axes.get_yaxis().set_ticklabels([])
#    ax.axes.axis('tight')
#    ax.axes.set_xlim(2,12)
#    ax.axes.set_xlabel('')
#    ax.axes.set_title('')
#
#axes[0, 0].set_title('Stable')
#axes[0, 1].set_title('Unstable')
#axes[0, 2].set_title('Neutral')
#
#plt.tight_layout()
#

#
#
#sns.set(style="dark")
#import pandas as pd
#
#seasons_stable_all = stable.groupby('time.season')
##seasons_unstable = unstable.groupby('time.season').mean('time')
##seasons_nutral = nutral.groupby('time.season').mean('time')
#
#
## Quick plot to show the results
#notnull = pd.notnull(seasons_stable_all.F[0])
#
#color_season = ['blue', 'green', 'red', 'orange']
#fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(12,10))
#for i, season in enumerate(('DJF', 'MAM', 'JJA', 'SON')):
#
#    # Create a cubehelix colormap to use with kdeplot
#    cmap = 'Reds'
#
#    sns.kdeplot(seasons_stable_all['F'].sel(season=season).where(notnull), 
#                seasons_stable_all['mz_stable'].sel(season=season).where(notnull), 
#                cmap=cmap, shade=True, cut=5, ax=axes[i, 0])
#
#    axes[i, 0].set_ylabel(season)
#    axes[i, 1].set_ylabel('')
#    axes[i, 2].set_ylabel('')
#
#
#
#
#
#for ax in axes.flat:
#    ax.axes.get_xaxis().set_ticklabels([])
#    ax.axes.get_yaxis().set_ticklabels([])
#    ax.axes.axis('tight')
#    ax.axes.set_xlim(2,12)
#    ax.axes.set_xlabel('')
#    ax.axes.set_title('')
#
#axes[0, 0].set_title('Stable')
#axes[0, 1].set_title('Unstable')
#axes[0, 2].set_title('Neutral')
#
##plt.tight_layout()
#f.tight_layout()