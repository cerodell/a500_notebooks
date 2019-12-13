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

from cr500.utils import plt_set, constants, wind_eq, read_cabauw
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
    m_z_i = wind_eq.loglin_stable(np.array(stable_i.UST),z_list[i],np.array(stable_i.L))
    m_z.append(m_z_i)
m_z = np.stack(m_z)
m_z[m_z < 0] = np.nan

## Add to Stable DataArray
stable_i.update({'mz_stable':(('time','z'), m_z.T)})

## Remove any odd negative values...some occur and I dont know why 
#stable = stable_i.where((stable_i.mz_stable > 0) & (stable_i.mz_stable < 40), drop=False)
stable = stable_i

sum_me = np.array(stable.F[:,-4])
stable_non_nans = (~np.isnan(sum_me)).sum()
print(f'We have {stable_non_nans} number of times of a Stable BL')

# %% [markdown]
#
# # Plot the mean wind speed at height for all stable conditions using Log Linear Wind Equation

# %%

###############################################################################################
"""####################### Plot Stable Conditions mean and do Stats  #########################"""
###############################################################################################

fig, ax = plt.subplots(1,1, figsize=(12,10))
fig.suptitle('Wind Profile in Stable Surface Layer \n Average 2001 - 2018', fontsize= plt_set.title_size, fontweight="bold")

#ax.scatter(stable.F.mean(dim='time'),z_list, color = 'red', marker='+', label = 'in-situ')
ax.plot(stable.F.mean(dim='time'),z_list, color = 'red', linestyle = '--', label = 'in-situ')
ax.plot(stable.mz_stable.mean(dim='time'),z_list, color = 'k', label = 'model')

ax.set_xlabel("Wind Speed $(ms^-1)$", fontsize= plt_set.label)
ax.set_ylabel("Height Above Ground Level  \n (AGL)", fontsize= plt_set.label)
ax.tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
ax.xaxis.grid(color='gray', linestyle='dashed')
ax.yaxis.grid(color='gray', linestyle='dashed')
ax.set_facecolor('lightgrey')
ax.legend(loc='best')
   
fig.savefig(save + 'Stable_Surface_Layer')


# %% [markdown]
#
# # Compare the Log Linear Wind Equation to Obervational Tower Data  in a Stable ABL

# %%
## Hexbin plot or 2D histogram
j = sns.jointplot(x= stable.F, y= stable.mz_stable, kind='hex', bins=15)
j.annotate(stats.pearsonr)
j.fig.set_size_inches(8,8)
plt.subplots_adjust(top=0.9)
j.fig.suptitle('Stable Conditions \n  Log Linear Wind Eq (model) to Observation (in-situ)')
j.savefig(save + "Stats_Stable.png")


## Solve for the root mean square error
rmse = wind_eq.rmse(stable.mz_stable , stable.F)
print(f"The root mean square error for Stable Surface Layer is {rmse}")

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
m_z[m_z < 0] = np.nan

## Add to Stable DataArray
nutral_i.update({'mz_neutral':(('time','z'), m_z.T)})

## Remove any odd negative values...some occur and I dont know why 
#nutral = nutral_i.where((nutral_i.mz_neutral > 0), drop=False)
nutral = nutral_i
sum_me = np.array(nutral.F[:,-4])
nutral_non_nans = (~np.isnan(sum_me)).sum()
print(f'We have {nutral_non_nans} number of times of a Nutral BL')


# %% [markdown]
#
# # Plot the mean wind speed at height for all neutral conditions using Log Wind Equation and tower observations


# %%

###############################################################################################
"""####################### Plot Neutral Conditions mean and do Stats  #########################"""
###############################################################################################
fig, ax = plt.subplots(1,1, figsize=(12,10))
fig.suptitle('Wind Profile in Nutrual Surface Layer', fontsize= plt_set.title_size, fontweight="bold")


#ax.scatter(nutral.F.mean(dim='time'),z_list, color = 'red', marker='+', label = 'in-situ')
ax.plot(nutral.F.mean(dim='time'),z_list, color = 'red', linestyle='--', label = 'in-situ')
ax.plot(nutral.mz_neutral.mean(dim='time'),z_list, color = 'k', label = 'model')

ax.set_xlabel("Wind Speed $(ms^-1)$", fontsize= plt_set.label)
ax.set_ylabel("Height Above Ground Level  \n (AGL)", fontsize= plt_set.label)
ax.tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
ax.xaxis.grid(color='gray', linestyle='dashed')
ax.yaxis.grid(color='gray', linestyle='dashed')
ax.set_facecolor('lightgrey')
ax.legend(loc='best')
   
fig.savefig(save + 'Nutrual_Surface_Layer')


# %% [markdown]
#
# # Compare the Log Wind Equation to Obervational Tower Data  in a Neutral ABL
# %%
## Hexbin plot or 2D histogram
j = sns.jointplot(x= nutral.F, y= nutral.mz_neutral, kind='hex', bins=15)
j.annotate(stats.pearsonr)
j.fig.set_size_inches(8,8)
plt.subplots_adjust(top=0.9)
j.fig.suptitle('Neutral Conditions \n Log Wind Eq (model) to Observation (in-situ)')
j.savefig(save + "Stats_Nutrual.png")


## Solve for the root mean square error
rmse = wind_eq.rmse(nutral.mz_neutral , nutral.F)
print(f"The root mean square error for Neutral Surface Layer is {rmse}")

# %% [markdown]
#
# # Apply conditional statements to var_ds to find a time of stable unstable surface layer. 
#
# ### $$M(z) = M_B * (({\zeta}^D)^A)*exp[A*(1-\zeta ^D)]$$
# ### $$\zeta = (1/C)*(z/z_i)*(w_*/u_*)^B$$


# %% [markdown]

# ### Make dummy data for w_str casue i dont have any :(
# ### Data range was taken from stulls Minessota field campaign

# %%

###############################################################################################
"""################## Unstable Conditions using Radix Wind Equation #######################"""
###############################################################################################


#######  THIS SECTION NEEDS WORK 
## Make dummy data for w_str casue i dont have any :(
## Data range was taken from stulls Minessota field campaign
lengeth = len(np.array(unstable_ii.P0))
w_str = np.random.uniform(low=1., high=3.5, size=(lengeth,))
z_i = np.random.uniform(low=1000., high=1500, size=(lengeth,))
unstable_ii.update({'w_str':(('time'), w_str)})
unstable_ii.update({'z_i':(('time'), z_i)})

#print(np.nanmax(unstable_ii.F[:,0]))

## Make U_str not equal zero causes calc to go to inf later :/
unstable_iii = unstable_ii.where((unstable_ii.UST > 0.000001), drop=False)
min_u_str, max_u_str = np.nanmin(unstable_iii.UST), np.nanmax(unstable_iii.UST)
print(f"Min UST {min_u_str} ms^-1 & Max UST {max_u_str} ms^-1")


# %% [markdown]

# # Solve for wind in the Radix Layer

# %%

dim = [] 

## I hate this loop its slow and stupid
for i in range(len(z_list)):
    dim_i = wind_eq.dimRxL(unstable_iii.w_str,unstable_iii.UST, z_list[i], unstable_iii.z_i)
    dim.append(dim_i)

dim_stack = np.stack(dim)

dim_stack[dim_stack < 0] = np.nan

dim_fl = dim_stack
dim_fl[dim_fl > 1] = np.nan


wsp_z = np.array(unstable_iii.F)

wsp_z[dim_fl.T > 1]

### Add dim to data array
unstable_iii.update({'dim':(('time','z'), dim_fl.T)})

mz_unstable = wind_eq.RxL(unstable_iii.dim,unstable_iii.F[:,0])

a = mz_unstable[:,4:]  
b = (wsp_z[:,0:4])
please_work = np.concatenate((b, a), axis=1)


unstable_iii.update({'mz_unstable':(('time','z'), mz_unstable)})
unstable = unstable_iii



# %% [markdown]
#
# # Plot the mean wind speed at height for all unstable conditions using Radix Wind Equation compared to tower observations


# %%

###############################################################################################
"""####################### Plot Unstable Conditions mean and do Stats  #########################"""
###############################################################################################
fig, ax = plt.subplots(1,1, figsize=(12,10))
fig.suptitle('Wind Profile in Unstable Radix Layer', fontsize= plt_set.title_size, fontweight="bold")


#ax.scatter(unstable.F.mean(dim='time'),z_list, color = 'red', marker='+', label = 'in-situ')
ax.plot(unstable.F.mean(dim='time'),z_list, color = 'red', linestyle = '--', label = 'in-situ')
ax.plot(unstable.mz_unstable.mean(dim='time'),z_list, color = 'k', label = 'model')

ax.set_xlabel("Wind Speed $(ms^-1)$", fontsize= plt_set.label)
ax.set_ylabel("Height Above Ground Level  \n (AGL)", fontsize= plt_set.label)
ax.tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
ax.xaxis.grid(color='gray', linestyle='dashed')
ax.yaxis.grid(color='gray', linestyle='dashed')
ax.set_facecolor('lightgrey')
ax.legend(loc='best')


fig.savefig(save + 'Unstable_Wind_Profile')



# %% [markdown]
#
# # Compare the Radix Layer Wind Equation to Obervational Tower Data in Unstable  condtions

# %%
## Hexbin plot or 2D histogram

#sns.set(rc={'figure.figsize':(11.7,8.27)})
j = sns.jointplot(x= unstable.F, y= unstable.mz_unstable, kind='hex', bins=15)
j.annotate(stats.pearsonr)
j.fig.set_size_inches(8,8)
plt.subplots_adjust(top=0.9)
j.fig.suptitle('Unstable Conditions \n Radix Eq (model) to Observation (in-situ)')
j.savefig(save + "Stats_Unstable.png")




## Solve for the root mean square error
rmse = wind_eq.rmse(unstable.mz_unstable , unstable.F)
print(f"The root mean square error for unstable Radix Layer is {rmse}")


# %% [markdown]
#
# # With the Radix Equation Failing misserably we will instead solve a nonlinear function in variable (a) and to fit the observed tower data


# %%


def wind_func(z, a0, a1, a2, a3):
    'nonlinear function in a and to fit to data'
    fit = a0 + a1*z + a2*z**2. + a3*np.log(z)
    return fit

rev_z = z_list[::-1]
rev_z = rev_z[1:]

## Make and Array for wind speed form the tower
wsp = np.array(unstable.F)
## Flip Array by its z dimiention
wsp = wsp[:,:-1]

## Create a mask of for any location of nan valuse
mask = np.any(np.isnan(wsp), axis=1)

## apply mask to remove nan values
wsp_for_fit = wsp[~mask]

## Initialize a list of wind speed as caluclated by the fit function
mz_fit_i = []

poptf, pcovf = [], []
## Loop and caluclate the fit function for every time
for i in range(len(wsp_for_fit[:,0])):
    popt, pcov = curve_fit(wind_func,z_list[:-1],wsp_for_fit[i,:]) 
    #### Create a z axis of the same shape as function and run funtion 
    poptf.append(popt)
    pcovf.append(pcov)
    
    fit_f = wind_func(z_list,*popt)
    mz_fit_i.append(fit_f)
    
## stack that shit 
mz_fit = np.stack(mz_fit_i)

# %% [markdown]
#
# # Plot the mean wind speed at height for all unstable conditions using WInd Fit Equation compared to tower observations


# %%
###############################################################################################
"""####################### Plot Unstable Conditions mean and do Stats  #########################"""
###############################################################################################
fig, ax = plt.subplots(1,1, figsize=(12,10))
fig.suptitle('Wind Profile in Unstable Radix Layer', fontsize= plt_set.title_size, fontweight="bold")


#ax.scatter(unstable.F.mean(dim='time'),z_list, color = 'red', marker='+', label = 'in-situ')
ax.plot(unstable.F.mean(dim='time'),z_list, color = 'red', linestyle = '--', label = 'in-situ')
ax.plot(mz_fit.mean(axis = 0),z_list, color = 'k', label = 'model')

ax.set_xlabel("Wind Speed $(ms^-1)$", fontsize= plt_set.label)
ax.set_ylabel("Height Above Ground Level  \n (AGL)", fontsize= plt_set.label)
ax.tick_params(axis='both', which='major', labelsize= plt_set.tick_size)
ax.xaxis.grid(color='gray', linestyle='dashed')
ax.yaxis.grid(color='gray', linestyle='dashed')
ax.set_facecolor('lightgrey')
ax.legend(loc='best')


fig.savefig(save + 'Unstable_Wind_Profile_Fitted')

# %% [markdown]

# # Compare the Wind Fit Function to Obervational Tower Data in Unstable condtions.....of no suprise its perfect :) 

# %%
## Hexbin plot or 2D histogram

sns.set(rc={'figure.figsize':(11.7,8.27)})
j = sns.jointplot(x=wsp_for_fit , y= mz_fit[:,:-1] , kind='hex', bins=15)
j.annotate(stats.pearsonr)
j.fig.set_size_inches(8,8)
plt.subplots_adjust(top=0.9)
j.fig.suptitle('Unstable Conditions \n  Wind Fit Function (model) to Observation (in-situ)')
j.savefig(save + "Stats_Unstable_Fitted.png")


## Solve for the root mean square error
rmse = wind_eq.rmse(mz_fit[:,:-1],wsp_for_fit)
print(f"The root mean square error for unstable fit is {rmse}")

# %% [markdown]

# # Make an Xarray of the new fitted fuction ...this could be better

# %% 

## declair you dimentions for fitted xr and make dict for each variable you will add
dims_fit = ('time', 'z')
fit_dict1 = {'wsp_fit' : (dims_fit,np.array(wsp_for_fit))}
fit_dict2 = {'mz_fit' : (dims_fit,np.array(mz_fit[:,:-1]))}

## place them in a list
fit_list = [fit_dict1, fit_dict2]

## make the xr
fit_ds = read_cabauw.xarray_like(fit_list)


# %% [markdown]
#
# # Look a the seasonality difference of the varied wind profile equations

# #### This shows mean wind speed profile of the observations versus the wind profile equations 

# %%


import pandas as pd

seasons_stable = stable.groupby('time.season').mean('time')
seasons_nutral = nutral.groupby('time.season').mean('time')
seasons_unstable = unstable.groupby('time.season').mean('time')
#seasons_fitted = fit_ds.groupby('time.season').mean('time')


# Quick plot to show the results
notnull = pd.notnull(seasons_stable.F[0])

color_season = ['blue', 'green', 'red', 'orange']
fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(12,10))
for i, season in enumerate(('DJF', 'MAM', 'JJA', 'SON')):
    seasons_stable['F'].sel(season=season).where(notnull).plot(
        ax=axes[i, 0], y= 'z', color = color_season[i], linestyle='--')
    
    seasons_stable['mz_stable'].sel(season=season).where(notnull).plot(
        ax=axes[i, 0], y= 'z', color = color_season[i])
    
    seasons_nutral['F'].sel(season=season).where(notnull).plot(
        ax=axes[i, 1], y= 'z', color = color_season[i], linestyle='--')
    
    seasons_nutral['mz_neutral'].sel(season=season).where(notnull).plot(
        ax=axes[i, 1], y= 'z', color = color_season[i])

    seasons_unstable['F'].sel(season=season).where(notnull).plot(
        ax=axes[i, 2], y= 'z', color = color_season[i], linestyle='--')
    
    seasons_unstable['mz_unstable'].sel(season=season).where(notnull).plot(
        ax=axes[i, 2], y= 'z', color = color_season[i])


## Thought i could ge the the fitted function to work in here but had issue with the time....like i ran out of time haha 
#    seasons_fitted['wsp_fit'].sel(season=season).where(notnull).plot(
#        ax=axes[i, 3], y= 'z', color = color_season[i], linestyle='--')
#    
#    seasons_fitted['mz_fit'].sel(season=season).where(notnull).plot(
#        ax=axes[i, 3], y= 'z', color = color_season[i])

    axes[i, 0].set_ylabel(season)
    axes[i, 1].set_ylabel('')
    axes[i, 2].set_ylabel('')
#    axes[i, 3].set_ylabel('')






for ax in axes.flat:
    ax.axes.get_xaxis().set_ticklabels([])
    ax.axes.get_yaxis().set_ticklabels([])
    ax.axes.axis('tight')
    ax.axes.set_xlim(2,12)
    ax.axes.set_xlabel('')
    ax.axes.set_title('')

axes[0, 0].set_title('Stable')
axes[0, 1].set_title('Neutral')
axes[0, 2].set_title('Unstable')
#axes[0, 3].set_title('Unstable Fitted')

plt.tight_layout()

fig.savefig(save + 'Seasonal_Mean_Wind_Profiles')

# %% [markdown]
#
# # Look a the seasonality difference of the varied wind profile equations

# #### This shows the distribution of wind speed profiles of the observations versus the wind profile equations 

# #### I feel I am close to getting this but again time time time isnt on my side.. sorry mic


# %%

#
#sns.set(style="dark")
#import pandas as pd

seasons_stable_all = list(stable.groupby('time.season'))
seasons_unstable_all = list(unstable.groupby('time.season'))
seasons_nutral_all = list(nutral.groupby('time.season'))


#tester = np.array(seasons_stable_all[0][1].mz_stable)

# %%

#color_season = ['blue', 'green', 'red', 'orange']
fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(12,10))
seasons_list = ['DJF', 'JJA', 'MAM', 'SON']

for i in range(len(seasons_list)):
    print(i)
#     Create a cubehelix colormap to use with kdeplot
    cmap = 'Reds'
##
    axes[i,0].scatter(np.array(seasons_stable_all[i][1].F) ,  np.array(seasons_stable_all[i][1].mz_stable), alpha=0.3)
#    
    axes[i,1].scatter(seasons_nutral_all[i][1].F,  seasons_nutral_all[i][1].mz_neutral , alpha=0.3)
    
    axes[i,2].scatter(seasons_unstable_all[i][1].F,  seasons_unstable_all[i][1].mz_unstable ,alpha=0.3)

#
#    axes[i, 0].set_ylabel(season)
#    axes[i, 1].set_ylabel('')
#    axes[i, 2].set_ylabel('')
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