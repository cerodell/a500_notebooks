# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all
#     cell_metadata_json: true
#     formats: ipynb,py:percent
#     notebook_metadata_filter: all,-language_info,-toc,-latex_envs
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Calculating fluxes for a tropical boundary layer

# %% [markdown]
# *I changed indexing below for the mean theta and wvel calculations and removed
# the newaxis complication.   The dimensions of the velocity and temperature
# are (time,z,y,x)  or axis=(0,1,2,3).  Hopefully the sensible heat flux profile is correct
# (it's much smaller than for a dry boundary layer, because most of the energy is being
# transported by vapor in the tropical marine boundary layer*

# %% [markdown]
# ## Reading a netcdf file

# %% [markdown]
# 1.  Fetch the netcdf file tropical.nc from my webserver (50 Mbytes)

# %%
import numpy as np
from netCDF4 import Dataset
import a500
from a500.utils.data_read import download

# %% [markdown]
# Only download the data once -- set this to True if you need
# this file, otherwise skip

# %%
the_root = "http://clouds.eos.ubc.ca/~phil/docs/atsc500"
the_file =  "tropical_subset.nc"
out=download(the_file,root=the_root, dest_folder=a500.data_dir)

# %% [markdown]
# 2\.  Identify the file structure using ncdump

# %%
from a500.utils.ncdump import ncdump

the_file = a500.data_dir / the_file
with Dataset(the_file) as nc_in:
    ncdump(nc_in)

# %% [markdown]
# 3\.  Read a variable using the netCDF4 module:  http://unidata.github.io/netcdf4-python/

# %%
with Dataset(the_file,'r') as nc_in:
    print(list(nc_in.variables.keys()))
    the_temp=nc_in.variables['TABS'][...]  
    #
    # remove the time dimension since we only have one timestep
    #   
    print('temp shape',the_temp.shape)
    the_height=nc_in.variables['z'][...]
    print('height shape',the_height.shape)
    the_press=nc_in.variables['p'][...]
    the_press=the_press*100.  #convert to Pa
    wvel=nc_in.variables['W'][...]

# %% [markdown]
# # Plotting three temperature profiles

# %%
# %matplotlib inline
#
# plot 3 arbitrary columns
#
from matplotlib import pyplot as plt
meter2km=1.e-3
random_xy=[(10,20),(80,40),(25,75)]
fig,ax=plt.subplots(1,1)
for x,y in random_xy:
    temp_profile=the_temp[0,:,x,y]
    out=ax.plot(temp_profile,the_height*meter2km)  
out=ax.set(xlabel='temperature (K)',ylabel='height (km)',title='3 profiles')

# %% [markdown]
# # Calculate the horizontal average temperature

# %%
xmean=the_temp.mean(axis=3)
print(xmean.shape)
xymean=xmean.mean(axis=2)
xymean=xymean.squeeze()
print(xymean.shape)

# %% [markdown]
# # Calculate the temperature pertubation

# %%
fig,ax=plt.subplots(1,1)
for x,y in random_xy:
    temp_profile=the_temp[0,:,x,y] - xymean
    ax.plot(temp_profile,the_height*meter2km)  
ax.set(xlabel='temperature (K)',ylabel='height (km)',title='3 perturbation profiles')


# %% [markdown]
# # Now look at the fluxes

# %% [markdown]
# 1.  Write a function that takes TABS and p and returns $\theta$, the potential temperature
#
# 1.  Make a plot of the vertical profile of $\overline{w^\prime \theta^\prime}$ vs. z for
#     this model timestep
#

# %%
def theta(p,T):
    """
      input: p (Pa), scalar pressure
             T (K),  scalar or array of temperatures
      output: thetaval (K) potential temperature
    """
    Rd=287  #dry air gas constant  J/kg/K
    cp=1004.  #dry air heat capacity at constant pressure (J/kg/K)
    p0=1.e5  #reference pressure (Pa)
    the_exp=Rd/cp
    thetaval=T*(p0/p)**the_exp
    return thetaval



# %% [markdown]
# ### Calculate $\theta$ at every level

# %%
thetavals=np.empty_like(the_temp)
for lev,press_lev in enumerate(the_press):
    thetavals[0,lev,:,:]=theta(press_lev,the_temp[0,lev,:,:])   

# %% [markdown]
# ### Plot some profiles

# %%
fig,ax=plt.subplots(1,1)
for x,y in random_xy:
    theta_profile=thetavals[0,:,x,y]
    ax.plot(theta_profile,the_height*meter2km)  
ax.set(xlabel='potential temperature (K)',ylabel='height (km)',title='3 profiles');

# %% [markdown]
# ### Subtract average using numpy broadcasting

# %% [markdown]
# Find the perturbation wrt the mean theta at each height.  Because of Python's broadcasting
# rules ([broadcasting](http://docs.scipy.org/doc/numpy/user/basics.broadcasting.html)) we can't
# subtract a (1,112,110,100) and a (112,) array.  Numpy will let us subtract a
# (100,110,112) and a (112,) array  so we:
#
# 1) squeeze wvel to (112,110,100)
#
# 2) transpose to (100,110,112)
#
# 3) subtract the mean theta 
#
# 4) transpose back

# %%
w_bar=wvel[0,...].mean(axis=2).mean(axis=1)
wvel_squeeze=wvel.squeeze()
w_prime= wvel_squeeze.T - w_bar
w_prime=w_prime.T
print(w_prime.shape)

# %%
theta_bar=thetavals[0,...].mean(axis=2).mean(axis=1)
thetavals_squeeze=thetavals.squeeze()
theta_prime= thetavals_squeeze.T - theta_bar
theta_prime=theta_prime.T
print(theta_prime.shape)

# %%
flux=w_prime*theta_prime
print(flux.shape)
flux_profile=flux.mean(axis=2).mean(axis=1)
print(flux_profile.shape)

# %%
fig,ax=plt.subplots(1,1)
out=ax.plot(flux_profile,the_height*meter2km)
out=ax.set(xlabel='potential temperature flux (K m/s)',ylabel='height (km)',title='theta flux')

# %% [markdown]
# ### convert kinematic flux to W/m^2
#
# 1.  Calculate the the density profile using the ideal gas law
#
# $$p = \rho R_d T$$

# %%
Rd=287 #J/kg/K
rho=the_press/(Rd*xymean)
fig,ax=plt.subplots(1,1,figsize=(6,6))
out=ax.plot(rho,the_height*meter2km)
out=ax.set(xlabel=r'$\rho\ (kg/m^3)$',ylabel='height (km)',title='density of dry air')

# %% [markdown]
# 2\.  Calculate the enthalpy flux in $W m^{-2}$ using:
#
# $$\rho \overline{w^\prime h^\prime} = \rho c_{pd} \left( \overline{w^\prime T^\prime} \right ) \approx \rho c_{pd} \left( \overline{w^\prime \theta^\prime} \right ) $$
#
#
#
#

# %%
cpd=1004.  #J/kg/K
flux=cpd*rho*flux_profile
fig,ax=plt.subplots(1,1)
ax.plot(flux,the_height*meter2km)
ax.set(xlabel='sensible heat flux (W/m^2)',ylabel='height (km)',title='heat flux')
out=ax.set(ylim=(0,2))
