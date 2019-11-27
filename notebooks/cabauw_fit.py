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

# %% [markdown] {"toc": true}
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#2am-UTC" data-toc-modified-id="2am-UTC-1">2am UTC</a></span></li><li><span><a href="#Q1-2:00-UTC----fill-in-the-the-cell-above-to-plot-the-curve-fit-on-top-of-the-hourly-average" data-toc-modified-id="Q1-2:00-UTC----fill-in-the-the-cell-above-to-plot-the-curve-fit-on-top-of-the-hourly-average-2">Q1 2:00 UTC -- fill in the the cell above to plot the curve fit on top of the hourly average</a></span></li><li><span><a href="#Q2--repeat-this-for-10:00--UTC-above" data-toc-modified-id="Q2--repeat-this-for-10:00--UTC-above-3">Q2  repeat this for 10:00  UTC above</a></span></li><li><span><a href="#Q3:-repeat-this-for-14:00--UTC-above" data-toc-modified-id="Q3:-repeat-this-for-14:00--UTC-above-4">Q3: repeat this for 14:00  UTC above</a></span></li><li><span><a href="#Buoyancy-flux-and-L" data-toc-modified-id="Buoyancy-flux-and-L-5">Buoyancy flux and L</a></span></li></ul></div>

# %% [markdown]
# # Replace my month with yours and fill in the curve fits for 2am, 10am and 2pm

# %%
import glob
from netCDF4 import Dataset
from dateutil.parser import parse
import datetime
import numpy as np
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pytz
from datetime import datetime as dt

# %% [markdown]
# # read in the july 2014 profiles

# %%
the_file='cabauw_ubc.nc'
group='m201407'
with Dataset(the_file,'r') as nc_ubc:
    jul_nc=nc_ubc.groups[group]
    z=nc_ubc.variables['z'][...]
    jul_speed=jul_nc.variables['F'][...]
    jul_ta002 = jul_nc.variables['TA002']

jul_speed.shape

# %% [markdown]
# # calculate hourly averages

# %%
hourly_wind_avg=jul_speed.mean(axis=2)

# %%
hourly_wind_avg.shape

# %%
the_month='jul, 2014'
hour=2
fig,ax=plt.subplots(1,3,figsize=(8,6))
ax[0].plot(hourly_wind_avg[:,hour,:].T,z)
ax[0].set(title='hour: {} UTC'.format(hour))
hour=10
ax[1].plot(hourly_wind_avg[:,hour,:].T,z)
ax[1].set(title='hour: {} UTC'.format(hour))
fig.suptitle('{} hourly avg winds'.format(the_month))
hour=14
ax[2].plot(hourly_wind_avg[:,hour,:].T,z)
ax[2].set(title='hour: {} UTC'.format(hour))
fig.suptitle('{} hourly avg winds'.format(the_month));


# %% [markdown]
# # Fit the wind profiles to a modified log(z) function 

# %% [markdown]
# Follow [Verkaik and Holtslag, 2007](http://ezproxy.library.ubc.ca/login?url=http://link.springer.com/10.1007/s10546-006-9121-1) and fit the windspeed to a modified log profile using scipy.optimize.curve_fit (see their page 710 below equation 1)
#
# ```
# S=a0 + a1*z + a2*z**2 + a3*np.log(z)
# direc=b0 + b1*z + b2*z**2
# theta=c0 + c1*z + c2*z**2. + c3*np.log(z)
# ```

# %%
def wind_func(z, *coeffs):
    'nonlinear function in a and to fit to data'
    fit = coeffs[0] + coeffs[1]*z + coeffs[2]*z**2. + coeffs[3]*np.log(z)
    return fit


# %%
import numpy as np
from scipy.optimize import curve_fit
import scipy
from scipy.stats.distributions import  t
import warnings
warnings.filterwarnings("ignore",category=scipy.optimize.OptimizeWarning)

# %% [markdown]
# ## 2am UTC

# %%
# 2 UTC
hour1=2 # early morning
day=19 
fig,ax=plt.subplots(1,3,figsize=(10,6))
sample1=jul_speed[day,hour,:,:]
fig.suptitle('U profiles for day {} at {} UTC'.format(day,hour1))
ax[0].plot(sample1.T,z)
ax[0].set(title='6 10 minute samples',xlabel='wind speed (m/s)',
         ylabel='height (m)')
ax[1].plot(hourly_wind_avg[day,hour1,:],z)
ax[1].set(title='hourly average')
ax[2].plot(hourly_wind_avg[day,hour1,:],z)
ax[2].set(title='hourly average plus interpolated values')
#
# flip tower data so it goes from bottom to top
# and get rid of the lowest level, which doesn't
# have a measurement
#
rev_z=z[::-1]
rev_z=rev_z[1:]
test=hourly_wind_avg[day,hour1,::-1]
#
# lose the bottom value
#
test=test[1:]

# %% [markdown]
# ## Q1 2:00 UTC -- fill in the the cell above to plot the curve fit on top of the hourly average

# %%
#
# use curve fit to find the interpolated wind speed and plot it
#
#ax[2].plot(speed_interp,zinterp,'ro',alpha=0.5)
#print(f'fit coefficients {pars}')

# %%
# 10 UTC
hour2=10 # later morning
fig,ax=plt.subplots(1,3,figsize=(10,6))
fig.suptitle('U profiles for day {} at {} UTC'.format(day,hour2))
sample2=jul_speed[day,hour2,:,:]
ax[0].plot(sample2.T,z)
ax[0].set(title='6 10 minute samples',xlabel='wind speed (m/s)',
         ylabel='height (m)')
ax[1].plot(hourly_wind_avg[day,hour2,:],z)
ax[1].set(title='hourly average')
ax[2].plot(hourly_wind_avg[day,hour2,:],z)
ax[2].set(title='hourly average plus interpolated values')
#
# flip tower data so it goes from bottom to top
# and get rid of the lowest level, which doesn't
# have a measurement
#
rev_z=z[::-1]
rev_z=rev_z[1:]
test=hourly_wind_avg[day,hour2,::-1]
test=test[1:]



# %% [markdown]
# ## Q2  repeat this for 10:00  UTC above

# %%
# 14 UTC
hour3=14 # afternoon
fig,ax=plt.subplots(1,3,figsize=(10,6))
fig.suptitle('U profiles for day {} at {} UTC'.format(day,hour3))
sample3=jul_speed[day,hour3,:,:]
ax[0].plot(sample3.T,z)
ax[0].set(title='6 10 minute samples',xlabel='wind speed (m/s)',
         ylabel='height (m)')
ax[1].plot(hourly_wind_avg[day,hour3,:],z)
ax[1].set(title='hourly average')
ax[2].plot(hourly_wind_avg[day,hour3,:],z)
ax[2].set(title='hourly average plus interpolated values')
#
# flip tower data so it goes from bottom to top
# and get rid of the lowest level, which doesn't
# have a measurement
#
rev_z=z[::-1]
rev_z=rev_z[1:]

test=hourly_wind_avg[day,hour3,::-1]
test=test[1:]

# %% [markdown] {"trusted": true}
# ## Q3: repeat this for 14:00  UTC above

# %% [markdown]
# # calculate L

# %%
with Dataset(the_file,'r') as nc_ubc:
    jul_nc=nc_ubc.groups[group]
    H=jul_nc.variables['H'][...]
    LE = jul_nc.variables['LE'][...]
    USTAR = jul_nc.variables['UST'][...]
    TA002 = jul_nc.variables['TA002'][...]
    Q002 = jul_nc.variables['Q002'][...]
    P0 = jul_nc.variables['P0'][...]
    timevec = jul_nc.variables['time'][...]
    timevec = [dt.fromtimestamp(item,pytz.utc) \
               for item in timevec.flat]
    
Rd=287.  #J/kg/K
# cp = 1004.  #J/kg/K
k = 0.4
g=9.8
rho = P0*1.e2/(Rd*(TA002 + 273.15))



# %% [markdown]
# ## Buoyancy flux and L

# %%
#fleagle and bussinger eq. 6.31
Eb = H + 0.02*LE
#virtural temperature 
Tv = TA002 + 273.15  + 0.61*Q002*1.e-3
#Fleagle and Businger 6.47
L = - Tv*cp*rho*USTAR**3./(k*g*Eb)
good = np.abs(Eb) > 1

# %%
fig,ax=plt.subplots(1,1)
out=plt.hist(L[good].flatten(),bins=np.linspace(-150,150,40))

# %%
fig,ax=plt.subplots(1,1,figsize=(8,6))
fig.autofmt_xdate()
ax.plot(timevec,L.flatten())
title='Obukhov length L {}'.format(the_month)
out=ax.set(title=title,ylabel='L $(m)$',ylim=[-150,150])

# %%
