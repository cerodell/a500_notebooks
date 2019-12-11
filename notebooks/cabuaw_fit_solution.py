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

# %% [markdown]
# # calculate monthly averages

# %%
monthly_wind_avg=jul_speed.mean(axis=(0,2))
monthly_wind_avg.shape

# %%
the_month='jul, 2014'
alpha=0.4
hour=2
fig,ax=plt.subplots(1,3,figsize=(8,6))
ax[0].plot(hourly_wind_avg[:,hour,:].T,z,alpha=alpha)
ax[0].plot(monthly_wind_avg[hour,:].T,z,'k-',lw=2)
ax[0].set(title='hour: {} UTC'.format(hour))
hour=10
ax[1].plot(hourly_wind_avg[:,hour,:].T,z,alpha=alpha)
ax[1].plot(monthly_wind_avg[hour,:].T,z,'k-',lw=2)
ax[1].set(title='hour: {} UTC'.format(hour))
fig.suptitle('{} hourly avg winds'.format(the_month))
hour=14
ax[2].plot(hourly_wind_avg[:,hour,:].T,z,alpha=alpha)
ax[2].plot(monthly_wind_avg[hour,:].T,z,'k-',lw=2)
ax[2].set(title='hour: {} UTC'.format(hour))
fig.suptitle('{} hourly avg winds with monthly average'.format(the_month));


# %% [markdown]
# # Fit the wind profiles to a modified log(z) function

# %% [markdown]
# Follow [Verkaik and Holtslag, 2007](http://ezproxy.library.ubc.ca/login?url=http://link.springer.com/10.1007/s10546-006-9121-1) and fit the windspeed to a modified log profile using scipy.optimize.curve_fit
#
# ```
# S=a0 + a1*z + a2*z**2 + a3*np.log(z)
# direc=b0 + b1*z + b2*z**2
# theta=c0 + c1*z + c2*z**2. + c3*np.log(z)
# ```

# %%
def wind_func(z, *coeffs):
    'nonlinear function in a and to fit to data'
    fullfit = coeffs[0] + coeffs[1]*z + coeffs[2]*z**2. + coeffs[3]*np.log(z)
    return fullfit

def wind_terms(z, *coeffs):
    'nonlinear function in a and to fit to data'
    fullfit = coeffs[0] + coeffs[1]*z + coeffs[2]*z**2. + coeffs[3]*np.log(z)
    term1 = coeffs[0] + coeffs[3]*np.log(z)
    term12 = term1 + coeffs[1]*z
    term123 = term12 + coeffs[2]*z**2.
    return fullfit, term1, term12, term123


# %%
import numpy as np
from scipy.optimize import curve_fit
import scipy
from scipy.stats.distributions import  t
import warnings
warnings.filterwarnings("ignore",category=scipy.optimize.OptimizeWarning)

# %% [markdown]
# ## 2am UTC   -- most stable

# %%
# 2 UTC
hour1=2 # early morning
fig,ax=plt.subplots(1,2,figsize=(8,6))
fig.suptitle(f'July 2014 U profiles at {hour1} UTC')
ax[0].plot(monthly_wind_avg[hour1,:],z)
ax[0].set(title='monthly average')
#
# flip tower data so it goes from bottom to top
# and get rid of the lowest level, which doesn't
# have a measurement
#
rev_z=z[::-1]
rev_z=rev_z[1:]
#test=hourly_wind_avg[day,hour1,::-1]
test=monthly_wind_avg[hour1,::-1]
#
# lose the bottom value
#
test=test[1:]
initial_guess=[0,0,0,0]
pars, pcov = curve_fit(wind_func, rev_z, test, p0=initial_guess)
zinterp=np.linspace(10.,200,50.)
speed_interp=wind_func(zinterp,*pars)
ax[0].plot(speed_interp,zinterp,'ro',alpha=0.5)
print(f'fit coefficients {pars}')
fit, term1, term12, term123 = wind_terms(zinterp, *pars)
ax[1].plot(fit,zinterp,'ro',alpha=0.8,label='fullfit')
ax[1].plot(term1,zinterp,label='term1')
ax[1].plot(term12,zinterp,label='term12')
ax[1].plot(term123,zinterp,label='term123')
ax[1].legend();

# %% [markdown]
# ## 10am  -- convection starting

# %%
# 10 UTC
hour2=10 # later morning
fig,ax=plt.subplots(1,3,figsize=(10,6))
fig.suptitle(f'U profiles for  at {hour2} UTC')

test=monthly_wind_avg[hour2,::-1]
#
# lose the bottom value
#
test=test[1:]
initial_guess=[0,0,0,0]
pars, pcov = curve_fit(wind_func, rev_z, test, p0=initial_guess)
zinterp=np.linspace(10.,200,50.)
speed_interp=wind_func(zinterp,*pars)
ax[0].plot(speed_interp,zinterp,'ro',alpha=0.5)
print(f'fit coefficients {pars}')
fit, term1, term12, term123 = wind_terms(zinterp, *pars)
ax[1].plot(fit,zinterp,'ro',alpha=0.8,label='fullfit')
ax[1].plot(term1,zinterp,label='term1')
ax[1].plot(term12,zinterp,label='term12')
ax[1].plot(term123,zinterp,label='term123')
ax[1].legend();

# %%
# 14 UTC
hour3=14 # afternoon
fig,ax=plt.subplots(1,2,figsize=(10,6))
fig.suptitle(f'U profiles July at {hour3} UTC')
test=monthly_wind_avg[hour3,::-1]
#
# lose the bottom value
#
test=test[1:]
initial_guess=[0,0,0,0]
pars, pcov = curve_fit(wind_func, rev_z, test, p0=initial_guess)
zinterp=np.linspace(10.,200,50.)
speed_interp=wind_func(zinterp,*pars)
ax[0].plot(speed_interp,zinterp,'ro',alpha=0.5)
print(f'fit coefficients {pars}')
fit, term1, term12, term123 = wind_terms(zinterp, *pars)
ax[1].plot(fit,zinterp,'ro',alpha=0.8,label='fullfit')
ax[1].plot(term1,zinterp,label='term1')
ax[1].plot(term12,zinterp,label='term12')
ax[1].plot(term123,zinterp,label='term123')
ax[1].legend();


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

# %% [markdown]
# # check the u profile at 20 and 80 meters
#
# How good do equations 18 and 19 of our similarity layer notes do in
# fitting the Cabauw profiles?  The functions `find_gradient` and
# `find_scaling` calculate the actual scaled gradient phi from
# the tower profiles.  Remember the conclusion in Figure 5 of the V-H paper.  The hedgerows,etc. cause more mixing, which decreases the gradient in stable conditions compared to completely flat terrain.
#
# <img src="images/figure5_vh.png" width="60%" />
#
# Can we reproduce these results with their dataset?
#
# Here are equations 18 and 19 for reference:

# %% [markdown]
# ## the measurements

# %%
# get indexes for z = 20 m and z = 80 m
level_80 = np.where(z == 80)[0][0]
level_20 = np.where(z == 20)[0][0]

# helper functions for calculating phis
def find_gradient(hourly_wind_avg, z, day, hour, level):
    "Function to quickly calculate dUdz for a certain hour and level"
    "Uses 'forward euler' to numerically calculate gradient"
    
    return (hourly_wind_avg[day,hour,level] - 
            hourly_wind_avg[day,hour,level+1])/(z[level] - z[level+1])

def find_scaling(hourly_wind_avg, z, L, ustars, day, hour, level):
    """Calculate phi and zeta at a given level and time"""
    
    k = 0.40
    height = z[level]
    ustar = ustars[day, hour, level]
    dUdz = find_gradient(hourly_wind_avg, z, day, hour, level)
    phi = (k*height/ustar)*dUdz
    zeta = z[level]/L[day, hour, level]
    return phi,zeta


# %% [markdown]
# ## The B-Dyer fits
#
# The two functions `dyer` and `fit_dyer` calculate the phi predicted from the
# textbook equations and then using the V-Holtzlag fitted coefficients

# %%
def dyer(zeta):
    'calculate curve from Dyer equations (Surface layer similarity 18 + 19)'
    #gradient_fit = (zeta**(1))*(coeffs[0] + 2*coeffs[1]*zeta + coeffs[2]*1/zeta) # derivative, dU/dz as estimated from Nieuwstadt (1984)
    gamma_1 =16
    beta = 5
    if zeta < 0:
        # unstable
        phi_dyer = (1 - gamma_1*zeta)**(-1/4)
    else:
        phi_dyer = 1 + beta*zeta
    return phi_dyer

def fit_dyer(zetas, *coeffs):
    'calculate curve from Dyer equations (Surface layer similarity 18 + 19)'
    #gradient_fit = (zeta**(1))*(coeffs[0] + 2*coeffs[1]*zeta + coeffs[2]*1/zeta) # derivative, dU/dz as estimated from Nieuwstadt (1984)
    phi_dyers = []
    for zeta in zetas:
        if zeta < 0:
            # unstable
            phi_dyer = (1 - coeffs[0]*zeta)**(-1/4)
        else:
            phi_dyer = 1 + coeffs[1]*zeta
        phi_dyers.append(phi_dyer)
    return phi_dyers


# %% [markdown]
# # Plot these
#
# Compare the B-dyer equations with the V-H local scaling fit and the Cabauw observations
#
# Note that the 

# %%
# calculate and plot cabauw data with the Dyer equations + fit
level_80 = np.where(z == 80)[0][0]
level_20 = np.where(z == 20)[0][0]
day = 19

hour1, hour2, hour3 = 2, 10, 14
hours = [hour1,hour2,hour3]
levels = [level_80, level_20]

phi_len = len(hours)*len(levels)
phis = np.empty([phi_len])
zetas = np.empty([phi_len])

phis80 = []
zetas80 = []
phis20 = []
zetas20 = []

# calculate phi from data
for hour in hours:
    for level in levels:
        phi,zeta = find_scaling(hourly_wind_avg, z, L, USTAR, day, hour, level)
        if level==levels[0]:
            phis80.append(phi)
            zetas80.append(zeta)
        else:
            phis20.append(phi)
            zetas20.append(zeta)

fig,ax=plt.subplots(1,2,figsize=(10,6))
fig.suptitle(r"$\phi_m\ vs.\ \zeta\ for\ day\ {}$".format(day,hour3),fontsize=20)
ax[0].plot(zetas80, phis80, '.b',label='cabauw')
ax[1].plot(zetas20, phis20, '.g',label='cabauw')
ax[0].set(xlabel=r'$\zeta$',ylabel=r'$\phi_m$',title='z = {} m'.format(z[levels[0]]))
ax[1].set(xlabel=r'$\zeta$',ylabel=r'$\phi_m$',title='z = {} m'.format(z[levels[1]]))

# calculate Dyer values for phi
dyer_zetas80 = np.arange(min(zetas80),max(zetas80),0.001)
dyer_zetas20 = np.arange(min(zetas20),max(zetas80),0.001)

phis_dyer80 = [dyer(dyer_zeta80) for dyer_zeta80 in dyer_zetas80]   
phis_dyer20 = [dyer(dyer_zeta20) for dyer_zeta20 in dyer_zetas20] 

ax[0].plot(dyer_zetas80, phis_dyer80 , '-m',label='dyer')
ax[1].plot(dyer_zetas20, phis_dyer20, '-m',label='dyer')

# fit curve for z = 80 m
initial_guess=[0,0]
pars80, pcov80 = curve_fit(fit_dyer, zetas80, phis80, p0=initial_guess)
zetainterp=np.linspace(min(zetas80),max(zetas80),50.)
grad_interp=fit_dyer(zetainterp,*pars80)
fit80=ax[0].plot(zetainterp,grad_interp,'c--',alpha=0.5,label='fit')

# fit curve for z = 20 m
initial_guess=[0,0]
pars20, pcov20 = curve_fit(fit_dyer, zetas20, phis20, p0=initial_guess)
zetainterp=np.linspace(min(zetas20),max(zetas20),50.)
grad_interp=fit_dyer(zetainterp,*pars20)
fit20=ax[1].plot(zetainterp,grad_interp,'c--',alpha=0.5,label='fit')

ax[0].legend(loc='best')
ax[1].legend(loc='best')


# %%
