# -*- coding: utf-8 -*-
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
# <div class="toc"><ul class="toc-item"><li><span><a href="#use-interp1d-to-find-the-zero-crossing--for-each-of-the-48-profiles" data-toc-modified-id="use-interp1d-to-find-the-zero-crossing--for-each-of-the-48-profiles-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>use interp1d to find the zero crossing  for each of the 48 profiles</a></span></li><li><span><a href="#look-at-variablity-among-ensemble-members-at-timestep-20" data-toc-modified-id="look-at-variablity-among-ensemble-members-at-timestep-20-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>look at variablity among ensemble members at timestep 20</a></span></li><li><span><a href="#Does-the-inversion-height-grow-as-sqrt(time)?" data-toc-modified-id="Does-the-inversion-height-grow-as-sqrt(time)?-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Does the inversion height grow as sqrt(time)?</a></span></li><li><span><a href="#How-much-does-entrainment-accelerate-boundary-layer-growth?" data-toc-modified-id="How-much-does-entrainment-accelerate-boundary-layer-growth?-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>How much does entrainment accelerate boundary layer growth?</a></span></li></ul></div>

# %%
import urllib,os
from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np

the_vars=np.load('flux_profs.npz')
avg_flux = the_vars['avg_flux']
z= the_vars['z']
the_time=the_vars['the_time']

# %%
from netCDF4 import Dataset
with Dataset('case_60_10.nc','r') as f:
    print(f.gamma_K_km)

# %% [markdown]
# ### use interp1d to find the zero crossing  for each of the 48 profiles

# %%
# %matplotlib inline
from scipy.interpolate import interp1d
plt.style.use('ggplot')
fig,ax=plt.subplots(1,1,figsize=(10,8))
zstore=[]
for i in np.arange(0,avg_flux.shape[0],1):
    negheight=np.argmin(avg_flux[i,:])
    z_trun = z[2:negheight]
    flux=avg_flux[i,2:negheight]
    ax.plot(flux,z_trun)
    the_spline= interp1d(flux,z_trun)
    zero_cross=the_spline(0.)
    zstore.append(zero_cross)
    ax.text(0,zero_cross,'{:d}'.format(i))
    ax.plot(0,zero_cross,'b+')
   
    
out=ax.set(xlabel=r'$\overline{\theta}$ (K)',ylabel='height (m)',
           title='LES dry run:  surface flux=60 $W\,m^{-2}$, $\Gamma$=10 K/km',
           ylim=[0,800])


# %% [markdown]
# ### look at variablity among ensemble members at timestep 20

# %% [markdown] {"collapsed": true}
# ### Does the inversion height grow as sqrt(time)?

# %% [markdown]
# Since the LES was run with a constant surface flux of $F_0$ = 60 $W\,m^{-2}$ the inversion
# height should follow Stull 11.2.2f, p. 456:
#
# $$z_i^2 - z_{i0}^2 = \frac{2}{\gamma} \left [ \overline{w^\prime \theta^\prime}_s 
# - \overline{w^\prime \theta^\prime}_i \right ] \left ( t - t_0 \right ) $$
#
# With $\overline{w^\prime \theta^\prime}_s = F_0/c_p$ = 60/1004. = 0.06 $K\,m\,s^{-1}$
#
# So plotting log(height) vs. log(time) should yield a slope of 1/2.  The actual result is a slope
# of 0.52 (see below)

# %%
fig,ax=plt.subplots(1,1,figsize=(10,8))
day_frac=the_time - 169
ax.loglog(day_frac,zstore,'ro',label='les zero-cross')
ax.set(ylim=[100,1000],xlim=[5.e-3,1])

xdata=np.log10(day_frac)
ydata=np.log10(zstore)
slope, intercept=np.polyfit(xdata,ydata,1)
tvals=np.linspace(-2,np.log10(0.3),20)
yvals= intercept + slope*tvals
ax.loglog(10**tvals,10**yvals,'b-',label='fit')
out=ax.legend()
print("depth growing as t**{:<4.2f}".format(slope))

# %% [markdown]
# ### How much does entrainment accelerate boundary layer growth?

# %% [markdown]
# The turbulent plumes are doing work on the inversion, causing warm air to be
# pulled into the mixed layer and accelerating the rate at which the boundary
# layer is growing.  
#
# According Garratt eq. 6.20 if $\beta = \frac{ \overline{w^\prime \theta^\prime}_i}{\overline{w^\prime \theta^\prime}_s}$, then 
#
# $$\frac{\partial}{\partial t} \frac{h^2}{2} = \gamma^{-1} \left (1 + 2 \beta \right )\overline{w^\prime \theta^\prime}_s$$
#
#

# %%
fig,ax=plt.subplots(1,1,figsize=(10,8))
day_frac=the_time - 169
day_frac_sec = day_frac*24.*3600.  #time in seconds
h2 = np.array(zstore)**2.
ax.plot(day_frac,h2*1.e-6,'ro',label='les zero-cross')
slope, intercept=np.polyfit(day_frac_sec,h2,1)
cp=1004.  #J/kg/K
F0=60 # W/m^2
F0=F0/1004.  #convert to kinematic units assuming density=1 kg/m^3
gamma = 1.e-2  #K/m
#
# convert time from days to seconds
#
tvals=np.linspace(0,0.35,20)*24.*3600.
yvals= intercept + slope*tvals
ax.plot(tvals/24./3600.,yvals*1.e-6,'b-',label='fit')
ax.set(ylabel=r"squared inversion height $(m^{2} \times 10^{-6})$",xlabel="time (days)")
out=ax.legend(loc='best')

# %%
print(f'slope {slope:5.2f} m**2/second')

# %% [markdown]
# Solve for $\beta$:

# %%
F0=60./1004.
beta2=slope*gamma/(2.*F0) - 1.
print('β = {:<5.2f}'.format(beta2/2.))

# %% [markdown]
# So the LES is entraining more efficiently than predicted, probably due to the 25 meter resolution.

# %%
