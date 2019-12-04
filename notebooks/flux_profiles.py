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
# # Introduction
#
# This notebook picks off the point at which the buoyancy flux crosses the 0 point -- i.e. the bottom of the entrainment layer

# %%
import urllib,os
from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
import a500
import context

from  a500.utils.data_read import download
the_root="http://clouds.eos.ubc.ca/~phil/docs/atsc500/data/"
the_file='case_60_10.nc'
out=download(the_file,root=the_root,dest_folder=a500.data_dir)



def make_theta(temp,press):
    """
      temp in K
      press in Pa
      returns theta in K
    """
    p0=1.e5
    Rd=287.  #J/kg/K
    cpd=1004.  #J/kg/K
    theta=temp*(p0/press)**(Rd/cpd)
    return theta


# %% [markdown]
# ### accumulate the vertical flux profile for each of the 10 ensemble members  into the avg_flux variable

# %%
case_name= a500.data_dir / 'case_60_10.nc'
#
#  look at the first ensemble member
#

with Dataset(str(case_name),'r') as ncin:
    ensemble='c1'
    group = ncin.groups['c1']
    the_shape = group.variables['TABS'].shape
    avg_flux=np.zeros(the_shape[:2],dtype=np.float64)
    grouplist=np.arange(1,11)
    groupnames=['c{:d}'.format(groupnum) for groupnum in grouplist]
    for count,the_group in enumerate(groupnames):
        print('working on {}'.format(the_group))
        group = ncin.groups[the_group]
        temp=group.variables['TABS'][...]
        wvel = group.variables['W'][...]
        horiz_mean_temp = temp.mean(axis=3).mean(axis=2)
        horiz_mean_wvel = wvel.mean(axis=3).mean(axis=2)
        temp = temp.T - horiz_mean_temp.T
        wvel = wvel.T - horiz_mean_wvel.T
        flux =(temp*wvel).T
        horiz_mean_flux = flux.mean(axis=3).mean(axis=2)
        avg_flux += horiz_mean_flux
    avg_flux = avg_flux/len(groupnames)
    press=ncin.variables['press'][...]
    z=ncin.variables['z'][...]
    the_time=ncin.variables['time'][...]

# %% [markdown]
# ### use searchsorted and argmin to get the first zero crossing for each of the 48 profiles

# %%
# %matplotlib inline
plt.style.use('ggplot')
fig,ax=plt.subplots(1,1,figsize=(10,8))
# cp = 1004.
rho = 1
zstore=[]
for i in np.arange(0,avg_flux.shape[0],1):
    #ax.plot(avg_flux[i,:],z,'ro')
    ax.plot(avg_flux[i,:]*rho*cp,z)
    negheight=np.argmin(avg_flux[i,:])
    zero_cross = np.searchsorted(-avg_flux[i,2:negheight],0.)
    z_trun = z[2:negheight]
    zstore.append(z_trun[zero_cross])
    ax.text(0,z_trun[zero_cross],'{:d}'.format(i))
    ax.plot(0,z_trun[zero_cross],'b+')
out=ax.set(xlabel=r'$\overline{\theta flux}$ ($W\,m^{-2}$)',ylabel='height (m)',
           title='LES dry run:  surface flux=60 $W\,m^{-2}$, $\Gamma$=10 K/km',
           ylim=[0,800])

dump_profs=True
if dump_profs:
    save_dict=dict(z=z,avg_flux=avg_flux,the_time=the_time)
    np.savez('flux_profs',**save_dict)



# %% [markdown]
# ### look at variablity among ensemble members at timestep 20

# %%
timestep = 20
hold_vals=np.empty([len(groupnames),the_shape[1]])
with Dataset(case_name,'r') as ncin:
    for count,the_group in enumerate(groupnames):
        group = ncin.groups[the_group]
        temp=group.variables['TABS'][timestep,...]
        wvel = group.variables['W'][timestep,...]
        horiz_mean_temp = temp.mean(axis=2).mean(axis=1)
        horiz_mean_wvel = wvel.mean(axis=2).mean(axis=1)
        temp = temp.T - horiz_mean_temp.T
        wvel = wvel.T - horiz_mean_wvel.T
        flux =(temp*wvel).T
        horiz_mean_flux = flux.mean(axis=2).mean(axis=1)
        hold_vals[count,:]=horiz_mean_flux

fig,ax=plt.subplots(1,1,figsize=(10,8))
for i in range(hold_vals.shape[0]):
    ax.plot(hold_vals[i,:],z)
out=ax.set(xlabel=r'$\overline{\theta}$ (K)',ylabel='height (m)',
       title='LES dry run:  surface flux=60 $W\,m^{-2}$, $\Gamma$=10 K/km')


# %% [markdown] {"collapsed": true}
# ### Does the inversion height grow as sqrt(time)?

# %%
fig,ax=plt.subplots(1,1,figsize=(10,8))
day_frac=the_time - 169
ax.loglog(day_frac,zstore,'ro',label='les zero-cross')
ax.loglog(day_frac,1000.*np.sqrt(day_frac),label='sqrt(time)')
ax.set(ylim=[100,1000],xlim=[5.e-3,1])
out=ax.legend()




# %% [markdown]
# # Save the results to a json file

# %%
# Dump this to a json file
import json
#
# turn float32 into float64 for json and remove mask
# by filling missing values
# json only understands lists, strings and dictionaries
# but dictionaries can be nested
#
missing = -999.
day_frac_filled = day_frac.filled(missing).astype(np.float64)
day_frac_filled = list(day_frac_filled)
zstore_64 = [float(item) for item in zstore]
keep_dict=dict(heights=zstore_64,times=day_frac_filled)
keep_dict['units']={"heights":"m","times":"fractions of a 24 hour day",
                    "surface_flux":"W/m^2","stability":"K/m"}
keep_dict['history']="written by flux_profiles.ipynb"
keep_dict['surface_flux']=60.
keep_dict['stability']=0.01
keep_dict['missing_value'] = missing
outfile = a500.data_dir / "dry_les_flux_profiles.json"
with open(outfile,'w') as output:
    json.dump(keep_dict,output,indent=4)

# %%
