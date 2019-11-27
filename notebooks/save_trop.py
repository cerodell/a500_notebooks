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
#   latex_metadata:
#     affiliation: University of British Columbia
#     author: Philip Austin
#     bib: notebook.bib
#     course: ATSC 500
#     title: Vapor transport budget
# ---

# %% [markdown] {"toc": true}
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#Reading-a-netcdf-file" data-toc-modified-id="Reading-a-netcdf-file-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Reading a netcdf file</a></span></li><li><span><a href="#By-default-netCDF4-converts-netcdf-variables-to-masked-arrays" data-toc-modified-id="By-default-netCDF4-converts-netcdf-variables-to-masked-arrays-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>By default netCDF4 converts netcdf variables to masked arrays</a></span><ul class="toc-item"><li><span><a href="#How-much-liquid-water-is-in-the-domain?" data-toc-modified-id="How-much-liquid-water-is-in-the-domain?-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>How much liquid water is in the domain?</a></span></li><li><span><a href="#As-expected,-vapor-transport-dominates-the-energy-flux-in-the-warm-marine-boundary-layer" data-toc-modified-id="As-expected,-vapor-transport-dominates-the-energy-flux-in-the-warm-marine-boundary-layer-2.2"><span class="toc-item-num">2.2&nbsp;&nbsp;</span>As expected, vapor transport dominates the energy flux in the warm marine boundary layer</a></span></li></ul></li><li><span><a href="#Now-look-at-the-individual-budget-terms" data-toc-modified-id="Now-look-at-the-individual-budget-terms-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Now look at the individual budget terms</a></span></li><li><span><a href="#Save-all-this-in-a-dictionary" data-toc-modified-id="Save-all-this-in-a-dictionary-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>Save all this in a dictionary</a></span></li><li><span><a href="#check-some-values" data-toc-modified-id="check-some-values-5"><span class="toc-item-num">5&nbsp;&nbsp;</span>check some values</a></span></li><li><span><a href="#Do-a-checkpoint-save-for-all-the-variables" data-toc-modified-id="Do-a-checkpoint-save-for-all-the-variables-6"><span class="toc-item-num">6&nbsp;&nbsp;</span>Do a checkpoint save for all the variables</a></span></li><li><span><a href="#What-are-the-rms-perturbation-profiles?" data-toc-modified-id="What-are-the-rms-perturbation-profiles?-7"><span class="toc-item-num">7&nbsp;&nbsp;</span>What are the rms perturbation profiles?</a></span></li><li><span><a href="#set-scaling-parameters-q*-and-w*-from-rms-values" data-toc-modified-id="set-scaling-parameters-q*-and-w*-from-rms-values-8"><span class="toc-item-num">8&nbsp;&nbsp;</span>set scaling parameters q* and w* from rms values</a></span></li><li><span><a href="#Look-at-mean-gradient-terms-in-XI" data-toc-modified-id="Look-at-mean-gradient-terms-in-XI-9"><span class="toc-item-num">9&nbsp;&nbsp;</span>Look at mean gradient terms in XI</a></span></li><li><span><a href="#Plot-the-three-equations-for-term-XI" data-toc-modified-id="Plot-the-three-equations-for-term-XI-10"><span class="toc-item-num">10&nbsp;&nbsp;</span>Plot the three equations for term XI</a></span></li></ul></div>

# %% [raw]
# \pagestyle{first}

# %% [markdown]
# # Including water vapor in the tropical analysis

# %% [markdown]
# 1.  What do the fluxes of sensible and latent heat look like for the tropical case?
#
# 1.  Drill down into the data using
#
#     a.  Two dimensional histograms
#     
#     b.  images with normalized color pallettes

# %% [markdown]
# ## Reading a netcdf file

# %% [markdown]
# 1.  Fetch the netcdf file tropical.nc from my webserver (50 Mbytes)

# %%
import context
import urllib.request
import numpy as np
from netCDF4 import Dataset
from matplotlib import pyplot as plt
from a500.utils.ncdump import ncdump
import a500
from a500.thermo import calc_theta
import copy

# %%
from a500.utils.data_read import download
the_root = "http://clouds.eos.ubc.ca/~phil/docs/atsc500"
the_file =  "tropical_subset.nc"
out=download(the_file,root=the_root,dest_folder=a500.data_dir)
g2kg=1.e-3
kg2g=1.e3


# %%
def calc_thetav(theta, wv, wl):
    """
    Thetav using stull 1.5.1a
    
    Parameters
    ----------
    
    theta: float
       potential temperature(K)

    wv: float
       vapor mixing ratio (kg/kg)
       
    wl: float
       liquid mixing ratio (kg/kg)
       
    Returns
    -------
    
    theta_v: float
      virtual potential temperature
    """
    theta_v = theta*(1 + 0.61*wv - wl)
    return theta_v


# %%
def xy_avg(var):
    """
    assuming var is 3d with axis (z,y,x)
    take the avg over x followed by the axg over y
    """
    return array3d.mean(axis=2).mean(axis=1)

def do_reynolds(array3d):
    """
        do a spatial-mean reynolds average of a 3d field array3d
        needs dimensions arranged as (z,y,x)
        returns avg(z),perturb(z,y,x)
    """
    avg=array3d.mean(axis=2).mean(axis=1)
    perturb=array3d.T - avg
    perturb=perturb.T
    return avg,perturb


# %% [markdown]
# 2\.  Identify the file structure using ncdump

# %%
with Dataset(a500.data_dir / 'tropical_subset.nc') as ncin:
    ncdump(ncin)

# %% [markdown]
# 3\.  Read a variable using the netCDF4 module:  http://unidata.github.io/netcdf4-python/

# %%
with Dataset(a500.data_dir / 'tropical_subset.nc','r') as nc_in:
    print(list(nc_in.variables.keys()))
    the_height=nc_in.variables['z'][...]
    hit = the_height < 2000.
    zvals = the_height[hit]
    xvals=nc_in.variables['x'][...]
    yvals=nc_in.variables['y'][...]
    the_temp=nc_in.variables['TABS'][0,hit,:,:]  
    #
    # remove the time dimension since we only have one timestep
    #   
    print('temp shape',the_temp.shape)
    xvals=nc_in.variables['x'][...]
    yvals=nc_in.variables['y'][...]
    print('height shape',zvals.shape)
    the_press=nc_in.variables['p'][hit]
    the_press=the_press*100.  #convert to Pa
    wvel=nc_in.variables['W'][0,hit,:,:]  #m/s
    uvel=nc_in.variables['U'][0,hit,:,:]  #m/s
    vvel=nc_in.variables['V'][0,hit,:,:]  #m/s
    qv=nc_in.variables['QV'][0,hit,:,:]*g2kg  #vapor kg/kg
    ql=nc_in.variables['QN'][0,hit,:,:]*g2kg  #liquid  kg/kg
    pp = nc_in.variables['PP'][0,hit,:,:]  #Pressure perturbation (Pa)
    theta = (calc_theta(the_temp.T,the_press)).T
    thetav = calc_thetav(theta,qv,ql)


# %% [markdown]
# ## By default netCDF4 converts netcdf variables to masked arrays

# %%
type(ql), type(qv), type(wvel)

# %% [markdown]
# See [the masked array docs](https://docs.scipy.org/doc/numpy/reference/maskedarray.generic.html)

# %% [markdown]
# ### How much liquid water is in the domain?

# %%
bins=np.linspace(0.05,1.5,20)
out=plt.hist(ql.compressed()*kg2g,bins=bins)
plt.title('histogram of cloud lwc (g/kg)');

# %% [markdown]
# # Calculate the temperature and vapor fluxes

# %%
Rd=287  #J/kg/K
g2kg = 1.e-3
kg2g= 1.e3
cpd=1004.
lv = 2.5e6
meter2km = 1.e-3
qv_avg,qv_perturb=do_reynolds(qv)
temp_avg,temp_perturb=do_reynolds(the_temp)
rho_avg=the_press/(Rd*temp_avg)
w_avg,w_perturb = do_reynolds(wvel)
T_flux=((w_perturb*temp_perturb).T*rho_avg).T*cpd  #W/m^2
qv_flux=((w_perturb*qv_perturb).T*rho_avg).T*lv #W/m^2
T_flux_1d=(T_flux).mean(axis=2).mean(axis=1)
qv_flux_1d=(qv_flux).mean(axis=2).mean(axis=1)

# %%
fig,ax=plt.subplots(1,2,figsize=(10,6))
ax[0].plot(T_flux_1d,zvals*meter2km)  
ax[0].set(xlabel='sensible heat flux (W/m^2)',ylabel='height (km)',title='temperature flux')
ax[1].plot(qv_flux_1d,zvals*meter2km)  
out=ax[1].set(xlabel='latent heat flux (W/m^2)',ylabel='height (km)',title='vapor flux')
out=[the_ax.set(ylim=(0,1.5)) for the_ax in ax]

# %% [markdown]
# ### As expected, vapor transport dominates the energy flux in the warm marine boundary layer

# %% [markdown]
# ## Now look at the individual budget terms

# %%
theta = calc_theta(the_temp.T,the_press)
theta = np.array(theta.T)
keys=['u','v','w','pp','theta','qv','thetav']
values = [uvel,vvel,wvel,pp,theta,qv,thetav]
var_dict=dict(zip(keys,values))
var_dict['p0']=the_press
var_dict['z'] = zvals
var_dict['x'] = xvals
var_dict['y'] = yvals


# %%
def save_reynolds(key,var):
    """
    given a variable name (key), and the variable
    make two new keys "key_avg" and "key_pert"
    and return the reynolds average and perturbation
    as a dictionary with those keys
    
    Parameters
    ----------
    
    key: string
       name of variable
       
    var: ndarray or masked array
       data array
    """
    new_keys=[f"{key}_avg",f"{key}_pert"]
    avg, perturb = do_reynolds(var)
    return dict(zip(new_keys,(avg,perturb)))


# %% [markdown]
# ## Save all this in a dictionary

# %%
keep_dict = copy.copy(var_dict)
for key, value in var_dict.items():
    #
    # pressure perturbation already done
    #
    if key in ['pp','x','y','z','p0']:
        continue
    new_dict = save_reynolds(key, value)
    keep_dict.update(new_dict)
keep_dict['qv_pert_grams']=keep_dict['qv_pert']*1.e3
print(keep_dict.keys())


# %% [markdown]
# ## check some values  
#
# loop through the dictionary and make some histograms

# %%
def hist_prep(var):
    """
    compress and unflatten arrays for histogramming
    """
    try:
        #
        # masked array
        #
        var=var.compressed().flat
    except AttributeError:
        #
        # regular ndarray
        #
        var = var.flat
    return var
    


# %%
fig, axvec = plt.subplots(4,2,figsize=(10,10),constrained_layout=True)
the_vars = ['u_pert', 'v_pert', 'w_pert', 'theta_pert','thetav_pert','qv_pert_grams','pp']
axvec_flat = list(axvec.flat)
for count,key in enumerate(the_vars):
    var=hist_prep(keep_dict[key])
    axvec_flat[count].hist(var)
    axvec_flat[count].set_title(key)
fig.delaxes(axvec_flat[-1])

# %% [markdown]
# ## Do a checkpoint save for all the variables

# %%
the_file = a500.data_dir / 'tropical_vapor.npz'
np.savez_compressed(the_file,**keep_dict)

# %%
list(keep_dict.keys())

# %%
a=np.load(the_file)


# %% [markdown]
# \begin{align}
# \frac{\partial(\overline{q^{\prime} u^\prime_i})}{\partial t}+
# \underbrace{\bar{U_{j}} 
# \frac{\partial(\overline{q^{\prime} u^\prime_{i}})}{\partial x_{j}}}_{II}
# =-\underbrace{\overline{q^\prime u^\prime_j} \frac{\partial \bar{U}_{i}}{\partial x_{j}}}_{III}
# -\underbrace{\overline{u^\prime_{i} u^\prime_{j}} \frac{\partial \bar{q}}{\partial x_{j}}}_{XI}
# + \underbrace{\frac{\partial ( \overline{ q^{\prime} u^\prime_j u^\prime_i } ) }{\partial x_j}}_{IV}
# \end{align}
#
# %
#
# %% [markdown]
# \begin{align}
# +\underbrace{\delta_{i 3}\left( \frac{\overline{\mathrm{q}^{\prime} \theta_{v}^{\prime}}}{\bar{\theta}_{v}}\right ) g}_V
# + \underbrace{f_c \epsilon_{i j 3} \left ( \overline{u_j^\prime q^\prime} \right )}_{VI}
# -  \left ( \frac{1}{\bar{\rho}}\right )
# \left [ \underbrace { \frac{ \partial ( \overline{p^\prime q^\prime} ) } {\partial x_i}}_{VII} 
# + \underbrace{ \overline{ p^\prime \frac{\partial q^\prime}{\partial x_i}} }_{VIII} \right ]
# %\right ] }
# \end{align}
# %% [markdown]
# \begin{align}
# +\underbrace{\frac{v \partial^{2}(\overline{q^{\prime} u^\prime_i})}{\partial x_{j}^{2}}}_{IX}
# -2 v \overline { \underbrace{\left ( \frac{\partial u_i^\prime }{ \partial x_j} \right ) \left ( \frac{\partial q^\prime} {\partial x_j} \right )}_{X} }
# \end{align}

# %% [markdown]
# ## What are the rms perturbation profiles?

# %%
def plot_pert(ax,keep_dict,key):
    """
    make a quick horizontal rms vs. height plot of variable key in keep_dict
    
    Pararmeters
    -----------
    
    ax: matplotlib axis
    
    keep_dict: dict
      dictionary with variables
        
    key: string
      name of variable to plot
      
    Returns
    -------
    
    ax: matplotlib axis
      same axis with plot
    """
    the_var = (keep_dict[key]**2.).mean(axis=(1,2))
    ax.plot(the_var**0.5,keep_dict['z'],'r-')
    ax.set_title(key)
    
fig, axes = plt.subplots(4,2,figsize=(10,10),constrained_layout=True)
flat_axes = list(axes.flat)
for ax, varname in zip(flat_axes,the_vars):
    plot_pert(ax,keep_dict,varname)
fig.delaxes(flat_axes[-1])



# %%
#help(axes[0,0].plot)

# %% [markdown]
# ## set scaling parameters q* and w* from rms values

# %%
#
# take zi = 1500 m and find scaling variables
#
zi = 1500.
hit = zvals < zi
w_pert = (keep_dict['w_pert'][hit,:,:]**2.).mean(axis=(1,2))**0.5
wstar = w_pert.mean()
print(f"wstar={wstar}")
q_pert = (keep_dict['qv_pert'][hit,:,:]**2.).mean(axis=(1,2))**0.5
qstar = q_pert.mean()
print(f"qstar={qstar*kg2g}")
term_scale = wstar**2.*qstar/zi
print(f"term_scale = {term_scale}")


# %% [markdown]
# ## Look at mean gradient terms in XI
#
# \begin{equation*}
# \underbrace{\overline{u^\prime_{i} u^\prime_{j}} \frac{\partial \bar{q}}{\partial x_{j}}}_{XI}
# \end{equation*}

# %%
def do_partial(keep_dict,varname,axis_val):
    """
    calculate partial deriv along axis_val using
    simple differencing
    
    Paramteters
    -----------
    
    keep_dict: dict
       dictionary with variables
       
    varname: string
       name of variable to differentiate
       
    axis_val: int
       axis to di
    """
    axis_dict={0:'x',1:'y',2:'z'}
    the_var = keep_dict[varname]
    numerator = np.diff(the_var,axis=axis_val)
    coord_var = axis_dict[axis_val]
    denominator = np.diff(keep_dict[coord_var])
    return numerator/denominator
    


# %%
fig, ax = plt.subplots(1,1)
ax.plot('qv_avg','z','b-',data=keep_dict)
ax.set_title('average vertical qv profile');


# %% [markdown]
# ## Plot the three equations for term XI

# %%
def calc_middle(vector):
    middle = (vector[1:] + vector[:-1])/2.
    return middle

def do_xi(keep_dict,hit):
    qv=keep_dict['qv_avg'][hit]
    z=keep_dict['z'][hit]
    qvavg_grad = np.diff(qv)/np.diff(qv)
    wperturb = keep_dict['w_pert'][hit,:,:]
    uperturb = keep_dict['u_pert'][hit,:,:]
    vperturb = keep_dict['v_pert'][hit,:,:]
    w_flux = ((wperturb.T*wperturb.T)[:,:,1:]*qvavg_grad).mean(axis=(0,1))
    u_flux = ((uperturb.T*wperturb.T)[:,:,1:]*qvavg_grad).mean(axis=(0,1))
    v_flux = ((vperturb.T*wperturb.T)[:,:,1:]*qvavg_grad).mean(axis=(0,1))
    return w_flux, v_flux, u_flux

hit = keep_dict['z']< 1500
w_flux, v_vlux, uflux = do_xi(keep_dict,hit)
z_coord = calc_middle(keep_dict['z'][hit])
fig,ax = plt.subplots(1,1,figsize=(10,10))
ax.plot(w_flux,z_coord);

