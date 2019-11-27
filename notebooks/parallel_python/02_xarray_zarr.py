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

# # 02: xarray, netcdf and zarr
#
# Motivation:  how you store your data can an enormous effect on performance.
#
# Four issues:
#
# 1) Compression vs. cpu time to uncompress/compress
#
# 2) Multithreaded read/writes
#
# 3) Performance for cloud storage (Amazon, Google Compute, Azure)
#
# 4) Throttling data reads to fit in available memory and avoid swapping ("chunking")
#
# ## The current defacto standard in atmos/ocean science
#
#
# * [netcdf/hdf5](ttps://www.unidata.ucar.edu/software/netcdf/docs/index.html)
#

# %% [markdown]
# ## Some challenges with netcdf
#
# * [The Pangeo project](http://pangeo-data.org/)
#
# * [Cloud challenges](https://medium.com/pangeo)

# %% [markdown]
# ## create an xarray

# %%
import glob
import numpy as np
import pdb
import xarray
from matplotlib import pyplot as plt

# %%
# #!conda install -y xarray

# %%
import context
from a500.utils.data_read import download
from pathlib import Path

# %% [markdown]
# ## Download toy model data

# %%
# get 10 files, each is the same timestep for a member of a
# 10 member ensemble



root = "http://clouds.eos.ubc.ca/~phil/docs/atsc500/small_les"
for i in range(10):
    the_name = f"mar12014_{(i+1):d}_15600.nc"
    print(the_name)
    url = "{}/{}".format(root, the_name)
    download(the_name, root=root, dest_folder=context.data_dir)

# %%
# get 10 files, each is the same timestep for a member of a
# 10 member ensemble


the_files = context.data_dir.glob("mar12*nc")
the_files = list(the_files)
the_files

# %% [markdown]
# ## Sort in numeric order

# %%



def sort_name(pathname):
    """
      sort the filenames so '10' sorts
      last by converting to integers
    """
    str_name = str(pathname.name)
    front, number, back = str_name.split("_")
    return int(number)


# %% [markdown]
# ## Make an xarray

# %% [markdown]
# Now use xarray to stitch together the 10 ensemble members along a new "virtual dimenstion".
# The variable "ds"  is an xray dataset, which controls the reads/writes from the
# 10 netcdf files
#

# %%
the_files.sort(key=sort_name)

#
#  put the 10 ensembles together along a new "ens" dimension
#  using an xray multi-file dataset
#
ds = xarray.open_mfdataset(the_files, engine="netcdf4", concat_dim="ens")
# dump the structure
print(ds)
#
#  3-d ensemble average for temp
#
x = ds["x"]
y = ds["y"]
z = ds["z"]
temp = ds["TABS"]
mean_temp = temp[:, :, :, :].mean(dim="ens")
#
# same for velocity
#
wvel = ds["W"]
mean_w = wvel[:, :, :, :].mean(dim="ens")
#
# now look at the perturbation fields for one ensemble member
#
wvelprime = wvel[0, :, :, :] - mean_w
Tprime = temp[0, :, :, :] - mean_temp
flux_prime = wvelprime * Tprime
flux_profile = flux_prime.mean(dim="x").mean(dim="y")
keep_dict = dict(
    flux_prof=flux_profile,
    flux_prime=flux_prime.values,
    wvelprime=wvelprime.values,
    Tprime=Tprime.values,
    x=x,
    y=y,
    z=z,
)

# %% [markdown]
# ## Dump to a zarr file

# %%
ds.to_zarr("zarr_dir", "w")

# %%
fig1, ax1 = plt.subplots(1, 1)
ax1.plot("flux_prof", "z", data=keep_dict)
ax1.set(
    title="Ens 0: vertically averaged kinematic heat flux",
    ylabel="z (m)",
    xlabel="flux (K m/s)",
)

fig2, ax2 = plt.subplots(1, 1)
z200 = np.searchsorted(keep_dict["z"], 200)
flux_200 = keep_dict["flux_prime"][z200, :, :].flat
ax2.hist(flux_200)
ax2.set(title="histogram of kinematic heat flux (K m/s) at z=200 m")

fig3, ax3 = plt.subplots(1, 1)
wvel200 = keep_dict["wvelprime"][z200, :, :].flat
ax3.hist(wvel200)
ax3.set(title="histogram of wvel' at 200 m")

fig4, ax4 = plt.subplots(1, 1)
Tprime200 = keep_dict["Tprime"][z200, :, :].flat
ax4.hist(Tprime200)
ax4.set(title="histogram ot T' at z=200 m")
