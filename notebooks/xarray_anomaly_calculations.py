# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all
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

# %% [markdown] toc=true
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#Calculating-Climatologies-and-Anomalies-with-Xarray-and-Dask:" data-toc-modified-id="Calculating-Climatologies-and-Anomalies-with-Xarray-and-Dask:-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Calculating Climatologies and Anomalies with Xarray and Dask:</a></span><ul class="toc-item"><li><span><a href="#A-Workaround-for-a-Longstanding-Problem" data-toc-modified-id="A-Workaround-for-a-Longstanding-Problem-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>A Workaround for a Longstanding Problem</a></span><ul class="toc-item"><li><span><a href="#The-Dataset:-MERRA2-Daily-Surface-Temprature" data-toc-modified-id="The-Dataset:-MERRA2-Daily-Surface-Temprature-1.1.1"><span class="toc-item-num">1.1.1&nbsp;&nbsp;</span>The Dataset: MERRA2 Daily Surface Temprature</a></span></li><li><span><a href="#Default,-No-Rechunking" data-toc-modified-id="Default,-No-Rechunking-1.1.2"><span class="toc-item-num">1.1.2&nbsp;&nbsp;</span>Default, No Rechunking</a></span></li><li><span><a href="#With-Rechunking" data-toc-modified-id="With-Rechunking-1.1.3"><span class="toc-item-num">1.1.3&nbsp;&nbsp;</span>With Rechunking</a></span></li><li><span><a href="#The-Workaround:-Xarray.map_blocks" data-toc-modified-id="The-Workaround:-Xarray.map_blocks-1.1.4"><span class="toc-item-num">1.1.4&nbsp;&nbsp;</span>The Workaround: <code>Xarray.map_blocks</code></a></span></li></ul></li><li><span><a href="#Compute-the-climatology-and-anomalies-as-2D-maps" data-toc-modified-id="Compute-the-climatology-and-anomalies-as-2D-maps-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>Compute the climatology and anomalies as 2D maps</a></span></li></ul></li></ul></div>

# %% [markdown]
# # Calculating Climatologies and Anomalies with Xarray and Dask:
#
# ## A Workaround for a Longstanding Problem
#
# Climatologies are anomalies are a core operation in climate science. Many workflows start with the following operations:
# - Group spatiotemporal data by month or dayofyear (determined by the resolution of the dataset)
# - Take a mean of each group to determine the "climatology"
# - Broadcast the climatology back to the original dataset and subtract it, producing the "anomaly"
#
# Xarray makes this easy. We often write code like
#
#     gb = ds.groupby('time.month')
#     clim = gb.mean(dim='time')
#     anom = gb - clim
#
# Unfortunately there are problems related to how dask deals with this operation 
#
# - https://github.com/pydata/xarray/issues/1832
# - https://github.com/dask/dask/issues/874
# - https://github.com/pangeo-data/pangeo/issues/271
#
# There have been many attempted fixes over the years (see linked PRs above). But none of them has been totally successful.
#
# Here we desribe a new approach.

# %%
import xarray as xr
from dask.distributed import Client
import gcsfs
# %matplotlib inline
xr.__version__

# %% [markdown]
# ### The Dataset: MERRA2 Daily Surface Temprature

# %%
gcs = gcsfs.GCSFileSystem(token = 'anon')
to_map = gcs.get_mapper("ivanovich_merra2/t2maxdaily.zarr/")
ds = xr.open_zarr(to_map)
ds

# %%
ds.t2mmax.data

# %% [markdown]
# ### Default, No Rechunking

# %%
gb = ds.groupby('T.dayofyear')
clim = gb.mean(dim='T')
anom = gb - clim
anom_std = anom.std(dim='T')
anom_std.t2mmax.data

# %%
from dask.distributed import Client

client = Client("tcp://10.32.5.32:43525")
client

# %% [markdown]
# We see we have balooned up to almost 100,000 tasks

# %%
# %time anom_std.load()

# %% [markdown]
# Two minutes is a really long time to process 12 GB of data. And the dask cluster almost choked in the process.
#
# The parallelism became too fine-grained, resulting in too much communication overhead.
#
# ### With Rechunking
#
# Since the operation is embarassingly parallel in the space dimension, but the data are chunked in the time dimension, one idea is that rechunking could help.

# %%
ds_rechunk = ds.chunk({'T': -1, 'Y': 3})
ds_rechunk.t2mmax.data

# %%
gb = ds_rechunk.groupby('T.dayofyear')
clim = gb.mean(dim='T')
anom = gb - clim
anom_std = anom.std(dim='T')
anom_std.t2mmax.data


# %% [markdown]
# This created **4.5 million tasks**! Clearly not the solution we were hoping for. I'm not even going to try to compute it. For whatever reason, the way these operation (mostly indexing and broadcasting) are interpreted by dask array does not allow them to leverage the parallelism we know is possible.
#
# ### The Workaround: `Xarray.map_blocks`
#
# Since the computation is embarassingly parallel in the space dimension, I could use `dask.array.map_blocks` to operate on each chunk in isolation. The problem is, I don't know how to write the groupby and broadcasting logic in pure numpy. I need xarray and its indexes.
#
# The solution is xarray's new `map_blocks` function.

# %%
def calculate_anomaly(ds):
    # needed to workaround xarray's check with zero dimensions
    # https://github.com/pydata/xarray/issues/3575
    if len(ds['T']) == 0:
        return ds
    gb = ds.groupby("T.dayofyear")
    clim = gb.mean(dim='T')
    return gb - clim


# %%
t2mmax_anom = xr.map_blocks(calculate_anomaly, ds_rechunk.t2mmax)
t2mmax_anom.data

# %% [markdown]
# That seems great! Only 300 chunks! Let's see how it performs.

# %%
# %time t2mmax_std = t2mmax_anom.std(dim='T').load()

# %% [markdown]
# This was about twice as fast. Moreover, it feels like a more scalable approach.

# %%
t2mmax_std.plot()

# %% [markdown]
# ## Compute the climatology and anomalies as 2D maps
#
# The advantage of the `map_blocks` approach is that it doesn't create too many chuncks. That way we can lazily build more operations on top of the anomaly dataset.
#
# Below we count the number of "hot events" (anomaly > 1 degree for two consecutive days) per year.

# %%
rolling = t2mmax_anom.rolling(T = 2, center = True)
rolling_hot = rolling.max()
rolling_hot

# %%
yearly_events = (rolling_hot > 1).astype('int').resample(T='YS').sum()
yearly_events.data

# %%
yearly_events.load()

# %%
# skip 2019
yearly_events_mean = yearly_events[:-1].mean(dim='T')
yearly_events_anom = yearly_events[:-1] - yearly_events_mean

# %%
yearly_events_mean.plot()

# %%
yearly_events_anom[0].plot()

# %%
yearly_events_anom[-1].plot()

# %%
