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
# <div class="toc"><ul class="toc-item"><li><span><a href="#Using-the-fft-to-compute-correlation" data-toc-modified-id="Using-the-fft-to-compute-correlation-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Using the fft to compute correlation</a></span><ul class="toc-item"><li><span><a href="#Start-with-a-standard-autocorrelation" data-toc-modified-id="Start-with-a-standard-autocorrelation-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>Start with a standard autocorrelation</a></span></li><li><span><a href="#Now-do-this-using-Numerical-recipes-13.2.3" data-toc-modified-id="Now-do-this-using-Numerical-recipes-13.2.3-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>Now do this using Numerical recipes 13.2.3</a></span></li></ul></li></ul></div>

# %% [markdown]
# # Using the fft to compute correlation
#
# Below I use aircraft measurments of $\theta$ and wvel taken at 25 Hz.  I compute the 
# autocorrelation using numpy.correlate and numpy.fft and show they are identical, as we'd expect.  For discussion, See  [Numerical Recipes Chapter 13.2, p. 538](http://clouds.eos.ubc.ca/~phil/docs/atsc500/pdf_files/numerical_recipes_fft.pdf)
# (user: green, password: house)

# %%
from matplotlib import pyplot as plt
import context
from a500.utils.data_read import download
plt.style.use('ggplot')
import urllib
import os
download('aircraft.npz',root='http://clouds.eos.ubc.ca/~phil/docs/atsc500/data',
        dest_folder = context.data_dir)

# %% [markdown]
# ## Start with a standard autocorrelation

# %%
#http://stackoverflow.com/questions/643699/how-can-i-use-numpy-correlate-to-do-autocorrelation
import numpy as np
data = np.load(context.data_dir / 'aircraft.npz')
wvel=data['wvel'] - data['wvel'].mean()
theta=data['theta'] - data['theta'].mean()
autocorr = np.correlate(wvel,wvel,mode='full')
auto_data = autocorr[wvel.size:]
ticks=np.arange(0,wvel.size)
ticks=ticks/25.
fig,ax = plt.subplots(1,1,figsize=(10,8))
ax.set(xlabel='lag (seconds)',title='autocorrelation of wvel using numpy.correlate')
out=ax.plot(ticks[:300],auto_data[:300])


# %%
np.mean(wvel**2.)

# %% [markdown]
# ## Now do this using Numerical recipes 13.2.3
#
# $$
# \operatorname{Corr}(g, h)_{j} \Longleftrightarrow G_{k} H_{k}^{*}
# $$

# %%
import numpy.fft as fft
the_fft = fft.fft(wvel)
auto_fft = the_fft*np.conj(the_fft)
auto_fft = np.real(fft.ifft(auto_fft))

fig,ax = plt.subplots(1,1,figsize=(10,8))
ax.plot(ticks[:300],auto_fft[:300])
out=ax.set(xlabel='lag (seconds)',title='autocorrelation using fft')

# %%
data['readme']

# %%
len(data['wvel'])/25.

# %%
