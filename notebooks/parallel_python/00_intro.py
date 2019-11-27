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

# %% [markdown] {"slideshow": {"slide_type": "slide"}, "toc": true}
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#00---Introduction" data-toc-modified-id="00---Introduction-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>00 - Introduction</a></span><ul class="toc-item"><li><ul class="toc-item"><li><span><a href="#Concurrency-vs.-parallelism" data-toc-modified-id="Concurrency-vs.-parallelism-1.0.1"><span class="toc-item-num">1.0.1&nbsp;&nbsp;</span>Concurrency vs. parallelism</a></span></li><li><span><a href="#Threads-and-processes" data-toc-modified-id="Threads-and-processes-1.0.2"><span class="toc-item-num">1.0.2&nbsp;&nbsp;</span>Threads and processes</a></span></li><li><span><a href="#Thread-scheduling" data-toc-modified-id="Thread-scheduling-1.0.3"><span class="toc-item-num">1.0.3&nbsp;&nbsp;</span>Thread scheduling</a></span></li><li><span><a href="#Releasing-the-GIL" data-toc-modified-id="Releasing-the-GIL-1.0.4"><span class="toc-item-num">1.0.4&nbsp;&nbsp;</span>Releasing the GIL</a></span></li></ul></li></ul></li></ul></div>

# %%
from IPython.display import Image

# %% [markdown]
# # 00 - Introduction

# %% [markdown] {"slideshow": {"slide_type": "slide"}}
#
# %% [markdown]
# ### Concurrency vs. parallelism
#
# [Google says: ](https://www.google.ca/search?q=concurrency+vs.+parallelism&rlz=1C5CHFA_enCA698CA698&oq=conncurrency+vs.+parallel&aqs=chrome.1.69i57j0l5.6167j0j7&sourceid=chrome&ie=UTF-8)

# %% [markdown]
# ### Threads and processes
#
# > From [Wikipedia](https://en.wikipedia.org/wiki/Thread_(computing)):
#
# >> "In computer science, a thread of execution is the smallest sequence of programmed instructions that can be managed independently by a scheduler, which is typically a part of the operating system.[1] The implementation of threads and processes differs between operating systems, but in most cases a thread is a component of a process. Multiple threads can exist within one process, executing concurrently and sharing resources such as memory, while different processes do not share these resources. In particular, the threads of a process share its executable code and the values of its variables at any given time."
#

# %% [markdown]
# #### Threads and processes in Python
#
# [Reference: Thomas Moreau and Olivier Griesel, PyParis 2017 [Mor2017]](https://tommoral.github.io/pyparis17/#1)
#
# #### Python global intepreter lock
#
# 1. Motivation: python objects (lists, dicts, sets, etc.) manage their own memory by storing a counter that keeps track of how many copies of an object are in use.  Memory is reclaimed when that counter goes to zero.
#
# 1. Having a globally available reference count makes it simple for Python extensions to create, modify and share python objects.
#
# 1. To avoid memory corruption, a python process will only allow 1 thread at any given moment to run python code.  Any thread that wants to access python objects in that process needs to acquire the global interpreter lock (GIL).
#
# 1. A python extension written in C, C++ or numba is free to release the GIL, provided it doesn't create, destroy or modify any python objects.  For example: numpy, pandas, scipy.ndimage, scipy.integrate.quadrature all release the GIL
#
# 1. Many python standard library input/output routines (file reading, networking) also release the GIL
#
# 1. On the other hand:  hdf5, and therefore h5py and netCDF4, don't release the GIL and are single threaded.
#
# 1. Python comes with many libraries to manage both processes and threads.
#

# %% [markdown]
# ### Thread scheduling
#
# If multiple threads are present in a python process, the python intepreter releases the GIL at specified intervals (5 miliseconds default) to allow them to execute:

# %%
Image(filename='images/morreau1.png')  #[Mor2017]

# %% [markdown]
# #### Note that these three threads are taking turns, resulting in a computation that runs slightly slower (because of overhead) than running on a single thread

# %% [markdown]
# ### Releasing the GIL
#
# If the computation running on the thread has released the GIL, then it can run independently of other threads in the process.  Execution of these threads are scheduled by the operating system along with all the other threads and processes on the system.
#
# In particular, basic computation functions in Numpy, like (\__add\__ (+), \__subtract\__ (-) etc. release the GIL, as well as universal math functions like cos, sin etc.

# %%
Image(filename='images/morreau2.png')  #[Morr2017]

# %%
Image(filename='images/morreau3.png') #[Morr2017]
