#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 11:42:52 2019

@author: rodell
"""

from cr500.utils import read_cabauw
from matplotlib import pyplot as plt
import numpy as np
from pathlib import Path
import xarray as xr
import context
import os
import time


###############################################################################################
"""################# Merge and write Surface Met, FLux and Tower DataArray  ######################"""
###############################################################################################

### Path to Data
data_dir_insitu = str(context.pro_data_dir)+str('/Insitu/')
files = sorted(list(Path(data_dir_insitu).glob("*.nc")))
#print(files)


### Call read cabauw function to compile lists of dictionaries in DataArray Formate 
var_list, met_list, tower_list, surf_list, fill_var = read_cabauw.read(files)


## Combine Nested  DataArrays from Lists of Dictionaries
met_ds = read_cabauw.xarray_unlike(met_list)
tower_ds = read_cabauw.xarray_unlike(tower_list)
surf_ds = read_cabauw.xarray_unlike(surf_list)

## Merge DataArrays
xarray_list = [met_ds, tower_ds, surf_ds]
var_ds = xr.merge(xarray_list)


## Create a directory to save DataArray (.zarr) file
out_dir = str(context.pro_data_dir)
out_dir = Path(str(context.pro_data_dir)+str('/xr/'))
out_dir.mkdir(parents=True, exist_ok=True)


## Write and save DataArray (.zarr) file
full_dir = str(out_dir) + str(f"/var_ds.zarr")
var_ds.compute()
var_ds.to_zarr(full_dir, "w")
print(f"wrote {out_dir}")




###############################################################################################
"""######################## Merge and write Doppler DataArray ###############################"""
###############################################################################################

#### Path to Data
#data_dir_insitu = str(context.pro_data_dir)+str('/Doppler/')
#files = sorted(list(Path(data_dir_insitu).glob("*.nc")))
##print(files)
#
#
#### Call read cabauw function to compile lists of dictionaries in DataArray Formate 
#doppler_list  = read_cabauw.read_doppler(files)
#doppler_list_0 = doppler_list[0]
#
### Combine Nested  DataArrays from Lists of Dictionaries
#doppler_ds = read_cabauw.xarray_unlike(doppler_list)
#
#
### Create a directory to save DataArray (.zarr) file
#out_dir = str(context.pro_data_dir)
#out_dir = Path(str(context.pro_data_dir)+str('/xr/'))
#out_dir.mkdir(parents=True, exist_ok=True)
#
#
### Write and save DataArray (.zarr) file
#full_dir = str(out_dir) + str(f"/doppler_ds.zarr")
#doppler_ds.compute()
#doppler_ds.to_zarr(full_dir, "w")
#print(f"wrote {out_dir}")









