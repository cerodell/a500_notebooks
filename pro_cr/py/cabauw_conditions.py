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




#data_dir_insitu = str(context.pro_data_dir)+str('/Doppler/')
#files = sorted(list(Path(data_dir_insitu).glob("*.nc")))
##print(files)
#
#
#doppler_list  = read_cabauw.read_doppler(files)
#
#doppler_list_0 = doppler_list[0]
#doppler_ds = read_cabauw.xarray_unlike(doppler_list)
#
#
#
#out_dir = str(context.pro_data_dir)
#out_dir = Path(str(context.pro_data_dir)+str('/xr/'))
#out_dir.mkdir(parents=True, exist_ok=True)
#
#
#full_dir = str(out_dir) + str(f"/doppler_ds.zarr")
#doppler_ds.compute()
#doppler_ds.to_zarr(full_dir, "w")
#print(f"wrote {out_dir}")


#data_dir_insitu = str(context.pro_data_dir)+str('/Insitu/')
#files = sorted(list(Path(data_dir_insitu).glob("*.nc")))
##print(files)
#
#
#var_list, met_list, tower_list, surf_list, fill_var = read_cabauw.read(files)
#
#
#met_ds = read_cabauw.xarray_unlike(met_list)
#tower_ds = read_cabauw.xarray_unlike(tower_list)
#surf_ds = read_cabauw.xarray_unlike(surf_list)
#
#xarray_list = [met_ds, tower_ds, surf_ds]
#var_ds = xr.merge(xarray_list)
#
#
#
#out_dir = str(context.pro_data_dir)
#out_dir = Path(str(context.pro_data_dir)+str('/xr/'))
#out_dir.mkdir(parents=True, exist_ok=True)
#
#
#
#
#full_dir = str(out_dir) + str(f"/surf_ds.zarr")
#surf_ds.compute()
#surf_ds.to_zarr(full_dir, "w")
#print(f"wrote {out_dir}")


