#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 14:03:45 2019

@author: rodell
"""

from matplotlib import pyplot as plt
import numpy as np
import xarray as xr
import context




filein = str(context.pro_data_dir)+str('/xr/var_ds.zarr')
var_ds = xr.open_zarr(filein)


#wsp = var_ds.F.plot()

#var_ds.sel(time = slice('2001-01-02','2001-01-03'))

tmep_10m = var_ds.TA[:,0].mean()


hour = var_ds.groupby('time.hour').mean(dim='time')

tmep_10m_hour = hour.TA[:,0]


#inversion = hour.where(hour.TA[:,0] > hour.TA[:,6], drop=True)

inversion = var_ds.where(var_ds.TA[:,0] > var_ds.TA[:,6], drop=True)

#isel(time=0).plot(size=6)

inversion.TA[:,6].plot(size=6, x='time')


temp_2m = inversion.TA[:,6]



#hour.TA.where(hour.TA[:,0] > hour.TA[:,6], drop=True).plot.scatter()

#time = var_d