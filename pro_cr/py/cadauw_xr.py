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


#var_ds.where(var_ds. < 50, drop=True))



#time = var_d