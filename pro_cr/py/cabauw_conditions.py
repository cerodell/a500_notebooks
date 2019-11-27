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



data_dir= str(context.pro_data_dir)
files = sorted(list(Path(data_dir).glob("*.nc")))
#print(files)


var_list, doppler_list, fill_var = read_cabauw.read(files)



ds_insitu = read_cabauw.xarray_insitu(var_list)
ds_doppler = read_cabauw.xarray_doppler(doppler_list)


#for key in ds_doppler:
#    zzz = ds_doppler.where(key == fill_var, key ,  np.nan)
    
test = doppler_list[0]

#ds_hour_insitu  = ds_insitu.groupby('Datetime_Met.hour').mean('time')
#ds_hour_doppler = ds_doppler.groupby('Datetime_Doppler.hour').mean('time')


#ds_doppler.horizontal_wind_speed



#ds_insitu.UST



















