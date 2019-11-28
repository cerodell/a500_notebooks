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



data_dir_insitu = str(context.pro_data_dir)+str('/Insitu/')
files = sorted(list(Path(data_dir_insitu).glob("*.nc")))
#print(files)


var_list, fill_var = read_cabauw.read(files)



ds_insitu = read_cabauw.xarray_insitu(var_list)


#for key in ds_doppler:
#    zzz = ds_doppler.where(key == fill_var, key ,  np.nan)
    

ds_hour_insitu  = ds_insitu.groupby('Datetime_Met.hour').mean('time')
ds_height_insitu  = ds_insitu.groupby('z').mean('time')

#ds_hour_doppler = ds_doppler.groupby('Datetime_Doppler.hour').mean('time')


#ds_doppler.horizontal_wind_speed



wsp = np.array(ds_height_insitu.F)
wsp[wsp==-9999]=np.nan
height = np.array(ds_hour_insitu.z)


#wsp =


plt.plot(wsp.T,height)

#wsp



















