#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 11:42:52 2019

@author: rodell
"""





from cr500.utils import read_cabauw
from pathlib import Path
import xarray as xr
import context
import os
import time



data_dir= str(context.pro_data_dir)
files = sorted(list(Path(data_dir).glob("*.nc")))
#print(files)

#glob('cesar*.nc')



#print(files)

#all_files = []
#for data_file in sorted(os.listdir(files)):
#    print(data_file)

var_dict = read_cabauw.read(files)

#for key in doppler_dict:
#    print(key)
#doppler_dict.keys

ds_insitu  = xr.Dataset(var_dict)

#ds_doppler = xr.Dataset(doppler_dict)


#test = ds_insitu.where(ds_insitu.time<100)


#ds_doppler



#xr.combine_nested(xr,'time')