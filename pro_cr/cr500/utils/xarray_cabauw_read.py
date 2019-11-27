#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 13:27:58 2019

@author: rodell
"""

import xarray
from netCDF4 import Dataset
from pathlib import Path
import context
from datetime import datetime
from pytz import utc
import matplotlib
from matplotlib import pyplot as plt
import pprint
import xarray as xr
matplotlib.use("Agg")
from cr500.utils import ncdump
import re
#
# the following regular expression captures one group
# of exactly 3 characters, all numbers between 0-9
# the filenames look like ncep_gec00.t00z.pgrb2f006_SA.nc
#
find_hour = re.compile(r".*grb2.*(\d{3,3}).*_SA.nc")
#
#
#
data_dir = "/scratch/paustin/ewicksteed/mydata/2016060100"
all_files = list(Path(data_dir).glob("*.nc"))
pprint.pprint(all_files)


def sort_hour(the_file):
    """
    sort the files by converting the 3 digit time to
    an integer and returning that number
    """
    the_match = find_hour.match(str(the_file))
    return int(the_match.group(1))


if __name__ == "__main__":
    all_files.sort(key=sort_hour)
    xarray_files = []
    for item in all_files:
        with Dataset(str(item)) as nc_file:
            the_time = nc_file.variables['time'][...]
            print(datetime.fromtimestamp(the_time, tz=utc))
            ds = xr.open_dataset(item)
            xarray_files.append(ds)
    ds_big = xr.combine_nested(xarray_files, 'time')
    time_average = ds_big.mean('time')
    #
    # time_average.data_vars
    # time_average.coords
    varnames = list(ds_big.variables.keys())
    #
    #
    # create an xarray out of these files
    #
    vel_vals = [
        'VVEL_200mb', 'VVEL_250mb', 'VVEL_500mb', 'VVEL_700mb', 'VVEL_925mb',
        'VVEL_1000mb'
    ]
    vel_dict = {}
    for key in vel_vals:
        vel_dict[key] = ds_big.variables[key]
    ds_small = xr.Dataset(vel_dict, ds_big.coords)