#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 10:54:24 2019

@author: rodell
"""

import glob
from netCDF4 import Dataset
from dateutil.parser import parse
import datetime
import numpy as np
import matplotlib.dates as mdates
import urllib,os
import context
import pandas as pd
from pathlib import Path
from netCDF4 import Dataset
import cr500
from cr500.utils import download, plt_set, wind_eq, ncdump

import numpy as np
from bs4 import BeautifulSoup
import requests

from matplotlib import pyplot as plt
#import matplotlib.patches as mpatches

import math



filelist=['cesar_surface_meteo_lc1_t10_v1.0_201807.nc',
          'cesar_tower_meteo_lc1_t10_v1.0_201807.nc',
          'cesar_surface_flux_lc1_t10_v1.0_201807.nc']


def get_attrs(ncvar):
    """
        get every attribute of a netcdf variable
        
    Parameters
    ----------
    
    ncvar: netcdf variable object
    Returns
    -------
    
    attr_dict: dict
      dictionary with attribute names and values
    """
    attributes=ncvar.ncattrs()
    attr_dict={}
    for attr in attributes:
        attr_dict[attr]=getattr(ncvar,attr)
    return attr_dict
 
    


all_files=context.pro_data_dir.glob('cesar*.nc')

for the_file in all_files:
    if str(the_file).find('nubiscope') > -1:
        continue
    print(the_file)
    with Dataset(the_file,'r') as f:
        details=f.variables['iso_dataset']
        attr_dict=get_attrs(details)
        lon=attr_dict['westbound_longitude']
        lat=attr_dict['northbound_latitude']
        title=attr_dict['title'].split()
        filetype='{}_{}'.format(*title[1:3])

        print(filetype)
        if filetype == 'meteorological_surface':
            print('Surf Met')
            press_sfc=f.variables['P0'][...]      ## Atmospheric air pressure     units = (hPa)
            temp_2m=f.variables['TA002'][...]     ## Air temperature at 2 m       units = (degC)
            wsp_10m=f.variables['F010'][...]      ## Wind speed                   units = (m s^-1)
            wdir_10m=f.variables['D010'][...]     ## Wind Direction               units = (degs)
            precip=f.variables['RAIN'][...]       ## Rain amount                  units = (mm)
        elif filetype == 'tower_meteorological':
#            ncdump.ncdump(f)

            print('Tower')
            temp_tower=f.variables['TA'][...]     ## Air temperature              units = (K)
            wsp_tower=f.variables['F'][...]       ## Wind speed                   units = (m s^-1)
            wsp_tower.T
            wdir_tower=f.variables['D'][...]      ## Wind Direction               units = (degs)
            z=f.variables['z'][...]               ## Height above surfac          units = (m)
        elif filetype == 'surface_fluxes':
            print('Surf Flux')
            u_str=f.variables['UST'][...]         ## Friction velocity            units = (m s^-1)
            sen_heat_flux=f.variables['H'][...]   ## Surface sensible heat flux   units = (W m^-2)
        else:
            print('No')
            
#        ncdump.ncdump(nc_in)
        

       




        
        
        
        
        
        