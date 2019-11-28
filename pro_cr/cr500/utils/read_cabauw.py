#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 10:54:24 2019

@author: rodell
"""

#import glob
from netCDF4 import Dataset
from dateutil.parser import parse
import datetime
import numpy as np
#import urllib,os
#import context
#from pathlib import Path
#import cr500
from cr500.utils import ncdump
from pytz import utc

#import numpy as np
#from bs4 import BeautifulSoup
#import requests
import xarray as xr






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
 
    


def maketime(file):
    """
    Convert hour sense time to datetime np64
    Parameters
    ----------
    
    file: netcdf file
    Returns
    -------
    
    time_vec, dims: time array and dimentions of time
    
    """
    the_time=file.variables['time'][...]
    dims=file.variables['time'].dimensions
    start_date=file.variables['product'].date_start_of_data
    start_date = parse(start_date)
    time_vec=[]
    for the_hour in the_time:
        time_vec.append(start_date + datetime.timedelta(hours=float(the_hour)))
    time_vec=np.array(time_vec)
    return(time_vec, dims)
    


def read(files):
    """
    Reads cabuaw surface flux, surface met, tower and doppler data into a directory 
    
    Parameters
    ----------
    
    files: netcdf file
    Returns
    -------
    
    var_dict, doppler_list: dictoray of dims and array for the flux met and tower netcdfs
                            and a list dictorays of dims and array for the doppler netcdfs
    
    """
    ## Initialize dictionary
#    var_dict={}
    

    ## Loop files
    var_list, doppler_list = [], []
    for the_file in files:
        if str(the_file).find('nubiscope') > -1:
            continue
        print(the_file)
        ## Read Files
        with Dataset(the_file,'r') as f:
            details=f.variables['iso_dataset']
            attr_dict=get_attrs(details)
            title=attr_dict['title'].split()
            filetype='{}_{}'.format(*title[1:3])
#            print(details)
            print(filetype)
            

            if filetype == 'meteorological_surface':
                var_dict={}
                print('Surf Met')
                for var in ['P0','RAIN']:
                    var_array = f.variables[var][...]
                    dims = f.variables[var].dimensions
                    var_dict.update({var : (dims,np.array(var_array))})
                
                time_vec, dims = maketime(f)
                var_dict.update({'Datetime_Met' : (dims,time_vec)})
                time_vec=[]
                
                var_list.append(var_dict)

            elif filetype == 'tower_meteorological':
                var_dict={}
                print('Tower')
                for var in ['TA','F','D','z']:
                    fill_var = f.variables[var]._FillValue
                    var_array = f.variables[var][...]
                    dims = f.variables[var].dimensions
                    var_dict.update({var : (dims,np.array(var_array))})
                
                time_vec, dims = maketime(f)
                var_dict.update({'Datetime_Tower' : (dims,time_vec)})
                time_vec=[]
                
                var_list.append(var_dict)
                
            elif filetype == 'surface_fluxes':
                var_dict={}
                print('Surf Flux')
                for var in ['UST','H']:
                    var_array = f.variables[var][...]
                    dims = f.variables[var].dimensions
                    var_dict.update({var : (dims,np.array(var_array))})
                
                time_vec, dims = maketime(f)
                var_dict.update({'Datetime_Flux' : (dims,time_vec)})
                time_vec=[]
                
                var_list.append(var_dict)
                
            elif filetype == 'processed_multi-beam':
                dop_dict = {}
                for var in ['vertical_velocity','horizontal_wind_speed','horizontal_wind_direction']:
                    var_i=f.variables[var][...]
#                    fill_var = f.variables[var]._FillValue
#                    print(fill_var)
#                    var_i[var_i == fill_var] = fill_var
#                    var_ii = var_i
#                    print("After nan fill")
#                    print(var_i)
                    scale_factor=f.variables[var].scale_factor
#                    print(scale_factor)
                    var_array=(var_i*scale_factor) 
#                    var_array = np.where(var_ii == fill_var, var_ii, np.empty)
                    var_array[var_array == fill_var] = np.NAN

#                    print(np.max(var_array))
#                    print(np.min(var_array))

                    dims = f.variables[var].dimensions
                    dop_dict.update({var : (dims,np.array(var_array))})
                        
                for var in ['range_resolution','height_1st_interval']:
                    var_array = f.variables[var][...]
                    dims = f.variables[var].dimensions
                    dop_dict.update({var : (dims,np.array(var_array))})   
                                
                time_vec, dims = maketime(f)
    
                dop_dict.update({'Datetime_Doppler' : (dims,time_vec)})
                time_vec=[]

                doppler_list.append(dop_dict)
        
            else:
                raise ValueError("didn't recognize {}".format(filetype))
                

    return(var_list, fill_var)


def xarray_doppler(dict_list): 
    xarray_files = []
    for index in dict_list:
        ds  = xr.Dataset(index)
        xarray_files.append(ds)
    ds_final = xr.combine_nested(xarray_files, 'time')
    return(ds_final)

    
def xarray_insitu(dict_list):
    xarray_files = []
    for index in dict_list:
        ds  = xr.Dataset(index)
        xarray_files.append(ds)
    ds_final = xr.merge(xarray_files)    
    return(ds_final)



