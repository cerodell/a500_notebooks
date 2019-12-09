#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 10:54:24 2019

This script to read  NetCDF files for a variety of instruments from
the Cabauw study site and compile them into a list dictionary formatted
to be converted to Xarray.

@author: rodell
"""

#import glob
from netCDF4 import Dataset
from dateutil.parser import parse
import datetime
import numpy as np
from cr500.utils import ncdump
from pytz import utc
import xarray as xr



def nanfill(var, fill_var):
    var[var == fill_var] = np.nan
    return(var)


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
    var_list, met_list, tower_list, surf_list = [], [], [], []
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
                for var in ['P0','RAIN', 'TA002', 'Q002']:
                    fill_var = f.variables[var]._FillValue
                    var_array = f.variables[var][...]
                    array = np.array(var_array, dtype=float)
                    array[array == fill_var] = np.nan

                    dims = f.variables[var].dimensions
                    var_dict.update({var : (dims,np.array(array))})
                
                time_vec, dims = maketime(f)
                var_dict.update({'time' : (dims,time_vec)})
                time_vec=[]
                
                met_list.append(var_dict)
                var_list.append(var_dict)


            elif filetype == 'tower_meteorological':
                var_dict={}
                print('Tower')
                for var in ['TA','F','D','z']:
                    fill_var = f.variables[var]._FillValue
                    var_array = f.variables[var][...]
                    array = np.array(var_array, dtype=float)
                    array[array == fill_var] = np.nan

                    dims = f.variables[var].dimensions
                    var_dict.update({var : (dims,np.array(array))})
                
                time_vec, dims = maketime(f)
                var_dict.update({'time' : (dims,time_vec)})
                time_vec=[]
                
                tower_list.append(var_dict)
                var_list.append(var_dict)

                
            elif filetype == 'surface_fluxes':
                var_dict={}
                print('Surf Flux')
                for var in ['UST','H', 'LE']:
                    fill_var = f.variables[var]._FillValue
                    var_array = f.variables[var][...]
                    array = np.array(var_array, dtype=float)
                    array[array == fill_var] = np.nan

                    dims = f.variables[var].dimensions
                    var_dict.update({var : (dims,np.array(array))})
                
                time_vec, dims = maketime(f)
                var_dict.update({'time' : (dims,time_vec)})
                time_vec=[]
                
                surf_list.append(var_dict)
                var_list.append(var_dict)

            else:
                raise ValueError("didn't recognize {}".format(filetype))
                

    return(var_list, met_list, tower_list, surf_list, fill_var)


def xarray_unlike(dict_list): 
    xarray_files = []
    for index in dict_list:
        ds  = xr.Dataset(index)
        xarray_files.append(ds)
    ds_final = xr.combine_nested(xarray_files, 'time')
    return(ds_final)

    
def xarray_like(dict_list):
    xarray_files = []
    for index in dict_list:
        ds  = xr.Dataset(index)
        xarray_files.append(ds)
    ds_final = xr.merge(xarray_files)    
    return(ds_final)






def read_doppler(files): 
    doppler_list = []
    for the_file in files:
        if str(the_file).find('nubiscope') > -1:
            continue
        print(the_file)
        
        with Dataset(the_file,'r') as f:
                details=f.variables['iso_dataset']
                attr_dict=get_attrs(details)
                title=attr_dict['title'].split()
                filetype='{}_{}'.format(*title[1:3])
        #            print(details)
                print(filetype)
                
                if filetype == 'processed_multi-beam':
                    dop_dict = {}
                    for var in ['vertical_velocity','horizontal_wind_speed','horizontal_wind_direction']:
                        var_i=f.variables[var][...]
                        fill_var = f.variables[var]._FillValue
                        scale_factor=f.variables[var].scale_factor
                        var_array=(var_i*scale_factor) 
                        array = np.array(var_array, dtype=float)
                        array[array == fill_var] = np.nan

    
#                        dims = f.variables[var].dimensions
                        dims = ('time','z_doppler')
                        dop_dict.update({var : (dims,np.array(array))})
                            
                    for var in ['height_1st_interval']:
                        var_array = f.variables[var][...]
                        array = np.array(var_array, dtype=float)

                        dims = ('z_doppler')
                        dop_dict.update({var : (dims,np.array(array))})   
                        
                    time_vec, dims = maketime(f)        
                    dop_dict.update({'time' : (dims,time_vec)})
                    time_vec=[]
    
                    doppler_list.append(dop_dict)

    
                else:
                    raise ValueError("didn't recognize {}".format(filetype))  
    return(doppler_list)
                    
                    
