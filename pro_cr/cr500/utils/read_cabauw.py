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
import urllib,os
import context
from pathlib import Path
import cr500
from cr500.utils import ncdump

import numpy as np
from bs4 import BeautifulSoup
import requests







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
 
    


#all_files=context.pro_data_dir.glob('cesar*.nc')


def read(files):
    """
    Reads cabuaw surface flux, surface met, tower and doppler data into a directory 
    
    """
    ## Initialize dictionary
    var_dict={}
    
    ## Sort Files...really important for the processed_multi-beam files
#    all_files = sorted(files, key = lambda file: os.path.getctime(file))
#    all_files = all_files[::-1]
    
#    all_files = []
#    for data_file in sorted(os.listdir(files)):
#        print(data_file)
#        

    ## Loop files
    doppler_list = []
    for the_file in files:
        if str(the_file).find('nubiscope') > -1:
            continue
        print(the_file)
        ## Read Files
        with Dataset(the_file,'r') as f:
            details=f.variables['iso_dataset']
            attr_dict=get_attrs(details)
            lon=attr_dict['westbound_longitude']
            lat=attr_dict['northbound_latitude']
            title=attr_dict['title'].split()
            filetype='{}_{}'.format(*title[1:3])
#            print(details)
            print(filetype)
            

            if filetype == 'meteorological_surface':
                print('Surf Met')
                for var in ['P0','RAIN']:
                    var_array = f.variables[var][...]
                    dims = f.variables[var].dimensions
                    var_dict.update({var : (dims,np.array(var_array))})
                    
                the_time=f.variables['time'][...]
                dims=f.variables['time'].dimensions
                start_date=f.variables['product'].date_start_of_data
                start_date = parse(start_date)
                time_vec=[]
                for the_hour in the_time:
                    time_vec.append(start_date + datetime.timedelta(hours=float(the_hour)))
                time_vec=np.array(time_vec)
                var_dict.update({'Datetime_Met' : (dims,time_vec)})
                time_vec=[]
    
            elif filetype == 'tower_meteorological':
    #            ncdump.ncdump(f)
                print('Tower')
                for var in ['TA','F','D','z']:
                    var_array = f.variables[var][...]
                    dims = f.variables[var].dimensions
                    var_dict.update({var : (dims,np.array(var_array))})
                    
                the_time=f.variables['time'][...]
                dims=f.variables['time'].dimensions
                start_date=f.variables['product'].date_start_of_data
                start_date = parse(start_date)
                time_vec=[]
                for the_hour in the_time:
                    time_vec.append(start_date + datetime.timedelta(hours=float(the_hour)))
                time_vec=np.array(time_vec)
                var_dict.update({'Datetime_Tower' : (dims,time_vec)})
                time_vec=[]
    
            elif filetype == 'surface_fluxes':
                print('Surf Flux')
                for var in ['UST','H']:
                    var_array = f.variables[var][...]
                    dims = f.variables[var].dimensions
                    var_dict.update({var : (dims,np.array(var_array))})
                    
                the_time=f.variables['time'][...]
                dims=f.variables['time'].dimensions
                start_date=f.variables['product'].date_start_of_data
                start_date = parse(start_date)
                time_vec=[]
                for the_hour in the_time:
                    time_vec.append(start_date + datetime.timedelta(hours=float(the_hour)))
                time_vec=np.array(time_vec)
                var_dict.update({'Datetime_Flux' : (dims,time_vec)})
                time_vec=[]
                                
            elif filetype == 'processed_multi-beam':
                dop_dict = {}
                for var in ['vertical_velocity','horizontal_wind_speed','horizontal_wind_direction']:
                    var_i=f.variables[var][...]
                    scale_factor=f.variables[var].scale_factor
                    var_array=(var_i*scale_factor) 
                    dims = f.variables[var].dimensions
                    dop_dict.update({var : (dims,np.array(var_array))})
                        
                for var in ['range_resolution','height_1st_interval','time','product']:
                    var_array = f.variables[var][...]
                    dims = f.variables[var].dimensions
                    dop_dict.update({var : (dims,np.array(var_array))})   
                    
                the_time=f.variables['time'][...]
                dims=f.variables['time'].dimensions
                start_date=f.variables['product'].date_start_of_data
                start_date = parse(start_date)
                time_vec=[]
                for the_hour in the_time:
                    time_vec.append(start_date + datetime.timedelta(hours=float(the_hour)))
                time_vec=np.array(time_vec)
                dop_dict.update({'Datetime_Doppler' : (dims,time_vec)})
                time_vec=[]

                doppler_list.append(dop_dict)
        
            else:
                raise ValueError("didn't recognize {}".format(filetype))
                

#    d = {}
#    for k in doppler_list[0].keys():
#        d[k] = np.concatenate(list(d[k] for d in doppler_list))
        
#    final_dict = {**var_dict, **d}

    return(var_dict)
#    return(final_dict)  

       

#
#def read(files):
#    """
#    Reads cabuaw surface flux, surface met, tower and doppler data into a directory 
#    
#    """
#    ## Initialize dictionary
#    var_dict={}
#    
#    ## Sort Files...really important for the processed_multi-beam files
#    all_files = sorted(files, key = lambda file: os.path.getctime(file))
#    all_files = all_files[::-1]
#    
#    ## Loop files
#    doppler_list = []
#    for the_file in all_files:
#        if str(the_file).find('nubiscope') > -1:
#            continue
#        print(the_file)
#        ## Read Files
#        with Dataset(the_file,'r') as f:
#            details=f.variables['iso_dataset']
#            attr_dict=get_attrs(details)
#            lon=attr_dict['westbound_longitude']
#            lat=attr_dict['northbound_latitude']
#            title=attr_dict['title'].split()
#            filetype='{}_{}'.format(*title[1:3])
##            print(details)
#            print(filetype)
#            
#
#            if filetype == 'meteorological_surface':
#                print('Surf Met')
#                for var in ['P0','RAIN']:
#                    var_array = f.variables[var][...]
#                    dims = f.variables[var].dimensions
#                    var_dict.update({var : (dims,var_array)})
#                    
#                the_time=f.variables['time'][...]
#                start_date=f.variables['product'].date_start_of_data
#                start_date = parse(start_date)
#                time_vec=[]
#                for the_hour in the_time:
#                    time_vec.append(start_date + datetime.timedelta(hours=float(the_hour)))
#                time_vec=np.array(time_vec)
#                var_dict['Datetime_Met'] = time_vec
#                time_vec=[]
#    
#            elif filetype == 'tower_meteorological':
#    #            ncdump.ncdump(f)
#                print('Tower')
#                for var in ['TA','F','D','z']:
#                    var_dict[var] = f.variables[var][...]
#
#                    
#                the_time=f.variables['time'][...]
#                start_date=f.variables['product'].date_start_of_data
#                start_date = parse(start_date)
#                time_vec=[]
#                for the_hour in the_time:
#                    time_vec.append(start_date + datetime.timedelta(hours=float(the_hour)))
#                time_vec=np.array(time_vec)
#                var_dict['Datetime_Tower'] = time_vec
#                time_vec=[]
#    
#            elif filetype == 'surface_fluxes':
#                print('Surf Flux')
#                for var in ['UST','H']:
#                    var_dict[var] = f.variables[var][...]
#                    
#                the_time=f.variables['time'][...]
#                start_date=f.variables['product'].date_start_of_data
#                start_date = parse(start_date)
#                time_vec=[]
#                for the_hour in the_time:
#                    time_vec.append(start_date + datetime.timedelta(hours=float(the_hour)))
#                time_vec=np.array(time_vec)
#                var_dict['Datetime_Flux'] = time_vec
#                time_vec=[]
#                                
#            elif filetype == 'processed_multi-beam':
#                dop_dict = {}
#                for var in ['vertical_velocity','horizontal_wind_speed','horizontal_wind_direction']:
#                    var_i=f.variables[var][...]
#                    scale_factor=f.variables[var].scale_factor
#                    dop_dict[var]=(var_i*scale_factor) 
#                        
#                for var in ['range_resolution','height_1st_interval','time','product']:
#                    dop_dict[var] = f.variables[var][...]
#       
#                the_time=f.variables['time'][...]
#                start_date=f.variables['product'].date_start_of_data
#                start_date = parse(start_date)
#                time_vec=[]
#                for the_hour in the_time:
#                    time_vec.append(start_date + datetime.timedelta(hours=float(the_hour)))
#                time_vec=np.array(time_vec)
#                dop_dict['Datetime_Doppler'] = time_vec
#                time_vec=[]
#
#
##                dop_dict['start_time'] = f.variables['product'].date_start_of_data
##                dop_dict['stop_time'] = f.variables['product'].date_end_of_data
#
#
#                doppler_list.append(dop_dict)
#        
#            else:
#                raise ValueError("didn't recognize {}".format(filetype))
#                
#
#    d = {}
#    for k in doppler_list[0].keys():
#        d[k] = np.concatenate(list(d[k] for d in doppler_list))
#    final_dict = {**var_dict, **d}
#
#    return(var_dict,doppler_list, dop_dict,final_dict )
##    return(final_dict)  

       


        
        
        
        
        
        