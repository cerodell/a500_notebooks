#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 14:58:20 2019

@author: crodell
"""
'''

'''

import urllib,os
import context
import pandas as pd
from pathlib import Path
from netCDF4 import Dataset
from cr500.utils import download, plt_set, wind_eq

import numpy as np
from bs4 import BeautifulSoup
import requests

from matplotlib import pyplot as plt
#import matplotlib.patches as mpatches

import math 

#from  bl_pro.utils.data_read import download
#the_file='lidar/WC16_XPIA/WLS866-16_2015_06_07__00_00_00.rtd'


the_file='lidar/WLS7-61_2015_06_07__00_00_00.sta'

#out=download(the_file, dest_folder= context.pro_data_dir)
    
#####################################################################
'''######################  Read Files ############################'''
#################################################################### 

case_name=  context.pro_data_dir / the_file


df = pd.read_table(case_name, encoding ='ISO-8859-1',  skiprows=range(0, 56), sep = '\s+')





#####################################################################
'''######################  Make Plots ############################'''
#################################################################### 

heights = ['Vhm1','Vhm2','Vhm3','Vhm4','Vhm5','Vhm6','Vhm7','Vhm8','Vhm9','Vhm10']

agl = [40, 60, 80, 100, 120, 140, 160, 180, 200, 220]



fig,ax = plt.subplots(1,1,figsize=(16,8))
fig.suptitle('Wind Profile Lidar', fontsize=plt_set.fig_size, fontweight="bold")

for height in heights:  
    ax.plot(df[f'{height}'], agl, alpha = 0.5, label = f'{height}')   
ax.set_ylabel("Time", fontsize = plt_set.label)
ax.set_ylabel('Wind Speed ($ms^-1$)', fontsize = plt_set.label)
#ax.tick_params(axis='both', which='major', labelsize=plt_set.tick_size)
#ax.xaxis.grid(color='gray', linestyle='dashed')
#ax.yaxis.grid(color='gray', linestyle='dashed')
ax.legend()
ax.set_facecolor('lightgrey')
#plt.close('all')



#fig,ax = plt.subplots(1,1,figsize=(16,8))
#fig.suptitle('Wind Profile Lidar', fontsize=plt_set.fig_size, fontweight="bold")

#for height in heights:  
#    ax.plot(df['Date'],  df[f'{height}'],alpha = 0.5, label = f'{height}')   
#ax.set_ylabel("Time", fontsize = plt_set.label)
#ax.set_ylabel('Wind Speed ($ms^-1$)', fontsize = plt_set.label)
##ax.tick_params(axis='both', which='major', labelsize=plt_set.tick_size)
##ax.xaxis.grid(color='gray', linestyle='dashed')
##ax.yaxis.grid(color='gray', linestyle='dashed')
#ax.legend()
#ax.set_facecolor('lightgrey')
#plt.close('all')
  


z_i = 1000
m_bl = 8 


u_avg,u_perturb = wind_eq.do_reynolds(df['um1'])
v_avg,v_perturb = wind_eq.do_reynolds(df['vm1'])
w_avg,w_perturb = wind_eq.do_reynolds(df['wm1'])    
theta_avg,theta_perturb = wind_eq.do_reynolds(df['Tm']+273.15)    



u_str = wind_eq.u_str((u_perturb*w_perturb),(v_perturb*w_perturb))
w_str = wind_eq.w_str(theta_avg,z_i,(w_perturb*theta_perturb))
  

z = np.arange(0,1500,50)  


m_z = []
for i in range(len(z)):
    m_z_i = wind_eq.RxL(w_str,u_str,z[i],z_i,m_bl)
    m_z.append(m_z_i)
    
    
m_z = np.stack(m_z)


fig,ax = plt.subplots(1,1,figsize=(16,8))
fig.suptitle('Wind Profile Lidar', fontsize=plt_set.fig_size, fontweight="bold")


for i in range(len(u_str)):  
    ax.plot(m_z[:,i], z,alpha = 0.5)   
ax.set_ylabel("Time", fontsize = plt_set.label)
ax.set_ylabel('Wind Speed ($ms^-1$)', fontsize = plt_set.label)
#ax.tick_params(axis='both', which='major', labelsize=plt_set.tick_size)
#ax.xaxis.grid(color='gray', linestyle='dashed')
#ax.yaxis.grid(color='gray', linestyle='dashed')
ax.legend()
ax.set_facecolor('lightgrey')
#ax.set_ylim(0,10)
#plt.close('all')

#fig,ax = plt.subplots(1,1,figsize=(16,8))
#fig.suptitle('Wind Profile Lidar', fontsize=plt_set.fig_size, fontweight="bold")
#
#
#for i in range(len(z)):  
#    ax.plot(df['Date'],  m_z[i,:],alpha = 0.5)   
#ax.set_ylabel("Time", fontsize = plt_set.label)
#ax.set_ylabel('Wind Speed ($ms^-1$)', fontsize = plt_set.label)
##ax.tick_params(axis='both', which='major', labelsize=plt_set.tick_size)
##ax.xaxis.grid(color='gray', linestyle='dashed')
##ax.yaxis.grid(color='gray', linestyle='dashed')
#ax.legend()
#ax.set_facecolor('lightgrey')
#ax.set_ylim(0,10)
##plt.close('all')



