#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 14:58:20 2019

@author: crodell
"""


import urllib,os
import context
import bl_pro
import pandas as pd
from pathlib import Path
from netCDF4 import Dataset

#from bs4 import BeautifulSoup
#import requests





#from  bl_pro.utils.data_read import download
the_file='WC16_XPIA/WLS866-16_2015_06_07__00_00_00.rtd'
#out=download(the_file,dest_folder=bl_pro.lidar_dir)
    
#####################################################################
##Read files
#################################################################### 

case_name= bl_pro.lidar_dir / the_file


feilds = ['Timestamp','Position','Temperature','Wiper Count','Alpha','Beta','Gamma','40m CNR (dB)',
     '40m Radial Wind Speed (m/s)','40m Radial Wind Speed Dispersion (m/s)','40m Wind Speed (m/s)',
     '40m Wind Direction','40m X-wind (m/s)','40m Y-wind (m/s)','40m Z-wind (m/s)',
     '50m CNR (dB)','50m Radial Wind Speed (m/s)','50m Radial Wind Speed Dispersion (m/s)',
     '50m Wind Speed (m/s)','50m Wind Direction','50m X-wind (m/s)','50m Y-wind (m/s)',
     '50m Z-wind (m/s)','60m CNR (dB)','60m Radial Wind Speed (m/s)',
     '60m Radial Wind Speed Dispersion (m/s)','60m Wind Speed (m/s)','60m Wind Direction',
     '60m X-wind (m/s)','60m Y-wind (m/s)','60m Z-wind (m/s)','80m CNR (dB)','80m Radial Wind Speed (m/s)',
     '80m Radial Wind Speed Dispersion (m/s)','80m Wind Speed (m/s)',
     '80m Wind Direction','80m X-wind (m/s)','80m Y-wind (m/s)', '80m Z-wind (m/s)','100m CNR (dB)','100m Radial Wind Speed (m/s)',
     '100m Radial Wind Speed Dispersion (m/s)','100m Wind Speed (m/s)',
     '100m Wind Direction','100m X-wind (m/s)','100m Y-wind (m/s)','100m Z-wind (m/s)','120m CNR (dB)','120m Radial Wind Speed (m/s)',
     '120m Radial Wind Speed Dispersion (m/s)','120m Wind Speed (m/s)',
     '120m Wind Direction','120m X-wind (m/s)','120m Y-wind (m/s),120m Z-wind (m/s)',
     '140m CNR (dB)','140m Radial Wind Speed (m/s)','140m Radial Wind Speed Dispersion (m/s)',
     '140m Wind Speed (m/s)','140m Wind Direction','140m X-wind (m/s)','140m Y-wind (m/s)',
     '140m Z-wind (m/s)','150m CNR (dB),150m Radial Wind Speed (m/s)','150m Radial Wind Speed Dispersion (m/s)',
     '150m Wind Speed (m/s)','150m Wind Direction','150m X-wind (m/s)','150m Y-wind (m/s)','150m Z-wind (m/s)',
     '160m CNR (dB)','160m Radial Wind Speed (m/s)','160m Radial Wind Speed Dispersion (m/s)','160m Wind Speed (m/s)',
     '160m Wind Direction','160m X-wind (m/s)','160m Y-wind (m/s)','160m Z-wind (m/s),180m CNR (dB)','180m Radial Wind Speed (m/s)',
     '180m Radial Wind Speed Dispersion (m/s)','180m Wind Speed (m/s)','180m Wind Direction,180m X-wind (m/s)','180m Y-wind (m/s)',
     '180m Z-wind (m/s)','200m CNR (dB)'',200m Radial Wind Speed (m/s)','200m Radial Wind Speed Dispersion (m/s)','200m Wind Speed (m/s)',
     '200m Wind Direction','200m X-wind (m/s)','200m Y-wind (m/s)','200m Z-wind (m/s)']

df = pd.read_csv(case_name, encoding ='ISO-8859-1',  skiprows=range(0, 41), sep = '\s+' , names=feilds)
#lines  = handle.readlines()
#
#idx = 12
#for line in lines:
##    entries = line.decode("utf-8").split('\t') 
#    entries = line.split('\t')
#    if(entries[0][0]=='D'): 
#        print(entries[idx])
    

