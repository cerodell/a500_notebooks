#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 15:49:14 2019

@author: rodell
"""

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd

import metpy.calc as mpcalc
from metpy.plots import Hodograph, SkewT
from metpy.units import units

import context

#user      = 'rodell'

file = 'sounding/KDNR_20160831_00UTC_sharppy.txt'
#filein  = '/Users/'+user+'/Google Drive File Stream/Shared drives/Research/CRodell/BL/Data/'
#save    = '/Users/'+user+'/Google Drive File Stream/Shared drives/Research/CRodell/BL/Images/'



ylabel     = 14
fig_size   = 20
tick_size  = 12
title_size = 16

#df_sonde = pd.read_csv(filein+file)


"######################  Read File ######################"
df_sonde = pd.read_csv(context.pro_data_dir / file).replace(-9999, np.NaN, regex=True)

#df_sonde = df_sonde.replace(r'  ', np.NaN, regex=True)
##Make a list of foats of each variable 
height, press = list(df_sonde['HGHT']), list(df_sonde['LEVEL'])
temp, dew = list(df_sonde['TEMP']), list(df_sonde['DWPT'])
wsp, wdir = list(df_sonde['WSPD']), list(df_sonde['WDIR'])



"#################Solve for Dew Point ###################"
##Solve for Dew Point (source: http://irtfweb.ifa.hawaii.edu/~tcs3/tcs3/Misc/Dewpoint_Calculation_Humidity_Sensor_E.pdf)
##constant varibles needed
beta    = 17.62
lambda_ = 243.15 #deg C

#dew =[]   #initiate list for Dew-Point
#
#for i in range(len(temp)):
#    dew_i = round((lambda_*((np.log(rh[i]/100))+((beta*temp[i])/(lambda_+temp[i]))))\
#                  /(beta-(np.log((rh[i]/100)+((beta*temp[i])/(lambda_+temp[i]))))),3)
#    dew.append(dew_i)



"#################Solve for u, and LCL ###################"
##Add units so Metpy can use data 
index = 3
tempC, dew, the_press = temp*units.degC, dew*units.degC, press*units.hPa
wsp, wdir = wsp*units.knots, wdir*units.degrees

##Calcualte u and v
u, v = mpcalc.wind_components(wsp, wdir)

#print(tempC[index:])
##Calculate the LCL
lcl_pressure, lcl_temperature = mpcalc.lcl(the_press[index], tempC[index], dew[index])

#print(lcl_pressure, lcl_temperature)

##Calculate the parcel profile.
parcel_prof = mpcalc.parcel_profile(the_press, tempC[index], dew[index]).to('degC')



"###################### Make Plot ######################"

##Create a new figure. The dimensions here give a good aspect ratio
fig = plt.figure(figsize=(9, 9))
fig.suptitle('Plecian Mountian \nMay 11th 1500 MDT', fontsize=16, x=0.125, horizontalalignment='left', \
             fontweight="bold")
skew = SkewT(fig, rotation=30)

##Plot the data using normal plotting functions, in this case using
##log scaling in Y, as dictated by the typical meteorological plot
skew.plot(the_press, tempC, 'r',linewidth=2.5)
skew.plot(the_press, dew, 'g', linewidth=2.5)
skew.plot_barbs(the_press[0::5], u[0::5], v[0::5])
#skew.ax.set_ylim(1000, 50)
#skew.ax.set_xlim(-40, 60)

##Plot LCL as black dot
skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')

##Plot the parcel profile as a black line
skew.plot(the_press, parcel_prof, 'k', linewidth=2)

##Shade areas of CAPE and CIN
skew.shade_cin(the_press, tempC, parcel_prof)
#skew.shade_cape(press, temp, parcel_prof)

##Plot a zero degree isotherm
skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)

##Add the relevant special lines
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()

##Create a hodograph
##Create an inset axes object that is 40% width and height of the
##figure and put it in the upper right hand corner.
#ax_hod = inset_axes(skew.ax, '40%', '40%', loc=1)
#h = Hodograph(ax_hod, component_range=50.)
#h.add_grid(increment=20)
#h.plot_colormapped(u, v, wsp)  # Plot a line colored by wind speed

# Show the plot
plt.show()
#fig.savefig(save +'skew_t_zoom')
























