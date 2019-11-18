"""
model a bulk entraining plume with constant entrainment rate
"""
import numpy as np
import pandas as pd
from pprint import pformat
from .constants import constants as c
from .thermlib import find_Tmoist,find_thetaep,find_rsat,tinvert_thetae
from scipy.interpolate import interp1d
from ..soundings.wyominglib import write_soundings, read_soundings
import json


from scipy.integrate import ode
import matplotlib.pyplot as plt
from ..skewT.nudge import nudge

def derivs(t, y, entrain_rate, interpTenv, interpTdEnv, interpPress):
    """Function that computes derivative vector for ode integrator
       see http://clouds.eos.ubc.ca/~phil/courses/atsc405/docs/entrain.pdf for equations

    Parameters
    ----------
    
    t: float
       time (s)
    y: vector
       4-vector containing wvel (m/s), height (m), thetae (K), rT (kg/kg)
    entrain_rate: float
                  1/m dm/dt (s-1)
    interpTenv: func
                interp1d function for environmental temperature T(z) 
    interpTdEnv: func
                interp1d function for environmental dewpoint temperature Td(z)
    interpPress: func
                interp1d function for presusure  p(z)

    Returns
    -------

    yp: vector
       4-vector containing time derivatives of wvel (m/s^2), height (m/s), thetae (K/s), rT (kg/kg/s)
    """
    yp = np.zeros((4,1))
    velocity = y[0]
    height = y[1]
    thetae_cloud = y[2]
    rT_cloud = y[3]
    #yp[0] is the acceleration, in this case the buoyancy 
    yp[0] = calcBuoy(height, thetae_cloud, interpTenv, interpTdEnv, interpPress)
    press = interpPress(height)*100. #Pa
    Tdenv = interpTdEnv(height) + c.Tc #K
    Tenv = interpTenv(height) + c.Tc #K
    rTenv = find_rsat(Tdenv, press) #kg/kg
    thetaeEnv = find_thetaep(Tdenv, Tenv, press)
    #yp[1] is the rate of change of height
    yp[1] = velocity
    #yp[2] is the rate of change of thetae_cloud
    yp[2] = entrain_rate*(thetaeEnv - thetae_cloud)
    #yp[3] is the rate of change of rT_cloud
    yp[3] = entrain_rate*(rTenv - rT_cloud)
    return yp

def calcBuoy(height, thetae0, interpTenv, interpTdEnv, interpPress):
    """function to calculate buoyant acceleration for an ascending saturated parcel
       this version neglects liquid water loading
    
    Parameters
    ----------
    
    height: float
            parcel height (m)
    thetae0: float
            parcel thetae (K)

    interpTenv: func
                interp1d function for environmental temperature T(z) 
    interpTdEnv: func
                interp1d function for environmental dewpoint temperature Td(z)
    interpPress: func
                interp1d function for presusure  p(z)

    Returns
    -------

    buoy: float
          buoyancy (m/s/s)
    """
    #input: height (m), thetae0 (K), plus function handles for
    #T,Td, press soundings
    #output: Bout = buoyant acceleration in m/s^2
    #neglect liquid water loading in the virtual temperature
    
    press=interpPress(height)*100.#%Pa
    Tcloud=find_Tmoist(thetae0,press) #K
    rvcloud=find_rsat(Tcloud,press); #kg/kg
    Tvcloud=Tcloud*(1. + c.eps*rvcloud)
    Tenv=interpTenv(height) + c.Tc
    Tdenv=interpTdEnv(height) + c.Tc
    rvenv=find_rsat(Tdenv,press); #kg/kg
    Tvenv=Tenv*(1. + c.eps*rvenv)
    TvDiff=Tvcloud - Tvenv
    buoy = c.g0*(TvDiff/Tvenv)
    return buoy

def read_sounding(filename):
    """given an h5 file from wyominglib return a dataframe plus a units dictionary

    Parameters
    ----------
    filename: str
             filename of the h5 sounding filename

    Returns
    ------
    df_sounding: pandas dataframe 
               : sounding columns are temperature (deg C), dewpoint (deg C), height (m), press (hPa)

    units_dict : dict
               : units for each column 
    """
    print('reading file: %s\n' %filename)
    attributes={}
    with h5py.File(filename,'r') as f:
        keys=f.attrs.keys()
        for key in keys:
            try:
                attributes[key]=f.attrs[key]
            except IOError:
                print('empty key: ',key)
    print('\nread in these attributes: \n\n',pformat(attributes))
    sounding_dict={}
    with pd.HDFStore(filename,'r') as store:
        times=store.keys()
        for the_time in times:
            sounding_dict[the_time]=store[the_time]
    df_sounding=sounding_dict[times[3]]
    title_string=attributes['header']
    index=title_string.find(' Observations at')
    location=title_string[:index]
    title='{} at {}'.format(location,times[0][2:])
    print('title: :',title)
    units=attributes['units'].split(';')
    units_dict={}
    for count,var in enumerate(df_sounding.columns):
        units_dict[var]=units[count]
    return df_sounding,units_dict
    
def integ_entrain(df_sounding,entrain_rate):
    """integrate an ascending parcel given a constant entrainment rate
       this version hardwired to start parcel at 800 hPa with cloud base
       values of environment at 900 hPa

    Parameters
    ----------

    df_sounding: pandas dataframe 
               : cloumns are temperature, dewpoint, height, press

    entrain_rate: float
                  1/m dm/dt (s-1)

    Returns
    -------

    df_out: dataframe
          dataframe containing wvel (m/s) ,cloud_height (m) , thetae (K), rT (kg/kg) for assending parcel

   interpPress: func
              interp1d function for presusure  p(z) (used for plotting)
    """
    press = df_sounding['pres'].values
    height = df_sounding['hght'].values
    temp = df_sounding['temp'].values
    dewpoint = df_sounding['dwpt'].values
    envHeight= nudge(height)

    interpTenv = interp1d(envHeight,temp)
    interpTdEnv = interp1d(envHeight,dewpoint)
    interpPress = interp1d(envHeight,press)
    #
    # call this cloudbase
    #
    p900_level = len(press) - np.searchsorted(press[::-1],900.)
    thetaeVal=find_thetaep(dewpoint[p900_level] + c.Tc,temp[p900_level] + c.Tc,press[p900_level]*100.)
    rTcloud = find_rsat(dewpoint[p900_level] + c.Tc, press[p900_level]*100.)
    #
    # start parcel here
    #
    p800_level = len(press) - np.searchsorted(press[::-1],800.)
    height_800=height[p800_level]
    winit = 0.5 #initial velocity (m/s)
    yinit = [winit, height_800, thetaeVal, rTcloud]  
    tinit = 0  #seconds
    tfin = 2500  #seconds
    dt = 10   #seconds

    #want to integrate derivs using dopr15 runge kutta described at
    # http://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html
    #
    r = ode(derivs).set_integrator('dopri5')
    r.set_f_params(entrain_rate, interpTenv, interpTdEnv, interpPress)
    r.set_initial_value(yinit, tinit)
    
    #the while loop below  integrates every dt seconds
    #we stop tracking the parcel when the time runs out, or if the parcel stops moving/is desecnding
    #
    var_out = []
    time_out =[]
    while r.successful() and r.t < tfin and r.y[0] > 0:
        #find y at the next time step
        #(r.integrate(t) updates the fields r.y and r.t so that r.y = integral of derivs(t) and r.t = time 
        #where derivs is a vector with the variables to be integrated
        #
        # move ahead by dt
        #
        r.integrate(r.t+dt)
        #
        # stop if there is negative vertical velocity
        #
        if r.y[0] <= 0:
            break
        #
        #save values for dataframe
        #
        var_out.append(r.y)
        time_out.append(r.t)
    #
    # convert the output into a datafram
    #
    colnames=['wvel','cloud_height','thetae_cloud','rT_cloud']
    df_out=pd.DataFrame.from_records(var_out,columns=colnames)
    df_out['time'] = time_out
    return df_out,interpPress

def make_plot(df_result,df_sounding,interpPress):
    """ return two plots with the vertical velocity and temperature profiles

    Parameters
    ----------
    df_result: dataframe
          dataframe containing wvel (m/s) ,cloud_height (m) , thetae (K), rT (kg/kg) for assending parcel

    df_sounding: pandas dataframe 
               : sounding columns are temperature (deg C), dewpoint (deg C), height (m), press (hPa)

    interpPress: func
            :    interp1d function for presusure  p(z)

    Returns
    -------
    ax1,ax2: matplotlib axes
           : plotting axis for the two figures
    """
    plt.style.use('ggplot')
    fig,ax1 = plt.subplots(1,1)
    ax1.plot(df_result['wvel'], df_result['cloud_height'], label='wvel')
    ax1.set_xlabel('vertical velocity (m/s)')
    ax1.set_ylabel('height above surface (m)')


    Tcloud = np.zeros_like(df_result['cloud_height'])
    rvCloud = np.zeros_like(df_result['cloud_height'])
    rlCloud = np.zeros_like(df_result['cloud_height'])

    for i,row in enumerate(df_result.iterrows()):
        the_press = interpPress(row[1]['cloud_height'])*100.
        Tcloud[i], rvCloud[i], rlCloud[i] = tinvert_thetae(row[1]['thetae_cloud'], 
                                                        row[1]['rT_cloud'], the_press)
        
    Tadia= np.zeros_like(df_result['cloud_height'])
    rvAdia = np.zeros_like(df_result['cloud_height'])
    rlAdia = np.zeros_like(df_result['cloud_height'])
    

    heights=df_result['cloud_height']
    thetae_cloud0 = df_result.loc[0,'thetae_cloud']
    rT_cloud0 = df_result.loc[0,'rT_cloud']
    for i, height in enumerate(heights):
        the_press = interpPress(height)*100.
        Tadia[i], rvAdia[i], rlAdia[i] = tinvert_thetae(thetae_cloud0, 
                                                      rT_cloud0, the_press)
    
    fig,ax2=plt.subplots(1,1)
    ax2.plot(Tcloud - c.Tc, df_result['cloud_height'], 'r-',label='cloud')
    height = df_sounding['hght'].values
    temp = df_sounding['temp'].values
    ax2.plot(temp, height, 'g-',label='environment')
    ax2.plot(Tadia - c.Tc, df_result['cloud_height'], 'b-',label='moist aidabat')
    ax2.set_xlabel('temperature (deg C)')
    ax2.set_ylabel('height above surface (m)')
    ax2.legend(loc='best')
    return ax1,ax2

def get_sounding(get_data=True,soundingfile='linttlerock.csv',metafile='metadata.json'):
    if get_data:
        values=dict(region='naconf',year='2012',month='7',start='0100',stop='3000',station='72340')
        write_soundings(values, 'littlerock')
        soundings= read_soundings('littlerock')
        the_time=(2012,7,17,0)
        sounding=soundings['sounding_dict'][the_time]
        sounding.to_csv(soundingfile)
        title_string=soundings['attributes']['header']
        index=title_string.find(' Observations at')
        location=title_string[:index]
        print(f'location: {location}')
        units=soundings['attributes']['units'].split(';')
        units_dict={}
        for count,var in enumerate(sounding.columns[1:]):
            units_dict[var]=units[count]
        soundings['attributes']['location']=location
        soundings['attributes']['units']=units
        soundings['attributes']['time']=the_time
        soundings['attributes']['history']='written by mixing_line_calc'
        with open(metafile,'w') as outfile:
            json.dump(soundings['attributes'],outfile,indent=4)
        attributes=soundings['attributes']
    else:
        sounding = pd.read_csv(soundingfile)
        with open(metafile,'r') as infile:
            attributes=json.load(infile)
    return sounding,attributes


if __name__ == "__main__":
    #
    # get the location of the sounding_dir so we can find the littlerock sounding
    # use the pathlib module so that this works on osx, linux and windows
    #
    df_sounding,attributes = get_sounding()
    #
    # construct the full path to the file using the os dependent syntax
    #
    entrain_rate = 2.e-4
    df_result, interpPress=integ_entrain(df_sounding,entrain_rate)
    ax1, ax2 = make_plot(df_result,df_sounding,interpPress)
    ax1.set_title(('vertical velocity of a cloud parcel vs height' '\n'
                   'entrainment rate of {:4.1e} $s^{{-1}}$').format(entrain_rate))
    ax2.set_title(('temp. of rising cloud parcel vs height' '\n'
                        'entrainment rate of {:4.1e} $s^{{-1}}$').format(entrain_rate))
    ax2.set(xlim=[-50,25],ylim=[0,15000])
    plt.show()









