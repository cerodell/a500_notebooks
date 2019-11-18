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
