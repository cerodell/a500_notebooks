class constants:
    """ 
   A class containing constants relevant to atmospheric science.

   Parameters
   ----------

   none -- class variables only

   References
   ----------
  
   Emanuel appendix 2
   
   """
    pa2hPa=1.e-2
    hPa2pa=1.e2
    Tp = 273.16
    #triple point in kelvin
    Tc = 273.15
    #0 deg C in Kelvin
    eps = 0.622
    # Rd/Rv
    phi0 = 9155.80612
    #lv0/Tp  (J/kg/K
    es0 = 6.112 #hPa
    p0 = 1.e5  #Pa
    lv0 = 2.501e6  #J/kg
    # J/kg
    Rv = 461.50
    # J/kg/K
    Rd = 287.04
    # J/kg/K
    cpv = 1870.
    # Heat capacity of water vapor (J/kg/K)
    cl = 4190.
    # Heat capacity of liquid water (J/kg/K)
    cpd = 1005.7
    # Heat capacity of dry air (J/kg/K)
    delta = (1 - eps)/eps  #Thompkins 2.33 for Tv
    g0 = 9.8
    # m/s^2
    D = 2.36e-5
    # Diffusivity m^2/s^1
    # Note: fairly strong function of temperature
    #       and pressure -- this is at 100kPa, 10degC
    rhol = 1000.

    
