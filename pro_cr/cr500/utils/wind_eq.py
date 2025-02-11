from cr500.utils import constants
import numpy as np

def RxL(w_str,u_str,z,z_i,m_bl):  
    '''
    This Function is used for the wind profile in the Radix Layer

    '''
    dim = (1/constants.C)*(z/z_i)*((w_str/u_str)**constants.B)
#    if dim > 1:
#       m_z = m_bl
#   else:
    m_z = m_bl * ((dim**constants.D)**constants.A) * np.exp(constants.A*(1-dim**constants.D))
    m_z = np.array(m_z)
    return m_z



def w_str(thetav_avg,z_i,wtheta_peturb):
    '''
    This function solves for the Free-convection scaling velocity
    (aka Deardorff Velocity)
    '''
    w_str = ((constants.g/thetav_avg)*z_i*wtheta_peturb)**(1/3)
    w_str = np.array(w_str)
    return w_str




def u_str(uw_peturb,vw_peturb):
    '''
    This function solves for the Friction velocity
    '''
    u_str = (uw_peturb**2 + vw_peturb**2)**(1/4)
    u_str = np.array(u_str)
    return u_str

def do_reynolds(var):
    """
    Do time average and pertuabtion of declared variable
    """
    avg= np.mean(var)
    perturb= np.array(var) - avg
    return avg,perturb


def logwind_neutral(u_str,z):
    m_z = (u_str / constants.k) * np.log(z / constants.z_o)
    return (m_z)


def loglin_stable(u_str,z,L):
    m_z = (u_str / constants.k) * (np.log(z / constants.z_o) + (6 * (z / L)))
    return(m_z)



def obukhov_len(u_str,temp,wtheta_peturb):
    L = (- u_str ** 3) / (constants.k * (constants.g / temp) * wtheta_peturb)
    return(L)




