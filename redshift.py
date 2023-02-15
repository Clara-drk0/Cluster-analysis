"""
Function that reads a region object extracted using region on a ds9 region file.

Libraries:  import astropy.units as u
            import numpy as np
            from astropy.constants import c

Inputs: velocities -> array  with velocities in km/s
        
Outputs: z -> redshift

To improve -> ??? 
"""

def redshift(velocities):

    #Convert speed of light from m/s to km/s
    c_kms = c.to(u.km/u.s) 

    #Redshift
    z = (np.sqrt(1+velocities/c_kms)/np.sqrt(1-velocities/c_kms))-1

    return(z)