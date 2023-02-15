"""
Function that estimates the SFR (in M_sun/yr) from the sum of the selected images (in counts/sec).
It employs Kennicutt98 parameterization

Libraries: import numpy as np
           import pandas as pd/ From astropy.table import Table

Inputs: sum_signal -> array, table/dataframe column, numpy float,... with the signal in Ha to convert.
                      Signal as estimated by function It needs a column with their position in SkyCoord
          
Outputs: SFR_Ha -> original table + 'SFR' column

To improve: -> Introduce parameters as inputs
            -> build a table with all the data 
            -> allow to choose different parameterizations for different emissions.
            -> Connect to the Kitchen appliances and learn to make coffee
"""

def SFRates(sum_signal):
    #-------------------------------------------------------------------------------
    #Parameters
    pix_scale = 0.333
    #npix = 1.6E5  #Number of pixels in NGC 3312-> This would be our table1['Sum']
    dust_corr_mwg = 1.15  #Galactic extinction from milky way gas (NED)
    dhc = 58.6  #Mpc
    abmag_zpt = 48.6    #-2.5*np.log10(3731E-23)=48.57043687739793 
    ZeroMag_AB = 3731  #Jy
    ZeroMag_AB_2 = 3631  #Jy
    zpt = 27
    k = 7.9E-42
    lamb = 6563 #A
    filter_width_A = 145.4
    ha_line = 6563
    #--------------------------------------------------------------------------------
    Jy2erg=1e-23
    c_in_m_s = 3E8
    c_in_A_s = 2.998E18
    pc2cm = 3e18
    erg2jy = 1e-23
    as_kpc = .278 
    Mpc2pc= 1e6
    
    photons = u.def_unit('photons')

    #Sum of the signal: photons/(s*pix)
    #Flux density from magnitude: erg/(s*cm^2)
    flux_Ha = (sum_signal*10**((-zpt)/2.5)*ZeroMag_AB_2*Jy2erg*c_in_A_s/(ha_line**2)*filter_width_A)  
    #Luminosity from flux density: erg/s
    Lum = flux_Ha*4*np.pi*(dhc*pc2cm*Mpc2pc)**2 
    #Star formation rates from Luminosity (Kennicut 98+dust correction): M_sun/yr
    SFR_Ha = (Lum*7.9E-42)*dust_corr_mwg
    
    return(SFR_Ha)