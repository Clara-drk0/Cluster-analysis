"""
Function that estimates the GF (in 1/M_sun) from the sum of the selected images (in counts/sec).
It employs a parameterization I took from somewhere :(

Libraries: import astropy
           import numpy as np
           import pandas as pd/ from astropy.table import Table

Inputs: hdu_path -> path to the fits file containing the data 
        HIfile_path -> path to the .txt containing the signal/flux? estimation for each source
          
Outputs: sources -> astropy table with HIfile_path data+ new GF column for each source

To improve: -> Aaaaaggh did this a long time ago...
            -> Find the notebook where I noted the physical dimension of measurements
                I think the ones provided from MK HI table are wrong :( 
            -> Make more general!
"""

def HI_mass(hdu_path, HIfile_path):

    #MeerKAT HI
    #--------------------------------------------------------------------------------------
    #Define cluster distance
    distance = 58.6

    #--------------------------------------------------------------------------------------
    #Number of pixels per beam
    hdu = fits.open(hdu_path)
    hi_cellsize = hdu[0].header['CDELT2'] * 3600. * u.arcsec
    bmajor = hdu[0].header['BMAJ']
    bminor = hdu[0].header['BMIN']
    pix_per_beam = bmajor / hi_cellsize * bminor / hi_cellsize * np.pi / (4 * np.log(2))
   
    #--------------------------------------------------------------------------------------
    #signal sum/flux??? from HI source catalog
    sources =  astropy.io.ascii.read(HIfile_path,data_start=0, delimiter=" ")
    f_sum = sources['col15'] #mJy?

    #--------------------------------------------------------------------------------------
    # Signal in JyHz
    SJyHz = f_sum / pix_per_beam

    #--------------------------------------------------------------------------------------
    #Estimate masses (where is this parameterization from?Kennicut?)
    M_HI = 49.7 * SJyHz * distance**2
    sources['M_HI'] = M_HI

    return(sources)