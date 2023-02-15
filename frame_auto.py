"""
Function that finds the cut containing a galaxy's data employing Cutout2D.

Libraries: import numpy
           import pandas / from astropy.table import Table
           from astropy.wcs import WCS 
           from astropy.nddata import Cutout2D
           from astropy.coordinates import SkyCoord
           import astropy.units as u

Inputs: hdu -> contains all info from file. hdu[1] has the data we are after
        names -> string array containing the names of the bodies contained in the region object
        ind -> index of the galaxy inside the region file
        ra,dec (deg) -> coordinates of the centre of the cut
        lim_bot,lim_top,lim_left,lim_right -> number of pixels in each direction from the centre
        
        
Outputs: cutout2D -> cutout object of the area

To improve: -> Give the choice to use ind or not
            -> Give the choice to employ pixels or celestial coordinates. Or any combination of both.
"""
def frame_auto(ind,hdu,names,ra,dec,lim_bot,lim_top,lim_left,lim_right): 

    #convert coordinates to pixels
    ra_c = ra[ind]
    dec_c = dec[ind]
    wcs = WCS(hdu_r[1].header)
    c = SkyCoord(ra_c,dec_c, unit="deg")
    x,y = astropy.wcs.utils.skycoord_to_pixel(c,wcs)

    #estimate size of cutout
    position = c
    size_y = (int(y)+int(lim_top[ind]))-(int(y)-int(lim_bot[ind]))
    size_x = (int(x)+int(lim_right[ind]))-(int(x)-int(lim_left[ind]))
    size = (size_y, size_x)     # pixels
    #cutout
    cutout = Cutout2D(hdu_r[1].data, position, size,wcs=wcs)
    
    
    return(cutout)