"""
Function that finds the cut containing a galaxy's data

Libraries: import numpy as np
           import pandas / from astropy.table import Table
           from astropy.wcs import WCS 
           from astropy.nddata import Cutout2D
           from astropy.coordinates import SkyCoord
           import astropy.units as u

Inputs: ind -> index of the galaxy inside the table/dataframe
        hdu -> contains all info from file. hdu[1] has the data we are after.
        names -> string array containing the names of the objects contained in the region object
        ra -> string/column of table/dataframe containing ra coordinates in deg for every object
        dec -> string/column of table/dataframe containing dec coordinates in deg for every object
        lim_lrbt -> string/column of table/dataframe containing limits in pixels to perform CutOut2D on hdu    
        
Outputs: cutout -> cutout2D object of the area comprising the object (better if it does not fit tighly so the
                    bkg subs can give account of spatial inhomogeneities)

To improve: -> Give the choice to use ind or not
            -> Give the choice to employ pixels or celestial coordinates. Or any combination of both.

"""
def frame_auto(ind,hdu,names,ra,dec,lim_bot,lim_top,lim_left,lim_right): 
    #-----------------------------------------------------------------------------------------------
    #convert coordinates to pixels
    ra_c = ra[ind]
    dec_c = dec[ind]
    wcs = WCS(hdu[1].header)
    c = SkyCoord(ra_c,dec_c, unit="deg")
    x,y = astropy.wcs.utils.skycoord_to_pixel(c,wcs)

    #-----------------------------------------------------------------------------------------------
    #estimate size of cutout
    position = c
    size_y = (int(y)+int(lim_top[ind]))-(int(y)-int(lim_bot[ind]))
    size_x = (int(x)+int(lim_right[ind]))-(int(x)-int(lim_left[ind]))
    size = (size_y, size_x)     # pixels
    #cutout
    cutout = Cutout2D(hdu[1].data, position, size,wcs=wcs)
    
    return(cutout)


"""
Function that measures the bkgsub signal of a galaxy inside a Cutout2D using photutils aperture_photometry

Packages:   Packages:   import numpy as np
            import pandas as pd
            from astropy.nddata import Cutout2D
            from astropy import units as u
            from astropy.table import Table
            from astropy.io import fits
            from astropy.coordinates import SkyCoord
            from astropy.wcs import WCS
            from astropy.stats import sigma_clipped_stats
            from photutils.aperture import CircularAperture,EllipticalAperture,RectangularAperture
            from photutils.aperture import aperture_photometry,ApertureStats
            from photutils.isophote import Ellipse, EllipseGeometry
            from photutils import DAOStarFinder

Call to other functions: frame_auto(ind,hdu,names,ra,dec,lim_bot,lim_top,lim_left,lim_right)
                         Warning: check table format

Inputs: hdu -> contains all info from file. hdu[1] has the data we are after.
        table -> astropy table object/pandas dataframe. Minimum columns:
                ·table['Name'] : contains objects' names
                ·table['ra'],table['dec'] : contains objects' ra/dec coordinates in deg
                ·table['frame_bottom'],table['frame_top'],table['frame_left'],table['frame_right'] :
                 limits in pixels to perform CutOut2D on hdu
                ·table['x0_pix'], table['y0_pix'], table['sma'], table['ellipticity'], table['tilt'] :
                 parameters of EllipseGeometry object framing the galaxy (consult photutils docs)
        newcol_name -> string contaning the name of the new column that will be created in the table with
                        the background subtracted signal measured in the galaxy
           
Outputs: table -> same table with one extra column containing photometry of galaxy. Photometry performed with
                  photutils aperture_photometry

To improve: -> Give the choice to use ind or not
            -> Give the choice to employ pixels or celestial coordinates. Or any combination of both.
            -> The fit of the galaxy should be performed automatically. Here one needs to provide parameters
                of a predefined ellipse fit (new function on the way!)
            -> Sources inside the galaxy's ellipse are not corrected (source of error!)
            -> Provide errors associated to galaxy fit, bkg sub, SE, photometry (maybe separate this function into 
               many others and combine results afterwards) 
"""
def signal_photutils(hdu,table,newcol_name): 
    
    for ind in range(len(table)):
        #-----------------------------------------------------------------------------------------------
        #AREA AND GALAXY
        cut = frame_auto(ind,hdu,table['Name'],table['ra'],table['dec'],table['frame_bottom'],table['frame_top'],table['frame_left'],table['frame_right'])
        data=cut.data
        wcs=cut.wcs
        
        geometry = EllipseGeometry(x0=table['x0_pix'][ind], y0=table['y0_pix'][ind], sma= table['sma'][ind], eps= table['ellipticity'][ind], pa=table['tilt'][ind] * np.pi / 180.0)
        aper_g = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,geometry.sma * (1 - geometry.eps), geometry.pa)
        
        #Make a mask so the source extraction only occurs in the region outside the gal area
        mask_g = aper_g.to_mask()
        test_mask = mask_g.to_image(data.shape)
        mask_gals = np.ma.make_mask(test_mask)

        #-----------------------------------------------------------------------------------------------
        #SOURCE EXTRACTION
        ##### Calculate image statistics to identify threshold for point sources (using only pixels with positive values):
        mean,median,std=sigma_clipped_stats(data[data>0],sigma=3,maxiters=2,cenfunc=np.ma.mean)
        ##### Calculate image statistics using all pixels to estimate how to replace point sources with background-like value:
        mean2,median2,std2=sigma_clipped_stats(data,sigma=3,maxiters=1,cenfunc=np.ma.mean)

        ##### Find point sources
        daofind=DAOStarFinder(fwhm=3,threshold=mean+1.*std)
        sources=daofind(data,mask = mask_gals)
        #print(sources)

        ##### Create a mask of the point sources:
        positions=[]
        for nbr in range(len(sources)): positions.append((sources['xcentroid'][nbr],sources['ycentroid'][nbr]))

        ap_s=CircularAperture(positions,r=10.)

        #Make a mask so the source extraction only occurs in the region outside the gal area
        #This mask will be employed for the aperture_stats
        maskt = mask_gals
        masks = ap_s.to_mask()
        for i in range(len(ap_s)):
            mask_im = masks[i].to_image(data.shape)
            maskt = maskt + mask_im
        array_mask = np.where(maskt,1,0)
        mask_tot = np.ma.make_mask(array_mask)

        #-----------------------------------------------------------------------------------------------
        #SIGNAL
        #Galaxy signal
        phot_gal = aperture_photometry(data, aper_g)
        # Apply photutils stats to estimate bckg signal
        #SkyAperture of the whole cut containing the galaxy
        positions = (cut.input_position_cutout[0],cut.input_position_cutout[1])
        aper_rect = RectangularAperture(positions,cut.shape[1],cut.shape[0])
        #Estimate background signal with new rect aperture and previously constructed mask
        phot_bkg = aperture_photometry(data,aper_rect,mask_tot)
        phot_stats = ApertureStats(data=data, aperture=aper_rect, mask=mask_tot)
        #Mean signal per pixel
        bkgxpix = phot_stats.mean

        #Save result 
        table[col_name][ind] = phot_gal['aperture_sum'][0]-bkgxpix*aper_g.area_overlap(data)
    
    return(table)


"""
Function that measures the bkgsub signal of a galaxy inside a Cutout2D using numpy sum()

Packages:   import numpy as np
            import pandas as pd
            from astropy.nddata import Cutout2D
            from astropy import units as u
            from astropy.table import Table
            from astropy.io import fits
            from astropy.coordinates import SkyCoord
            from astropy.wcs import WCS
            from astropy.stats import sigma_clipped_stats
            from photutils.aperture import CircularAperture,EllipticalAperture,RectangularAperture
            from photutils.aperture import aperture_photometry,ApertureStats
            from photutils.isophote import Ellipse, EllipseGeometry
            from photutils import DAOStarFinder

Call to other functions: frame_auto(ind,hdu,names,ra,dec,lim_bot,lim_top,lim_left,lim_right)
                         Warning: check table format

Inputs: hdu -> contains all info from file. hdu[1] has the data we are after.
        table -> astropy table object/pandas dataframe. Minimum columns:
                ·table['Name'] : contains objects' names
                ·table['ra'],table['dec'] : contains objects' ra/dec coordinates in deg
                ·table['frame_bottom'],table['frame_top'],table['frame_left'],table['frame_right'] :
                 limits in pixels to perform CutOut2D on hdu
                ·table['x0_pix'], table['y0_pix'], table['sma'], table['ellipticity'], table['tilt'] :
                 parameters of EllipseGeometry object framing the galaxy (consult photutils docs)
        newcol_name -> string contaning the name of the new column that will be created in the table with
                        the background subtracted signal measured in the galaxy
           
Outputs: table -> same table with one extra column containing photometry of galaxy. Photometry performed with
                  numpy sum()

To improve: -> Give the choice to use ind or not
            -> Give the choice to employ pixels or celestial coordinates. Or any combination of both.
            -> The fit of the galaxy should be performed automatically. Here one needs to provide parameters
                of a predefined ellipse fit (new function on the way!)
            -> Sources inside the galaxy's ellipse are not corrected (source of error!)
            -> Provide errors associated to galaxy fit, bkg sub, SE, photometry (maybe separate this function into 
               many others and combine results afterwards) 
            -> Compare to photutils method to provide error and speed!  
"""

def signal_npsum(hdu,table,col_name): 
    
    for ind in range(len(table)):
        #-----------------------------------------------------------------------------------------------
        #AREA AND GALAXY
        cut = frame_auto(ind,hdu,table['Name'],table['ra'],table['dec'],table['frame_bottom'],table['frame_top'],table['frame_left'],table['frame_right'])
        data=cut.data
        wcs=cut.wcs
        #print(np.sum(data))
        
        geometry = EllipseGeometry(x0=table['x0_pix'][ind], y0=table['y0_pix'][ind], sma= table['sma'][ind], eps= table['ellipticity'][ind], pa=table['tilt'][ind] * np.pi / 180.0)
        aper_g = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,geometry.sma * (1 - geometry.eps), geometry.pa)
        
        #Make a mask so the source extraction only occurs in the region outside the gal area
        mask_g = aper_g.to_mask()
        #mask_g.data.any()
        test_im = mask_g.to_image(data.shape)
        mask_gals = np.ma.make_mask(test_im)

        #-----------------------------------------------------------------------------------------------
        #SOURCE EXTRACTION
        ##### Calculate image statistics to identify threshold for point sources (using only pixels with positive values):
        mean,median,std=sigma_clipped_stats(data[data>0],sigma=3,maxiters=2,cenfunc=np.ma.mean)
        ##### Calculate image statistics using all pixels to estimate how to replace point sources with background-like value:
        mean2,median2,std2=sigma_clipped_stats(data,sigma=3,maxiters=1,cenfunc=np.ma.mean)

        ##### Find point sources
        daofind=DAOStarFinder(fwhm=3,threshold=mean+1.*std)
        sources=daofind(data,mask = mask_gals)
        #print(sources)

        ##### Create a mask of the point sources:
        positions=[]
        for nbr in range(len(sources)): positions.append((sources['xcentroid'][nbr],sources['ycentroid'][nbr]))

        ap_s=CircularAperture(positions,r=10.)

        #Make a mask so the source extraction only occurs in the region outside the gal area
        #This mask will be employed for the aperture_stats
        maskt = mask_gals
        masks = ap_s.to_mask()
        for i in range(len(ap_s)):
            mask_im = masks[i].to_image(data.shape)
            maskt = maskt + mask_im
        array_mask = np.where(maskt,1,0)
        mask_tot = np.ma.make_mask(array_mask)

        #-----------------------------------------------------------------------------------------------
        #SIGNAL
        #Galaxy signal
        mask_gals = np.ma.make_mask(test_im)
        mask_nogals = mask_tot == False
        #Number of pixels
        pix_g = np.sum(mask_gals)
        pix_g2 = np.sum(mask_nogals)
        pix_g,pix_g2
        #Bkg sub
        sign_nogal = cut.data[mask_nogals].sum()  #L(Halpha)
        new_data = data[mask_gals] - np.ones(cut.data[mask_gals].shape)*sign_nogal/pix_g2
        sign_gal = new_data.sum() 

        #Save result 
        #print(sign_gal)
        table[col_name][ind] = sign_gal
    
    return(table)


"""
Function that transforms the signal in 2 optical bands into magnitudes and then 2-1 colour.

Libraries: import numpy as np
           import pandas / from astropy.table import Table
           import astropy.units as u, astropy.constants as c

Inputs: sgn_b1,sgn_b2 -> array with signal in band 1;band 2. Units in erg/cm2sA
        zptb1,zptb2 -> zeropoint for band1;band2 (no units)
        zeromag1,zeromag2 -> AB magnitude zeropoint for band1;band2 (no units)
        wdth1,wdth2 -> filter width for band1;band2 
        wv1,wv2 -> wavelenght of the band1;band2
         
Outputs: colour -> b2-b1 colour magnitude (array)

To improve: -> Use astropy units!
            -> Include extinction
            -> Values do not match expectations. Maybe something more is missing...
            -> Study log units: https://docs.astropy.org/en/stable/units/logarithmic_units.html
"""
def colourb2_b1(sgn_b1,sgn_b2,zptb1,zptb2,wdth1,wdth2,zeromag1,zeromag2,wv1,wv2): 
    #-----------------------------------------------------------------------------------------------
    #Some lazy variables -_-
    Jy2erg=1e-23
    #c_in_m_s = 3E8
    c_in_A_s = 2.998E18
    #pc2cm = 3e18
    #erg2jy = 1e-23
    #as_kpc = .278 
    #Mpc2pc= 1e6
    #--------------------------------------------------------------------------------
    """
    #Example for Hydra cluster
        ZeroMag_AB_2 = 3631  #Jy  #zeromag1,zeromag2
        zpt = 27                  #zpt1,zpt2
        ha_filter_width_A = 170   #wdth
        u_filter_width_A = 880    #wdth
        r_filter_width_A = 1276   #wdth
        ha_line = 6563            #wv
        r_line = 7000             #wv
        u_line = 3650             #wv
        fl0_u = 4.27E-9 #erg/cm2sA
        fl0_r = 1.74E-9 #erg/cm2sA
    """

    #-----------------------------------------------------------------------------------------------
    #Flux from signal
    #Sum of the signal: photons/(s*pix)
    #Signal dimension: erg/(A*s*cm^2)
    flux1 = (sgn_b1*10**((-zptb1)/2.5)*zeromag1*Jy2erg*c_in_A_s/(wv1**2)*wdth1)  
    flux2 = (sgn_b2*10**((-zptb2)/2.5)*zeromag2*Jy2erg*c_in_A_s/(wv2**2)*wdth2)  

    #-----------------------------------------------------------------------------------------------
    #Magnitude from flux density(erg/(s*cm^2))
    colour = 2.5*(np.log10(flux2)-np.log10(flux2))
    
    return(colour)


"""
Function that plots a colour-stellar mass diagram using plotly.

Libraries: import numpy as np
           import pandas / from astropy.table import Table
           import astropy.units as u, astropy.constants as c

Inputs: x -> (array) stellar mass for each object
        y -> (array) colour for each object
         
Outputs: no output. prints a CM map

To improve: -> Save?
            -> Include extinction
            -> Values do not match expectations. Maybe something more is missing...
            -> Study log units: https://docs.astropy.org/en/stable/units/logarithmic_units.html
"""

def colmass_plot(x,y)

    fig = go.Figure(go.Histogram2dContour(
            x = x, y = y,
            colorscale = 'GnBu',
            contours_showlabels = True,
            #labelfont = dict(size = 12,color = 'white'),
    ))

    fig.add_trace(go.Scatter(
            x = x,y = y,
            xaxis = 'x',yaxis = 'y',
            mode = 'markers',
            marker = dict(symbol="diamond",color = 'crimson',size = 10,line=dict(width=2, color="black")),
        ))

    fig.update_layout(
        xaxis_title="log Ms",
        yaxis_title="u-r colour",
        title = 'u-r colour mass diagram'
    )
    fig.show()
    #Save png (needs to import kaleido)
    #fig.write_image("yourfile.png") 
    #Save html (needs to import plotly as ply)
    #plotly.offline.plot(fig, filename='C:/Users/Clara/Desktop/08_02/ur_shift.html')
