"""
Function that creates 2 new columns in a table with the distances from the cluster centre
to the galaxies in deg and rad.

Libraries: import astropy
           from astropy.table import Table/ import pandas as pd
           from astropy.coordinates import SkyCoord
           import numpy as np

Inputs: table_source -> table object (astropy or dataframe) that contains all the data from the galaxies.
                        It needs 2 columns with their position in ra,dec (deg)
          
Outputs: table_source -> original table + 'dist from centre' columns

To improve: ???
"""

def dist_ctr(table_source):
    
    #-----------------------------------------------------------------------------------
    #SkyCoord column
    coord = SkyCoord(table_source['ra']*u.deg,table_source['dec']*u.deg)
    #coord = SkyCoord(table_source['RA'],table_source['DEC'])
    
    #-----------------------------------------------------------------------------------
    #Angular distance to the centre
    hydra_ctr = SkyCoord.from_name('Hydra I Cluster')

    #-----------------------------------------------------------------------------------
    table_source['SkyCoord'] = SkyCoord(table_source['ra']*u.deg,table_source['dec']*u.deg)
    #table_source['SkyCoord'] = SkyCoord(table_source['RA'],table_source['DEC'])
    table_source['Dst_to_ctr_deg'] = (hydra_ctr.separation(table_source['SkyCoord']))*u.deg
    table_source['Dst_to_ctr_rad'] = np.zeros(len( table_source['Dst_to_ctr_deg']))
    for i in range(len(table_source['Dst_to_ctr_deg'])):
        table_source['Dst_to_ctr_rad'][i] = math.radians(table_source['Dst_to_ctr_deg'][i])
        #table_source['Dst_to_ctr_rad'][i] = 3.14*(table_source['Dst_to_ctr_deg'][i])/180
    
    return(table_source)


"""
Function whose output is the distance between 2 points in radians.

Libraries: import astropy
           from astropy.table import Table/ import pandas as pd
           from astropy.coordinates import SkyCoord
           import numpy as np

Inputs: table_source -> table object (astropy or dataframe) that contains all the data from the galaxies.
                        It needs 2 columns with their position in ra,dec (deg)
          
Outputs: dist_rad -> float? distance in rad 

To improve: ???
"""

def r_proj(table_source):
    
    #-----------------------------------------------------------------------------------
    #SkyCoord column
    coord = SkyCoord(table_source['ra'],table_source['dec'], unit='deg')
    
    #-----------------------------------------------------------------------------------
    #Angular distance to the centre
    hydra_ctr = SkyCoord.from_name('Hydra I Cluster')

    #-----------------------------------------------------------------------------------
    #table_source['SkyCoord'] = SkyCoord(table_source['RA'],table_source['DEC'])
    Skyc_arr = SkyCoord(table_source['ra'],table_source['dec'], unit='deg')
    dist = hydra_ctr.separation(Skyc_arr).deg
    dist_rad = np.radians(dist)
    
    return(dist_rad)