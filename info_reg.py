"""
Function that reads a region object extracted using region on a ds9 region file.

Libraries:  import regions

Inputs: region -> region object
        col_arr -> 1x2 array with the colours assigned for regions with Halpha and Halpha+Hi emission

Outputs: names -> string list containing the names of the bodies contained in the region object
         emission -> SkyCoords list containing coordinates in ra and dec (deg) of bodies
         coords -> SkyCoords list containing coordinates in ra and dec (deg) of bodies

To improve -> look whether there is a common format for regions (not just ds9) and generalise 
"""
def info_reg(region):
    #Declare empty lists
    names = []
    emission = []
    coord = []
    
    #Iterate for every galaxy
    for p in range(len(region)):
        coord.insert(p, region[p].center)
        names.insert(p, region[p].meta['text'])
        emission.insert(p, region[p].visual['edgecolor'])
        
    return (names,emission,coord)