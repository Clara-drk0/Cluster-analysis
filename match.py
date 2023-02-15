"""
**********Tercera versión. Si la distancia entre las dos galaxias emparejadas es demasiado grande,
se asigna una galaxia vacia.
Function that extracts the unmatched galaxies from our catalogue 

Inputs: table_1 -> astropy table object that contains (at least) ra and dec coord in deg
        table_2 -> astropy table object that contains (at least) ra and dec coord in deg
        id1 -> string identifying the first table
        id2 -> string identifying the first table
        x1 -> string with the name of the table that will contain those
              elements (and their data) from table_1 that had a match in table_2
        _1 -> string with the name of the table that will contain those
              elements (and their data) from table_1 that DID NOT had a match in table_2
        x2 -> string with the name of the table that will contain those
              elements (and their data) from table_2 that had a match in table_1
        _2 -> string with the name of the table that will contain those
              elements (and their data) from table_2 that DID NOT had a match in table_1
        
Outputs: table_1,table_2 -> tables with new columns indicating the index in the other table of the element
                            matched to each row
         x1,_1,x2,_2 -> tables containing the data detailed above       
         
To improve: Check the matching. So far it matches the galaxies, but I haven't checked if the matching
            is well done. 
            ->I don't think it is completely accurate. It might have to do with the maximum distance considered for a max.
              For some galaxies it's too much and for others it's not enough.
              Maybe reformulate: find all galaxies below the maxdist and choose among them
              This last idea will imply checking individually those matched poorly.
            -> Also it should return only some kind of pointer between tables, not a hundred new tables
"""

def Match(table_1,table_2,id1,id2,x1,_1,x2,_2):  #table_1a = glx; table_wise = wise_cat_aux;  table_field = glx_unmatched
    
   #------------------------------------------------------------------------------------
    #Add a new column of booleans that indicates whether the galaxy has already been asigned or not.
    name1 ='x'+id2
    name2 ='x'+id1
    name1b = name1+'_ind'
    name2b = name2+'_ind'

    MatchMax = 0.01

    newcol = np.zeros(len(table_1))
    newcol = [int(x) for x in newcol]
    table_1[name1] = newcol
    table_1[name1b] = [None]*len(table_1)
    newcolr = np.zeros(len(table_2))
    newcolr = [int(x) for x in newcolr]
    table_2[name2] = newcolr
    table_2[name2b] = [None]*len(table_2)
    #-------------------------------------------------------------------------------------
    #Create empty array for distances between galaxies
    dist = np.zeros(len(table_1))
    c = SkyCoord(0, 0, unit = "deg")
    #Create skycoord objects for each element in both tables
    sky_1 = SkyCoord(table_1['ra'],table_1['dec'], unit = 'deg')
    sky_2 = SkyCoord(table_2['ra'],table_2['dec'],unit = 'deg')
    #print(len(dist),len(sky_1),len(sky_2))     #Reference for the dimensions in console

    #-------------------------------------------------------------------------------------
    #ALGORITHM
    #Step 1: Match according to coordinates: we estimate the distance of every element in one table to
    #        every element in the other, in case they have not been matched already (in which case 
    #        table_1[name1][j] = 1)
    #Step 2: If it has already been matched, the value assigned to the distance is our threshold, 
    #        so it can be discarded.
    #Step 3: Choose the minimum distance
    #Step 4: If the minimum distance is lower than the threshold, it's a match.
    #        Save the galaxy as matched in both tables and add the id of their couple.
    
    for i in range(len(table_2)):                               #Step 1
        for j in range(len(table_1)):
            if table_1[name1][j] != 1:                          #Step 2
                dist[j] = (sky_2[i].separation(sky_1[j])).deg
            else: 
                dist[j] = MatchMax 


        dmin = dist.min()                                       #Step 3
        posmin = np.where(dist == dmin)[0]

        if dmin < MatchMax:                                    #Step 4
            table_1[name1][int(posmin[0])] = 1
            table_1[name1b][int(posmin[0])] = id2+str(i)
            table_2[name2][i] = 1 
            table_2[name2b][i] = id1+str(int(posmin[0])) 

    #Check matching: both sums should be the same
    #suma1 = sum(table_1[name1])
    #suma2 = sum(table_2[name2])
    #print(suma1,suma2)

    #-------------------------------------------------------------------------------------
    #Create output tables from masks
    mask1x = table_1[name1] == 1
    mask1_ = table_1[name1] == 0
    x1 = table_1[mask1x]
    _1 = table_1[mask1_]

    mask2x = table_2[name2] == 1
    mask2_ = table_2[name2] == 0
    x2 = table_2[mask2x]
    _2 = table_2[mask2_]
   
   
    return(table_1,table_2,x1,_1,x2,_2)
