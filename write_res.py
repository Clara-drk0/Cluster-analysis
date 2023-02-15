"""
Function that copies results to a file 

Libraries: from astropy.table import Table

Inputs: filename -> name of the file ('whatever.txt') where the data will be written. It needs to exist already.
        path -> location where the file will be saved.
        data -> 1x? array/Table with all the data to be saved.
        columns -> list with the names of the columns

Outputs: None. It writes in the file.

To improve -> It doesn't do that much, to be honest.
"""
def write_res(filename,path,data,columns):
    #Convert array to pandas dataframe to write in txt or csv
    df= pd.DataFrame(data,columns = columns)
    #In case data are in a table format, use this line
    #df = data.to_pandas()
    
    print(df)
    
    #String with complete path
    str_path = path + '\\' + filename
    
    #Open file
    file = open(str_path, "a+")
    
    #Write in file
    df.to_csv(str_path,sep = '\t',mode= 'a', index = False, header = False)
    
    #Close file
    file.close()