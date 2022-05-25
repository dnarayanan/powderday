import numpy as np
import re
import pdb


class DrainePAH:
    def __init__(self,size_list,size,lam,lum,cabs):
        self.size_list=size_list
        self.size=size
        self.lam=lam
        self.lum=lum
        self.cabs=cabs


def get_sizes(data,start_rows):
    sizes = []

    for start_row in start_rows:
        sizes.append(float(data[start_row-9].split()[0]))
        
    return sizes


#--------------------------------

#Compute the location of the starting
#rows for each table of spectra, as well as the number of rows
#in each table

#--------------------------------

def read_draine_file(filename):
    PAH_list = []

    with open(filename, 'r') as fh:
        data = fh.readlines()
        #start_rows = [ii for ii,row in enumerate(data) if 'lambda   nu*P_nu   C_abs' in row]
        start_rows = [ii for ii,row in enumerate(data) if '0.0912' in row]

    #find the table length
    table_length = start_rows[1]-start_rows[0]-11 #the number of lines in between each row

    
    sizes = get_sizes(data,start_rows)

    #--------------------------------
    
    #Loop through the table and generate the values for the objects
    
    #--------------------------------

    
    for counter,start_row  in enumerate(start_rows):
        lam = []
        lum = []
        cabs = []
        
        table_for_this_grain_size = data[start_row:table_length+start_row]
        
        for i in range(table_length):
            lam.append(float(table_for_this_grain_size[i].split()[0]))
            
            lum_str = table_for_this_grain_size[i].split()[1]
            lum_str = lum_str.replace('D','e')
            lum.append(float(lum_str))
            
            cabs_str = table_for_this_grain_size[i].split()[2]
            cabs_str = cabs_str.replace('D','e')
            cabs.append(float(cabs_str))
            

        PAH_list.append(DrainePAH(sizes,sizes[counter],lam,lum,cabs))
        

    return PAH_list
