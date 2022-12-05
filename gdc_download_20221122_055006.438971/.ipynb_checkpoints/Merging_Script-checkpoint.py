#!/usr/bin/python
import sys
import re
import fileinput
import pandas as pd 
import os 

#The function of this python script is to select the unstranded DNA column of each star-counts tsv file and bring that data into a new file. 

# directory name
dirname = '/Users/bigyambat/Documents/GitHub/TRGN510_FinalProject/gdc_download_20221122_055006.438971'
 
# extensions
ext = ('.tsv')
 
# scanning the directory to get required files

update_file = pd.read_csv("test3.txt", sep="\t")


#Searches through each tsv file and finds "unstranded" column. Appends it to the MANIFEST.txt (which is now converted to a dataframe). New file is CSV_Manifest

for files in os.scandir(dirname):
    if files.path.endswith(ext):
        #print(files)  # printing file name
        gene_file = pd.read_csv(files, sep="\t", header=1)
        unstranded_column = gene_file['unstranded']
        update_file[files.name] = unstranded_column

update_file.to_csv(r'/Users/bigyambat/Documents/GitHub/TRGN510_FinalProject/CSV_Manifest3.csv', index=False)






