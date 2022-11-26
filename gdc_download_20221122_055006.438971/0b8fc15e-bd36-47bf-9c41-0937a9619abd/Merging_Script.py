#!/usr/bin/python
import sys
import re
import fileinput
import pandas as pd 
import matplotlib.pyplot as plt

#The function of this python script is to select the unstranded DNA column of each star-counts tsv file and bring that data into a new file. 

# Usage of this script: 

gene_file = pd.read_csv(sys.argv[1], sep="\t", header=0)



unstranded_column = gene_file.columns[3]

unstranded_column_data = gene_file.loc[:unstranded_column]

print(unstranded_column_data)
