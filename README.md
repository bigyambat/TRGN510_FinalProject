# TRGN510 Final Project

## Title
Differential Expression Comparing Females vs Male Paragangliomas & Glomus Tumors from NCI GDC Datasets using DeSEQ2 package

## Author
Bigy Ambat 

## Description of Project
I will be comparing the differential expression of Paragangliomas & Glomus Tumors in Female vs Male patients. Analysis will using the DeSEQ2 R Package and will be using the following vignette: 
http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html. 

For this analysis, I will be data from the NCI GDC database from  Link to the dataset is below (with filters on)

https://portal.gdc.cancer.gov/exploration?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22content%22%3A%7B%22field%22%3A%22cases.demographic.vital_status%22%2C%22value%22%3A%5B%22alive%22%5D%7D%2C%22op%22%3A%22in%22%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.disease_type%22%2C%22value%22%3A%5B%22paragangliomas%20and%20glomus%20tumors%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.primary_site%22%2C%22value%22%3A%5B%22adrenal%20gland%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.samples.sample_type%22%2C%22value%22%3A%5B%22primary%20tumor%22%5D%7D%7D%5D%7D&searchTableTab=cases


The dataset includes 144 star-counts files (split into 61 Males and 83 Females). Datset is only using the TCGA cohort. Dataset is controlling for Primary Tumor and only considering those that are Alive. 

## Data
I will be using data from NCI GDC. As mentioned, there are 144 htc count files. However, the dataset below has not been filtered for the parameters above. Data wrangling will occur in R and be documented along with the vignette. 

Star-count files, clinical data, and other relavent files from GDC cart are uploaded to Github with this README

## Milestone 1
Due Date: Friday December 2nd
Data from each star-counts is formatted to fit MANIFEST.txt file. See Milestone 1 section below for full details once submitted. 

1) Star-count files from each folder were manually taken out of individual folders and placed into gdc_download_20221122_055006.438971 folder. 

2) Merging_Python.py script was created. This script took the unstranded columns from each star-counts.tsv file. I created a new txt file (test3.txt) that contained the ensemble genes column (gene_id column) from one of the star_counts tsv files. 

The script merged the unstranded columns to the test3.txt file

The test3.txt file was converted to CSV_Manifest3.csv to accomidate reading into RStudio. 


## Milestone 2
Due Date: Sunday December 4th 
Completion of data wrangling such that all dataframes are formatted in ideal manner. Dataframes can now be loaded into DESEQ2 vignette. 


## Deliverable 
Due Date: December 8th

An initial completion of the vignette including charts/data extracted as image files. Meet with Professor the following week to seek feedback and update

A complete repository with clear documentation an ddescription of analysis and results. Data presented in R' Markdown Script with block points indicating various section of project. 

**Deadline Extended with permission from Professor due to being sick**


