# TRGN510 Final Project

## Title
Differential Expression Comparing Low Grade vs High Grade Glioma/Astrocytoma from cBioPortal using DeSEQ2 package

## Author
Bigy Ambat 

## Description of Project
I will be comparing the differential expression of Low Grade & High Grade Glioma/Asctrocytom in Pediatriac Brain Cancer patients. Analysis will using the DeSEQ2 R Package and will be using the following vignette: 
http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html. 

For this analysis, I will be using the cBioPortal from Pediatric Brain Cancer (CPTAC/CHOP, Cell 2020). Link to the dataset is below. 
https://www.cbioportal.org/study/summary?id=brain_cptac_2020

The dataset includes 76 htc-count files (split into 66 Pediatric Low Grade Gliomas and 10 Pediatric High Grade Gliomas). Dataset is examinging pediatric patients in 3  groups: 0-5, 5-10, and 10-15 age groups. I will be controlling for the primary tumor type. 

## Data
I will be using data from cBioPortal. As mentioned, there are 76 htc count files. However, the dataset below has not been filtered for the parameters above. Data wrangling will occur in R and be documented along with the vignette. 

https://cbioportal-datahub.s3.amazonaws.com/brain_cptac_2020.tar.gz

## Milestone 1
Due Date: Friday November 18th
Data is fully loaded into vignette through the HT-Seq set. 

## Milestone 2
Due Date: Friday November 25th 
An initial completion of the vignette including charts/data extracted as image files. Meet with Professor the following week to seek feedback and update

## Deliverable 
Due Date: December 3rd

A complete repository with clear documentation an ddescription of analysis and results. Data presented in R' Script with block points indicating various section of project. 


