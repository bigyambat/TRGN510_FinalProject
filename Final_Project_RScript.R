
#Bioconductor Install for version 3.16
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install(version = "3.16")


BiocManager::install("vsm")

library(DESeq2)
library(dplyr)
library(apeglm)

setwd("/Users/bigyambat/Documents/GitHub/TRGN510_FinalProject")


#Data Wrangling

#Refer to comment below for reason behind setting header = FALSE. If issue does not apply, skip lines 15-17
Star_Matrix <- read.csv(file = 'CSV_Manifest3.csv', header = FALSE)

#Known Issue => R Studio messed up the header (with long file names). So, set header to false to avoid issue and make Row1 the column name
colnames(Star_Matrix) <- Star_Matrix[1,]
Star_Matrix <- Star_Matrix[-1,]


#Removing irrelavent unstranded reads (ie N_Unmapped, etc)
Star_Matrix2 <- Star_Matrix[-(1:4),]



#Reading in clinical.tsv file
clinical_sample <- read.csv(file = 'clinical.tsv', sep = "\t", header = TRUE)

#Clinical Data file has 2 rows for each sample. Selecting for every other row. 
clinical_sample2 <- clinical_sample[c(rep(FALSE,1),TRUE), ]

#Assigning row names starting from 1
rownames(clinical_sample2) <- 1:nrow(clinical_sample2)

#Selecting for relavent columns in clinical sample
clinical_sample3 <- clinical_sample2[, c(2,1,12)]

#Ordering clinical sample by gender
clinical_sample4 <- clinical_sample3[order(clinical_sample3$gender),]

#Making clinical sample ids as row names
rownames(clinical_sample4) <- clinical_sample4$case_submitter_id

#Removing irrelavent 1st column
clinical_sample4 <- clinical_sample4[-c(1)]




#####################################################################################



#Sample_ID & Case_ID correlation

#Reading Sample Name files which has all the different formats/names for each sample
Sample_Name_File <- read.csv(file = 'gdc_sample_sheet.2022-11-21.tsv', sep = "\t", header = TRUE)

#Selecing for File name & case.id columns
Sample_Name_File2 = subset(Sample_Name_File, select = -c(1, 3:5, 7:8) )

#Changing Case.ID Column to "case_submitter_id" to match Star_Matrix column
colnames(Sample_Name_File2) <- c("File.Name", "case_submitter_id")


##########################################################################################################################

#Changing the column names (file names) of Star_Matrix2 to their respective case_submitter_id counterpart
#Known Issue: Since case_submitter ids are out of order compared to the Star_Matrix file names, this has to be done manually by cross-referencing Sample_Name_File2 to Star_Matrix2 :(
#Easiest way to do this is to alphabetize  Sample_Name_File2 

colnames(Star_Matrix2) <- c("gene_id","TCGA-P7-A5NX","TCGA-SP-A6QC", "TCGA-S7-A7WL", "TCGA-RW-A7D0", "TCGA-WB-A817", "TCGA-W2-A7H7", "TCGA-WB-A81T", "TCGA-S7-A7WN", "TCGA-QR-A70O", "TCGA-XG-A823", "TCGA-WB-A81G", "TCGA-S7-A7WW", "TCGA-SR-A6MZ", "TCGA-QR-A707", "TCGA-QR-A70E", "TCGA-QR-A703", "TCGA-QR-A70V", "TCGA-QR-A70K", "TCGA-P8-A5KC", "TCGA-WB-A81D", "TCGA-WB-A81W", "TCGA-S7-A7X2", "TCGA-QR-A6ZZ", "TCGA-QR-A700", "TCGA-QR-A6GR", "TCGA-WB-A80L", "TCGA-QR-A70A", "TCGA-S7-A7WO", "TCGA-WB-A81R", "TCGA-P7-A5NY", "TCGA-RW-A68G", "TCGA-RW-A67W", "TCGA-WB-A81P", "TCGA-RW-A681", "TCGA-SR-A6N0", "TCGA-RW-A686", "TCGA-WB-A816", "TCGA-P7-A5NY", "TCGA-RW-A68F", "TCGA-RW-A68A", "TCGA-W2-A7H5", "TCGA-SR-A6MS", "TCGA-WB-A815", "TCGA-SR-A6MV", "TCGA-WB-A80Q", "TCGA-WB-A822", "TCGA-S7-A7WP", "TCGA-SP-A6QK", "TCGA-QR-A6GY", "TCGA-W2-A7HF", "TCGA-RW-A686", "TCGA-SR-A6MP", "TCGA-WB-A81F", "TCGA-WB-A81M", "TCGA-P8-A6RY", "TCGA-SR-A6MT", "TCGA-W2-A7HD", "TCGA-SQ-A6I4", "TCGA-WB-A80K", "TCGA-RW-A685", "TCGA-QR-A6GZ", "TCGA-WB-A80V", "TCGA-WB-A818", "TCGA-RW-A68D", "TCGA-S7-A7WM", "TCGA-RW-A68C", "TCGA-W2-A7HA", "TCGA-QT-A5XL", "TCGA-QR-A70M", "TCGA-WB-A81J", "TCGA-RW-A67X", "TCGA-QR-A6GO", "TCGA-QR-A70G", "TCGA-RW-A68B", "TCGA-QR-A70C", "TCGA-QR-A6GU", "TCGA-WB-A80M", "TCGA-WB-A80O", "TCGA-PR-A5PH", "TCGA-QR-A70W", "TCGA-TT-A6YJ", "TCGA-WB-A80N", "TCGA-QT-A5XP", "TCGA-S7-A7WX", "TCGA-SP-A6QI", "TCGA-RT-A6YC", "TCGA-SP-A6QD", "TCGA-QR-A70N", "TCGA-RW-A688", "TCGA-QR-A6GS", "TCGA-S7-A7WR", "TCGA-TT-A6YP", "TCGA-QT-A5XO", "TCGA-WB-A81K", "TCGA-WB-A81V", "TCGA-TT-A6YN", "TCGA-QT-A5XM", "TCGA-SP-A6QF", "TCGA-WB-A81H", "TCGA-SQ-A6I4","TCGA-SR-A6MY", "TCGA-S7-A7WV", "TCGA-QR-A6H4", "TCGA-QR-A6H2", "TCGA-W2-A7HE", "TCGA-QR-A6GX", "TCGA-SR-A6MU", "TCGA-P8-A5KC","TCGA-TT-A6YO", "TCGA-W2-A7HC", "TCGA-SP-A6QH", "TCGA-QT-A5XK", "TCGA-QR-A6H1", "TCGA-WB-A81S", "TCGA-QR-A70R", "TCGA-WB-A821", "TCGA-S7-A7X1", "TCGA-RW-A684", "TCGA-QR-A70H", "TCGA-QR-A6H5", "TCGA-QR-A70U", "TCGA-QR-A7IN", "TCGA-P8-A5KD", "TCGA-QR-A70J", "TCGA-WB-A81I", "TCGA-WB-A81N", "TCGA-P8-A5KD", "TCGA-RT-A6YA", "TCGA-QR-A708", "TCGA-QT-A69Q", "TCGA-S7-A7WT", "TCGA-PR-A5PF", "TCGA-QR-A6GW", "TCGA-SQ-A6I6", "TCGA-WB-A81A", "TCGA-RW-A67Y", "TCGA-WB-A819", "TCGA-QR-A70X", "TCGA-QR-A6GZ", "TCGA-RW-A67V", "TCGA-QT-A5XJ", "TCGA-W2-A7UY", "TCGA-SP-A6QJ", "TCGA-PR-A5PG", "TCGA-W2-A7HB", "TCGA-RT-A6Y9", "TCGA-S7-A7WU", "TCGA-S7-A7WQ", "TCGA-RW-A689", "TCGA-QT-A5XN")


##########################################################################################################################

#Extra Wrangling Steps Known Issue
#There are 150 samples in the Star_Matrix compared to 144 in the clinical_sample I cross-referenced the two to remove extra samples in Star_Matrix

#There are duplicated sample names in the Star_Matrix. I used the command below to remove duplicate column names. 
Star_Matrix3 <- Star_Matrix2[, !duplicated(colnames(Star_Matrix2))]

#Changing column 1 (gene_id) to row names and removing 1st column (which is now redundant)
Star_Matrix4 <- Star_Matrix3
row.names(Star_Matrix4) <- Star_Matrix3$gene_id
Star_Matrix4 <- Star_Matrix4[-c(1)]

#Ordering Star_Matrix to align with order of row names in clinical_sample
Star_Matrix5 <- Star_Matrix4[match(rownames(clinical_sample4), colnames(Star_Matrix4))]

#Converting Star_Matrix from datafram to a data matrix
Star_Matrix6 <- data.matrix(Star_Matrix5)


###############################################

#Add instructions on using Bioconductor and loading DESEQ2

#List of Packages that might be needed
library(BiocManager)
BiocManager::install("Named The Package You Need")


library(tibble)
library(tidyverse)
library(apeglm)
library(ggplot2)
library(vsn)
library(pheatmap)
library(ReportingTools)
library(DESeq2)

###########################################################################################################################

#Vignette start

#Load in the Star_Matrix6 & clinical_sample5 as cts and coldata respectively 

cts <- Star_Matrix6
coldata <- clinical_sample4



library(DESeq2)



#Loading all counts, clinical samples, and selecting for gender as variable of interest into dds object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ gender)

#DESeq command tests for the following:

#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 13396 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing

#Note: Command may take 5-10 minutes
dds <- DESeq(dds)




#Pre-Filtering

#Pre-Filtering is sometimes required to remove extraneous or unnecessary . Independent filtering is later on this dataset. 

#Number of total genes from dds object
nrow(dds)

#Filtering for genes that have at least 10 reads across each gene
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Number of total genes after pre-filtering
nrow(dds)

#It looks like all the genes have at least 10 reads. So, there was no change in the dataset. 

#Only need to run this you pre-filtered out genes.
dds <- DESeq(dds)


#Prints results of dds to a new variable. Only need to run this
res <- results( dds )
head(res)

#Meaning of Data
#Base Mean = Average of the normalized count values
#log2(FoldChange) = Change in gene expression between male and female
#lfcSE = Standard Error of the log2 fold change values
#stat = Wald's test to determine the weighted distance between gene expression
#pvalue = Hypothesis test to tell whether expression difference is significant
#padj Adjusted P values based on the Benjamini-Hochberg adjustment





#Log fold change shrinkage for visualization and ranking using lfcShrink function 

#Shrinkage does not change the fold changes of differential gene expression. Its function is to help with downstream assessments of results.
#LfcShrink adds a shrunken log2fold change and standard error results table from DESeq. Use the apeglm method for best shrinkage results. 

resultsNames(dds)

#Use the output of resultNames(dds) as the coef for the lfcShrink command (May take some time)
resLFC <- lfcShrink(dds, coef="gender_male_vs_female", type="apeglm")

resLFC

#Citation for using apeglm with lfcShrink command
#using 'apeglm' for LFC shrinkage. If used in published research, please cite:
  #Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
#sequence count data: removing the noise and preserving large differences.
#Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895




#P-Values and Adjusted P-Values

#Ordering P-Values by smallest
resOrdered <- res[order(res$pvalue),]

summary(res)


#P-Values less than 0.1
sum(res$padj < 0.1, na.rm=TRUE)


#Use results function below to generate results with alpha set to 0.05. 




#MA Plots

#MA plots show the log2(fold) differences with respect to the mean of normalized counts for all samples

#MA Plot of All Genes
plotMA(res, ylim=c(-2,2))

#MA Plot of Shrunken Genes
plotMA(resLFC, ylim=c(-2,2))


#Plot Counts

#Plots normalized counts with a pseudocount of 0.5. Variable of interest is specified as intgroup
plotCounts(dds, gene=which.min(res$padj), intgroup="gender")



#Extracting transformed values

#Two types of transformation methods: vst and rlog
#Transformed values are required to create PCA Plots and Heatmaps. Use blind = FALSE for faster runtime. 
#Note: Rlog may take a couple hours 



vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

head(assay(vsd), 3)

















 