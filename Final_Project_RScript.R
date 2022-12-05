
#BiocManager::install(c("DESeq2"))

library(DESeq2)
library(dplyr)

setwd("/Users/bigyambat/Documents/GitHub/TRGN510_FinalProject")


#Data Wrangling

#Refer to comment below for reason behind setting header = FALSE. If issue does not apply, skip lines 15-17
Star_Matrix <- read.csv(file = 'CSV_Manifest3.csv', header = FALSE)

#Known Issue => R Studio messed up the header (with long file names). So, set header to false to avoid issue and make Row1 the column name
colnames(Star_Matrix) <- Star_Matrix[1,]
Star_Matrix <- Star_Matrix[-1,]


#Removing irrelavent unstranded reads (ie N_Unmapped, etc)
Star_Matrix2 <- Star_Matrix[-(1:4),]


clinical_sample <- read.csv(file = 'clinical.tsv', sep = "\t", header = TRUE)

clinical_sample2 <- clinical_sample[c(rep(FALSE,1),TRUE), ]

rownames(clinical_sample2) <- 1:nrow(clinical_sample2)

#####################################################################################



#Sample_ID & Case_ID correlation

#Reading Sample Name files which has all the different formats/names for each sample
Sample_Name_File <- read.csv(file = 'gdc_sample_sheet.2022-11-21.tsv', sep = "\t", header = TRUE)

#Selecing for File name & case.id columns
Sample_Name_File2 = subset(Sample_Name_File, select = -c(1, 3:5, 7:8) )

#Changing Case.ID Column to "case_submitter_id" to match Star_Matrix column
colnames(Sample_Name_File2) <- c("File.Name", "case_submitter_id")

###############################################################################################################




#Changing the row names of Star_Matrix from tsv names to case_submiiter_id

Star_Matrix3(row.names) <- Star_Matrix2$gene_id

df.test <- data.frame(matrix(ncol = 1, nrow = 150))

Star_Matrix4 <- cbind(df.test, Star_Matrix3)

Sample_Name_File3 <- t(Sample_Name_File2)

##########################################################################################################################

#Changing the column names (file names) of Star_Matrix2 to their respective case_submitter_id counterpart
#Known Issue: Since case_submitter ids are out of order compared to the Star_Matrix file names, this has to be done manually by cross-referencing Sample_Name_File2 to Star_Matrix2 :(


colnames(Star_Matrix2) <- c("gene_id","TCGA-P7-A5NX","TCGA-SP-A6QC", "TCGA-S7-A7WL", "TCGA-WB-A817", "TCGA-W2-A7H7", "TCGA-WB-A81T", "TCGA-S7-A7WN", "TCGA-QR-A70O", "TCGA-XG-A823", "TCGA-WB-A81G", "TCGA-S7-A7WW", "TCGA-SR-A6MZ", "TCGA-QR-A707", "TCGA-QR-A70E", "TCGA-QR-A703", "TCGA-QR-A70V", "TCGA-QR-A70K", "TCGA-P8-A5KC", "TCGA-WB-A81D", "TCGA-WB-A81W", "TCGA-S7-A7X2", "TCGA-QR-A6ZZ", "TCGA-QR-A700", "TCGA-QR-A6GR", "TCGA-WB-A80L", "TCGA-QR-A70A", "TCGA-S7-A7WO", "TCGA-WB-A81R", "TCGA-P7-A5NY", "TCGA-RW-A68G", "TCGA-RW-A67W", "TCGA-WB-A81P", "TCGA-RW-A681", "TCGA-SR-A6N0", "TCGA-RW-A686", "TCGA-WB-A816", "TCGA-P7-A5NY", "TCGA-RW-A68F", "TCGA-RW-A68A", "TCGA-W2-A7H5", "TCGA-SR-A6MS", "TCGA-WB-A815", "TCGA-SR-A6MV", "TCGA-WB-A80Q", "TCGA-WB-A822", "TCGA-S7-A7WP", "TCGA-SP-A6QK", "TCGA-QR-A6GY", "TCGA-W2-A7HF", "TCGA-RW-A686", "TCGA-SR-A6MP", "TCGA-WB-A81F", "TCGA-WB-A81M", "TCGA-P8-A6RY", "TCGA-SR-A6MT", "TCGA-W2-A7HD", "TCGA-SQ-A6I4", "TCGA-WB-A80K", "TCGA-RW-A685", "TCGA-QR-A6GZ", "TCGA-WB-A80V", "TCGA-WB-A818", "TCGA-RW-A68D", "TCGA-S7-A7WM", "TCGA-RW-A68C", "TCGA-W2-A7HA", "TCGA-QT-A5XL", "TCGA-QR-A70M", "TCGA-WB-A81J", "TCGA-RW-A67X", "TCGA-QR-A6GO", "TCGA-QR-A70G", "TCGA-RW-A68B", "TCGA-QR-A70C", "TCGA-QR-A6GU", "TCGA-WB-A80M", "TCGA-WB-A80O", "TCGA-PR-A5PH", "TCGA-QR-A70W", "TCGA-TT-A6YJ", "TCGA-WB-A80N", "TCGA-QT-A5XP", "TCGA-S7-A7WX", "TCGA-SP-A6QI", "TCGA-RT-A6YC", "TCGA-SP-A6QD", "TCGA-QR-A70N", "TCGA-RW-A688", "TCGA-QR-A6GS", "TCGA-S7-A7WR", "TCGA-TT-A6YP", "TCGA-QT-A5XO", "TCGA-WB-A81K", "TCGA-WB-A81V", "TCGA-TT-A6YN", "TCGA-QT-A5XM", "TCGA-SP-A6QF", "TCGA-WB-A81H", "TCGA-SQ-A6I4","TCGA-SR-A6MY", "TCGA-S7-A7WV", "TCGA-QR-A6H4", "TCGA-QR-A6H2", "TCGA-W2-A7HE", "TCGA-QR-A6GX", "TCGA-SR-A6MU", "TCGA-P8-A5KC","TCGA-TT-A6YO", "TCGA-W2-A7HC", "TCGA-SP-A6QH", "TCGA-QT-A5XK", "TCGA-QR-A6H1", "TCGA-WB-A81S", "TCGA-QR-A70R", "TCGA-WB-A821", "TCGA-S7-A7X1", "TCGA-RW-A684", "TCGA-QR-A70H", "TCGA-QR-A6H5", "TCGA-QR-A70U", "TCGA-QR-A7IN", "TCGA-P8-A5KD", "TCGA-QR-A70J", "TCGA-WB-A81I", "TCGA-WB-A81N", "TCGA-P8-A5KD", "TCGA-RT-A6YA", "TCGA-QR-A708", "TCGA-QT-A69Q", "TCGA-S7-A7WT", "TCGA-PR-A5PF", "TCGA-QR-A6GW", "TCGA-SQ-A6I6", "TCGA-WB-A81A", "TCGA-RW-A67Y", "TCGA-WB-A819", "TCGA-QR-A70X", "TCGA-QR-A6GZ", "TCGA-RW-A67V", "TCGA-QT-A5XJ", "TCGA-W2-A7UY", "TCGA-SP-A6QJ", "TCGA-PR-A5PG", "TCGA-W2-A7HB", "TCGA-RT-A6Y9", "TCGA-S7-A7WU", "TCGA-S7-A7WQ", "TCGA-RW-A689", "TCGA-QT-A5XN")


##########################################################################################################################



#Vignette start
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ batch + condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_trt_vs_untrt")
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")