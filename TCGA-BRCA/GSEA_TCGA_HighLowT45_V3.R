library(dplyr)
library(tibble)
install.packages("caroline")
library(caroline)

#-------------------------------------------------------------------------------------
#Subsetting the TCGA dataset based on TRIM45 expression (High/Low)

expr <- read.csv2("gdc_brca_expr_edit.csv", sep = ";", as.is = T, check.names = F)
expr <- expr[,-2]

rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.data.frame(t(expr))

#Load info file 
info <- read.csv2("TCGA_BRCA_Updated_Clinical_Data.csv", sep = ";", as.is = T, check.names = F)
info <- info[,-c(2:44,46)]

#Merge expr and info based on patient ID
expr<- tibble::rownames_to_column(expr, "ID")
merged <- merge(expr, info, by.x = "ID", by.y = "TCGA_id" )

#Subset based on receptor status 
merged_ER <- subset(merged, merged$er_status_by_ihc == "Positive")
merged_HER <- subset(merged_ER, merged_ER$her2_status_by_ihc == "Negative")

#Define high low expression of TRIM45 (ENSG00000134253)
merged_ER$TRIM45_expression <- ifelse(merged_ER$ENSG00000134253 >= median(merged_ER$ENSG00000134253),
                                  'High', "Low")
merged_HER$TRIM45_expression <- ifelse(merged_HER$ENSG00000134253 >= median(merged_HER$ENSG00000134253),
                                      'High', "Low")

#Subset the df based on high and low TRIM45 expression 
HT45_ER <- subset(merged_ER, merged_ER$TRIM45_expression == "High")
LT45_ER <- subset(merged_ER, merged_ER$TRIM45_expression == "Low")

HT45_HER <- subset(merged_HER, merged_HER$TRIM45_expression == "High")
LT45_HER <- subset(merged_HER, merged_HER$TRIM45_expression == "Low")

#Transpose the df again and make into df again 
#Save the files 
rownames(HT45_ER) <- HT45_ER[,1]
HT45_ER <- HT45_ER[,-1]

rownames(LT45_ER) <- LT45_ER[,1]
LT45_ER <- LT45_ER[,-1]

rownames(HT45_HER) <- HT45_HER[,1]
HT45_HER <- HT45_HER[,-1]

rownames(LT45_HER) <- LT45_HER[,1]
LT45_HER <- LT45_HER[,-1]

HT45_ER <- as.data.frame(t(HT45_ER))
LT45_ER <- as.data.frame(t(LT45_ER))

HT45_HER <- as.data.frame(t(HT45_HER))
LT45_HER <- as.data.frame(t(LT45_HER))

#----------------------------------------------------------------------------------------------
#Merging the expression files for GSEA

#The patient id must include either High or low 
#High/low included so to compare these two groups
#This places High/Low before patient id 


colnames(HT45_ER) <- paste("High",colnames(HT45_ER), sep= "_")

colnames(LT45_ER) <- paste("Low",colnames(LT45_ER), sep= "_")

colnames(HT45_HER) <- paste("High",colnames(HT45_HER), sep= "_")

colnames(LT45_HER) <- paste("Low",colnames(LT45_HER), sep= "_")

library(tibble)

#Merge the files 
ER <- cbind(HT45_ER, LT45_ER)
HER <- cbind(HT45_HER, LT45_HER)

#Create a phenotype file-include only high/Low  
pheno_er <- ER[60486, ]
pheno_her <- HER[60486, ]

#Remove phenotype description (16090)
#Remove the NA genes by complete.cases function
#Add a description column
merged <- merged[-60484,]
HER_1 <- HER[complete.cases(HER),]

rm(ER_1,HER_1)
#merged <- add_column(merged, description = NA, .before = 2)

#save expression file 
#Edit the files in NotePad++ to fit the requirements for GSEA
write.csv2(HER, "HighLow_TRIM45_ER_pos_HER_neg_TCGA.csv")

#Save phenotype file 
write.csv2(pheno_her, "HighLow_TRIM45_ER_pos_HER_neg_pheno.csv")



















