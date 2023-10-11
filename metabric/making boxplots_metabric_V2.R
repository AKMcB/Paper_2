library(dplyr)
library(forcats)
library(biomaRt)
library(ggplot2)
library(ggpubr)
library(tidyverse)

setwd("C:/Users/abe186/UiT Office 365/O365-Bioinformatikk TRIM27 - General/Metabric")
expr <- read.csv2("Ensembl_Metabric_17724_duplicates removed .csv", sep=";", as.is = T, check.names = F)
expr <- expr[,-c(2,3)]

expr <- distinct(expr,Hugo_Symbol,.keep_all = T)
dup<-expr[duplicated(expr$Hugo_Symbol),] #add comma-when Undefined columns selected

rm(dup)

#Subset the gene(s) of interest

genes <- subset(expr, expr$Hugo_Symbol %in% c("TRIM27", "TRIM32", "TRIM45"))

genes <- subset(expr, expr$Hugo_Symbol %in% c("TAX1BP1", "CALCOCO2", "OPTN", "SQSTM1", "NBR1"))

genes <- subset(expr, expr$Hugo_Symbol %in% c("WDR41", "LAMP2", "KIAA1033", "WASHC1", "TOMM22", 
                                             "STAM")) #WASH1 does not exsist in dataset


genes<- subset(expr1, expr1$Hugo_Symbol %in% c("ATG2A", "ATG2B", "ATG3", "ATG4A", "ATG4B"
                                      , "ATG4C", "ATG4D", "ATG5", "ATG7","ATG9A",
                                      "ATG10", "ATG14", "ATG16L1",
                                      "ATG16L2", "ULK1", "BECN1", "RB1CC1",
                                      "WIPI1", "MAP1LC3A","MAP1LC3B", "MAP1LC3C","GABARAP", 
                                      "GABARAPL1","GABARAPL2", "TBK1", "PIK3C3", "PIK3R4",
                                      "SNX30", "SNX4", "TAX1BP1", "CALCOCO2", "OPTN", "SQSTM1", "NBR1"))

#Remove the default rownames and add the gene name as row name of the data
rownames(genes) <- genes[,1]
genes <- genes[,-1]

#First we need to transpose the TRIM27 file so that we can merge it with the patient information
#t-transpose data t(x)
t_genes <- t(genes)


#Convert it into a dataframe and add a column with the patient ids
t_genes <- as.data.frame(t_genes)
t_genes$Meta_id <- rownames(t_genes)

#write.csv(t_genes, "autophagy receptors_Metabric.csv")

##Import the info file
ids <- read.csv2("brca_metabric_clinical_data.csv", sep=";", as.is = T, check.names = F)
ids <- ids[,c(2,10)]

#Remove first column (if there is any without a column name) and keep distinct patients
##There should not be any duplicate in patient names if you are only taking the tumor samples
rownames(ids) <- ids[,1]
t_ids <- t(ids)
t_ids <- as.data.frame(t_ids)
t_ids <- t_ids[-1,]
ids<- t(t_ids)
ids <- as.data.frame(ids)
ids$Meta_id <- rownames(ids)
ids <- distinct(ids, Meta_id, .keep_all = T)
rm(t_ids)

#Merge gene expression with tumor subtype
#Merge-Merge two data frames by common columns or row names, or do other versions of database join operations.
#by, by.x, by.y	- specifications of the columns used for merging
merged <- merge(t_genes, ids, by.x = "Meta_id", by.y = "Meta_id")
rownames(merged) <- merged[,1]
merged<- merged[, -1]

merged$`Pam50 + Claudin-low subtype` <- gsub("claudin-low", "Basal", merged$`Pam50 + Claudin-low subtype`)


colnames(merged)[4] <- "PAM50"


#Made a simpler title for subtype coulumn 
#converted claudin-low to Basal
#removed NC not classified
merged <- read.csv2("info_Metabric.csv", sep=";", as.is = T, check.names = F)
#Make the boxplot

merged$PAM50    <- factor(merged$PAM50, 
                                       levels= c("LumA", "LumB", "Her2", "Basal", "Normal"), 
                                       labels = c("Luminal A", "Luminal B", "HER2", "Basal", "Normal"))
my_comparisons <- list( c("Luminal A", "HER2"),
                        c("Luminal A", "Basal"),
                        c("Luminal A", "Normal"),
                        c("Luminal B", "HER2"),
                        c("Luminal B", "Basal"),
                        c("Luminal B", "Normal"),
                        c("Luminal A", "Luminal B")) 
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("****", "***", "**", "*", "ns"))

# Assuming 'df' is your dataframe

#merged$`PAM50 Subtype` <- ifelse(merged$`PAM50 Subtype` == "NC", NA, merged$`PAM50 Subtype`)
merged <- merged[complete.cases(merged), ]

p <- ggboxplot(merged, x="PAM50", y="TRIM45",add="jitter", add.params = list(alpha=0.7,size=2),
               color = "black", shape = 21,
               fill="PAM50", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#A0D636","#DF2DE0","#333ED4"), 
               order = c("Luminal A", "Luminal B","HER2", "Basal", "Normal"),
               ylab = "Expression", xlab = "PAM50", title = "TRIM45",
               ggtheme = theme_pubr(legend = "right")) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic"), #Italic if it is a gene. 
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))+
  geom_signif(comparisons = my_comparisons,map_signif_level = T, y_position = c(10.0, 10.3,10.6, 9.2,9.5,9.8,9.2), textsize=4)
  

p 

pdf("TRIM45_exp_metabric_PAM50.pdf", height = 6, width = 6)
print(p)
dev.off()
ggexport(p, filename = "TRIM45_exp_PAM50_subtype_Metabric.png",res = 200, height = 2000, width = 2000)
