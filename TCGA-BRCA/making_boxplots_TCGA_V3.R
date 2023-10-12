library(dplyr)
library(forcats)
library(biomaRt)
library(ggplot2)
library(ggpubr)
library(tidyverse)

setwd("C:/Users/abe186/UiT Office 365/O365-Bioinformatikk TRIM27 - General/TCGA-BRCA/Raw data")
expr <- read.csv2("gdc_brca_expr_edit.csv", sep=";", as.is = T, check.names = F)
expr <- expr[,-1]
expr <- distinct(expr,gene,.keep_all = T)
dup<-expr[duplicated(expr$gene),] #add comma-when Undefined columns selected
rm(dup)

#Subset the gene(s) of interest

genes1 <- subset(expr, expr$gene %in% c( "TRIM27", "TRIM32", "TRIM45"))

genes_interactors <- subset(expr, expr$gene %in% c( "WDR41", "LAMP2", "KIAA1033", "WASH1", "TOMM22", 
                                             "STAM"))


genes<- subset(expr, expr$gene %in% c("ATG2A", "ATG2B", "ATG3", "ATG4A", "ATG4B"
                                      , "ATG4C", "ATG4D", "ATG5", "ATG7","ATG9A","ATG9B",
                                      "ATG10", "ATG12", "ATG13", "ATG14", "ATG16L1",
                                      "ATG16L2", "ULK1", "BECN1", "RB1CC1",
                                      "WIPI1", "MAP1LC3A","MAP1LC3B", "MAP1LC3C","GABARAP", 
                                      "GABARAPL1","GABARAPL2", "ULK2", "TBK1", "PIK3C3", "PIK3R4",
                                      "SNX30", "SNX4", "TAX1BP1", "CALCOCO2", "OPTN", "SQSTM1", "NBR1", 
                                      "TRIM45", "TRIM32", "TRIM27"))


#Remove the default rownames and add the gene name as row name of the data
rownames(genes1) <- genes1[,1]
genes1 <- genes1[,-1]

#First we need to transpose the TRIM27 file so that we can merge it with the patient information
#t-transpose data t(x)
t_genes1 <- t(genes1)


#Convert it into a dataframe and add a column with the patient ids
t_genes1 <- as.data.frame(t_genes1)
#t_genes$id <- rownames(t_genes)
t_genes1<- tibble::rownames_to_column(t_genes1, "id")

merged <- merge(t_genes_interactors, t_genes1, by.x = "id", by.y = "id")


#write.csv(merged, "TRIM_interactors_autophagy.csv")

##Import the info file
ids <- read.csv("TCGA_BRCA_Updated_Clinical_Data.csv", sep=";", as.is = T, check.names = F)
ids <- ids[,1:3]
#Remove first column (if there is any without a column name) and keep distinct patients
##There should not be any duplicate in patient names if you are only taking the tumor samples
ids <- ids[,-2]
ids <- distinct(ids, TCGA_id, .keep_all = T)


#Merge gene expression with tumor subtype
#Merge-Merge two data frames by common columns or row names, or do other versions of database join operations.
#by, by.x, by.y	- specifications of the columns used for merging
merged <- merge(t_genes1, ids, by.x = "id", by.y = "TCGA_id")
rownames(merged) <- merged[,1]
merged<- merged[, -1]

write.csv(merged, "TRIM27_autophagy_receptors_TCGA.csv")

#Make the boxplot
merged$BRCA_Subtype_PAM50    <- factor(merged$BRCA_Subtype_PAM50, 
                                                levels= c("LumA", "LumB", "Her2", "Basal", "Normal"), 
                                                labels = c("Luminal A", "Luminal B", "HER2-enriched", "Basal-like", "Normal-like"))

my_comparisons <- list( c("Luminal A", "HER2-enriched"),
                        c("Luminal A", "Basal-like"),
                        c("Luminal A", "Normal-like"),
                        c("Luminal B", "HER2-enriched"),
                        c("Luminal B", "Basal-like"),
                        c("Luminal B", "Normal-like"),
                        c("Luminal A", "Luminal B")) 
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("****", "***", "**", "*", "ns"))

colnames(merged)[4] <- "PAM50"

p <- ggboxplot(merged, x="PAM50", y="TRIM45",outlier.shape = NA,
               palette = c("#00AFBB", "#E7B800", "#FC4E07", "#A0D636","#DF2DE0","#333ED4"), 
               order = c("Luminal A", "Luminal B","HER2-enriched", "Basal-like", "Normal-like"),
               ylab = "TRIM45 Expression", xlab = "PAM50", title = "TCGA-BRCA",
               ggtheme = theme_pubr(legend = "right")) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face ="italic"), #Italic if it is a gene. 
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))+
  geom_signif(comparisons = my_comparisons,map_signif_level = T, y_position = c(7.0, 7.4,7.8, 5.9,6.3,6.7,5.9), textsize=10)+ 
  geom_jitter(shape = 21,stroke=0.1,size=3, aes(fill= PAM50, alpha=0.5), position = position_jitter(width = 0.3, height = 0.5))


p 

pdf("TRIM45_exp_TCGA_BRCA_PAM50_2.pdf", height = 6, width = 6)
print(p)
dev.off()


ggexport(p, filename = "TRIM45_exp_PAM50_subtype_BRCA.pdf", height = 15, width = 15)




