##################
# Load libraries #
##################

library(ggpubr)
library(tidyverse)

########################
# Load expression data #
########################
setwd("C:/Users/abe186/UiT Office 365/O365-PhD Anne - General/TCGA-BRCA/Raw data")
expr <- read.csv2("gdc_brca_expr_edit.csv", as.is = TRUE, check.names = FALSE)
expr <- expr[,-1]
expr <- distinct(expr,gene,.keep_all = T)
dup<-expr[duplicated(expr$gene),] #check for duplicates
rm(dup)

#Subset the gene(s) of interest

genes1 <- subset(expr, expr$gene == "TRIM45")

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

#merged <- merge(t_genes_interactors, t_genes1, by.x = "id", by.y = "id")
#write.csv(merged, "TRIM_interactors_autophagy.csv")

##Import the info file
ids <- read.csv("TCGA_BRCA_Updated_Clinical_Data.csv", sep=";", as.is = T, check.names = F)
ids <- ids[,1:3]
#Remove first column (if there is any without a column name) and keep distinct patients
##There should not be any duplicate in patient names if you are only taking the tumor samples
ids <- ids[,-2]
ids <- distinct(ids, TCGA_id, .keep_all = T)
?as.numeric

#Merge gene expression with tumor subtype
#Merge-Merge two data frames by common columns or row names, or do other versions of database join operations.
#by, by.x, by.y	- specifications of the columns used for merging
merged <- merge(t_genes1, ids, by.x = "id", by.y = "TCGA_id")
rownames(merged) <- merged[,1]
merged<- merged[, -1]

table(merged$BRCA_Subtype_PAM50)

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

colnames(merged)[2] <- "PAM50"


p <- ggplot(merged, aes(x = PAM50, y = TRIM45, fill = PAM50))+
  geom_point(alpha=0.5,position = position_jitter(width = 0.3, height = 0.5), shape= 21, size= 3)+
  geom_boxplot(fill = "white", alpha = 0.8, outlier.shape = NA) +
  labs(y = "TRIM45 (log2+1)", x = "PAM50", title = "TCGA-BRCA") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face ="bold"), 
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11, face ="italic"),#Italic if it is a gene. 
        axis.text.y = element_text(size = 10), 
        panel.background = element_rect(fill = "white",
                                        colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#A0D636","#DF2DE0","#333ED4"))+
  geom_signif(comparisons = my_comparisons,map_signif_level = T, y_position = c(7.0, 7.4,7.8, 5.9,6.3,6.7,5.9), textsize=10)
  
p 

pdf("TRIM45_exp_TCGA_BRCA_PAM50_3.pdf", height = 6, width = 6)
print(p)
dev.off()
