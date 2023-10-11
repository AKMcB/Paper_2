library(dplyr)
library(forcats)
library(biomaRt)
library(ggplot2)
library(ggpubr)
library(tidyverse)

#Dataset
expr <- read.csv2("BCCL_filtered_log2.csv", sep=";", as.is = T, check.names = F)

genes<- subset(expr, expr$Name %in% c("TRIM45", "TRIM27", "TRIM32"))

genes<- subset(expr, expr$Name %in% c("ATG2A", "ATG2B", "ATG3", "ATG4A", "ATG4B"
                                      , "ATG4C", "ATG4D", "ATG5", "ATG7","ATG9A","ATG9B",
                                      "ATG10", "ATG12", "ATG13", "ATG14", "ATG16L1",
                                      "ATG16L2", "ATG101", "ULK1", "BECN1", "RB1CC1",
                                        "WIPI1", "WIPI2", "WDR45B", "WDR45", "MAP1LC3A",
                                      "MAP1LC3B", "MAP1LC3C","GABARAP", "GABARAPL1", 
                                      "GABARAPL2"))

rownames(genes) <- genes[,1]
genes <- genes[,-1]
t_genes<- t(genes)
t_genes<- as.data.frame(t_genes)
t_genes <- tibble::rownames_to_column(t_genes, "id")

write.csv2(t_genes, "t_autophagy proteins.csv")
t_genes<- read.csv2("t_autophagy proteins.csv", sep=";", as.is = T, check.names = F)

#Sample info 
ids <- read.csv("sample_info.csv", sep=";", as.is = T, check.names = F)
ids<- ids[ ,c(2,4,21)]
info_BC<- subset(ids, ids$CCLE_Name %in% t_genes$id)

#Merging
merged<- merge(t_genes, info_BC, by.x= "id", by.y= "CCLE_Name" )
merged<- merged[-c(7, 11, 32), ] #removed beacuse subtype dont fit

rownames(merged)<-merged[,1]
merged<- merged[,-1]
str(merged)

merged$lineage_molecular_subtype      <- factor(merged$lineage_molecular_subtype, 
                                              levels= c("HER2_amp", "basal_A", "basal_B", "luminal"), 
                                              labels = c("HER2", "Basal A", "Basal B", "Luminal"))


my_comparisons <- list( c("Luminal", "HER2"),
                        c("Luminal", "Basal A"),
                        c("Luminal", "Basal B")) 
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("****", "***", "**", "*", "ns"))

colnames(merged)[5] <- "PAM50"

p <- ggboxplot(merged, x="PAM50", y="TRIM45",add="jitter", add.params = list(alpha=0.7,size=2),
               color = "black", shape = 21,
               fill="PAM50", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#A0D636","#DF2DE0","#333ED4"), 
               order = c("Luminal","HER2", "Basal A", "Basal B"),
               ylab = "Expression", xlab = "PAM50", title = "TRIM45",
               ggtheme = theme_pubr(legend = "right")) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic"), #Italic if it is a gene. 
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))+
  geom_signif(comparisons = my_comparisons,map_signif_level = T, y_position = c(4.0, 4.3,4.6), textsize=4)


p 

pdf("TRIM45_exp_BCCL_PAM50.pdf", height = 5, width = 5)
print(p)
dev.off()


ggexport(p, filename = "TRIM45_exp_BCCL_PAM50.pdf", res = 200, height = 2000, width = 2000)

