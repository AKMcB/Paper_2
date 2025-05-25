##################
# Load libraries #
##################

library(ggpubr)
library(tidyverse)

########################
# Load expression data #
########################

setwd("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/Cancer cell lines cyclopedia/Data files")
expr <- read.csv2("BCCL_filtered_log2.csv", sep=";", as.is = T, check.names = F)

genes<- subset(expr, expr$Name %in% c("TRIM45", "TRIM27", "TRIM32"))

rownames(genes) <- genes[,1]
genes <- genes[,-1]
t_genes<- t(genes)
t_genes<- as.data.frame(t_genes)
t_genes <- tibble::rownames_to_column(t_genes, "id")

####################
# Read sample info #
####################

ids <- read.csv("sample_info.csv", sep=";", as.is = T, check.names = F)
ids<- ids[ ,c(2,4,21)]
info_BC<- subset(ids, ids$CCLE_Name %in% t_genes$id)

#Merging
merged<- merge(t_genes, info_BC, by.x= "id", by.y= "CCLE_Name" )
merged<- merged[-c(7, 11, 32), ] #removed beacuse subtype dont fit

rownames(merged)<-merged[,1]
merged<- merged[,-1]
str(merged)

merged$lineage_molecular_subtype  <- factor(merged$lineage_molecular_subtype, 
                                              levels= c("luminal","HER2_amp", "basal_A", "basal_B"), 
                                              labels = c("Luminal","HER2-enriched", "Basal-like A", "Basal-like B" ))


my_comparisons <- list( c("Luminal", "HER2-enriched"),
                        c("Luminal", "Basal-like A"),
                        c("Luminal", "Basal-like B")) 
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("****", "***", "**", "*", "ns"))

colnames(merged)[5] <- "PAM50"

p <- ggplot(merged, aes(x = PAM50, y = TRIM45, fill = PAM50))+
  geom_point(alpha=0.5,position = position_jitter(width = 0.3, height = 0.5), shape= 21, size= 3)+
  geom_boxplot(fill = "white", alpha = 0.8, outlier.shape = NA) +
  labs(y = "TRIM45 (log2+1)", x = "PAM50", title = "BCCL") + 
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
  geom_signif(comparisons = my_comparisons,map_signif_level = T, y_position = c(4.3, 4.6,4.9), 
              textsize=10)

p 

pdf("TRIM45_exp_BCCL_PAM50_3.pdf", height = 6, width = 6)
print(p)
dev.off()
