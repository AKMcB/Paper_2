#################
# Load libaries #
#################
library(ggpubr)
library(tidyverse)

########################
# Load expression data #
########################

setwd("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/Metabric")
expr <- read.csv2("Ensembl_Metabric_17724_duplicates removed .csv", sep=";", as.is = T, check.names = F)
expr <- expr[,-c(2,3)]

expr <- distinct(expr,Hugo_Symbol,.keep_all = T)
dup<-expr[duplicated(expr$Hugo_Symbol),] 

rm(dup)

#Subset the gene(s) of interest
genes <- subset(expr, expr$Hugo_Symbol %in% c("TRIM27", "TRIM32", "TRIM45"))

rownames(genes) <- genes[,1]
genes <- genes[,-1]

t_genes <- t(genes)
t_genes <- as.data.frame(t_genes)
t_genes$Meta_id <- rownames(t_genes)


######################
# Read clinical info #
######################

ids <- read.csv2("brca_metabric_clinical_data.csv", sep=";", as.is = T, check.names = F)
ids <- ids[,c(2,10)]

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
merged <- merge(t_genes, ids, by.x = "Meta_id", by.y = "Meta_id")
rownames(merged) <- merged[,1]
merged<- merged[, -1]

merged$`Pam50 + Claudin-low subtype` <- gsub("claudin-low", "Basal", merged$`Pam50 + Claudin-low subtype`)

table(merged$`Pam50 + Claudin-low subtype`)

colnames(merged)[4] <- "PAM50"

####################
# Make the boxplot #
####################

merged$PAM50   <- factor(merged$PAM50, 
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

merged <- merged[complete.cases(merged), ]

p <- ggplot(merged, aes(x = PAM50, y = TRIM45, fill = PAM50))+
  geom_point(alpha=0.5,position = position_jitter(width = 0.3, height = 0.5), shape= 21, size= 3)+
  geom_boxplot(fill = "white", alpha = 0.8, outlier.shape = NA) +
  labs(y = "TRIM45 (log2+1)", x = "PAM50", title = "METABRIC") + 
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
  geom_signif(comparisons = my_comparisons,map_signif_level = T, y_position = c(10.0, 10.3,10.6, 9.2,9.5,9.8,9.2), 
              textsize=10)

p 

pdf("TRIM45_exp_metabric_PAM50_3.pdf", height = 6, width = 6)
print(p)
dev.off()