#############
# Libraries #
#############
library(tidyverse)
library(data.table)
library(gson)
library(matrixStats)
#########################
# Read in the edge file #
#########################
expr <- as.data.frame(fread("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/NCMM/TCGA-BRCA/lioness_output/lioness_filtered_for_genes_trim45.csv"))
head(expr)[1:5]

#Read in the gene list
list1 <- read.gmt("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/NCMM/TCGA-BRCA/gene_lists/HALLMARK_ESTROGEN_RESPONSE_EARLY.v2023.2.Hs.gmt")

list2 <- read.gmt("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/NCMM/TCGA-BRCA/gene_lists/HALLMARK_ESTROGEN_RESPONSE_LATE.v2023.2.Hs.gmt")

#Combine the gene files
genes <- rbind(list1, list2)
nrow(genes)

#Check for duplicates
genes <- genes[!duplicated(genes$gene),]
genes$term<- NULL
nrow(genes) #Removed 18 genes that were duplicated
head(genes)

expr <- subset(expr, expr$TF %in% genes$gene)

expr <- expr[,-2] #all targets are affect TRIM45
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.matrix(expr)

row_median <- as.data.frame(rowMedians(expr))
row_median

row_median <- rownames_to_column(row_median, "gene")
colnames(row_median)[2] <- "value"

row_median$type <- ifelse(row_median$value >=0, "Positive", "Negative")

p <-ggplot(row_median, aes(x = reorder(gene, value), y = value, fill = type), color= NA) +
  geom_bar(stat = "identity", )+
  scale_fill_manual(values = c("Positive" = "firebrick2", "Negative" = "steelblue2")) +  
  labs(title = "Top 10 and Bottom 10 TF affecting TRIM45",
       x = "Transcription Factor",
       y = "Mean value") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 18, color = "black"),
        axis.text.y = element_text( size = 18, color = "black"),
        panel.background = element_rect(fill = "white", colour = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))

p

png("edge_weight_TF_trim45_er_response_comb_median.png", width = 8, height = 6, units = "in", res = 300)  
plot(p)
dev.off()

pdf("edge_weight_TF_TRIM45_er_response_comb_median.pdf", width = 6, height = 4)  
plot(p)
dev.off()