library(tibble)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidytext)
library(Hmisc)
#Start here when working with new genes----------------------------------------------------------------------
setwd("C:/Users/abe186/UiT Office 365/O365-Bioinformatikk TRIM27 - General/pancanAtlasHub")
merged_1 <- read.delim("EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena")

df <- subset(merged_1, merged_1$sample == "TRIM45")
rownames(df) <- df$sample
df$sample <- NULL
df <- as.data.frame(t(df))
df <- tibble::rownames_to_column(df, "id")
df$id <- gsub('\\.',"-",df$id)


ann <- readr::read_tsv("TCGA_phenotype_denseDataOnlyDownload.tsv")
colnames(ann)[1] <- "id"

merged <- merge(df, ann, by.x="id", by.y = "id")
merged <- merged[complete.cases(merged),]
merged$TRIM45 <- as.numeric(merged$TRIM45)
merged <- merged[order(merged$TRIM45, decreasing = TRUE),]
#head(merged)

colnames(merged)[5] <- "primary_disease"

fac <- with(merged, reorder(primary_disease, TRIM45, median, order = T))
merged$primary_disease <- factor(merged$primary_disease, levels = levels(fac))
#label(merged$primary_disease) <- "Primary Disease"

library(tools)
merged$primary_disease <- as.character(merged$primary_disease)
#change the first character on each word to a capital character
merged$primary_disease <- tools::toTitleCase(merged$primary_disease)
str(merged)



p <- ggboxplot(merged, x="primary_disease", y="TRIM45",outlier.shape = NA,
          fill ="primary_disease",
          palette = c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFD92F", "#E0F3F8", "#ABD9E9",
                      "#74ADD1", "#4575B4", "#8DA0CB", "#E78AC3", "#762A83", "#9970AB", "#C2A5CF", "#E7D4E8",
                      "#F7F7F7", "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#4DAC26", "#80CDC1", "#35978F",
                      "#7FC97F", "#008837", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
                      "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
                      "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A51A3", "#FFFF99", "#B15928", "#0571B0", "#CC4C02"), 
          ylab = "TRIM45 Expression", xlab = "Primary Disease") +
          guides(fill=guide_legend(title="Primary Disease"))+
          scale_x_reordered()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face ="italic"), #Italic if it is a gene. 
        axis.text.x = element_text(vjust = 0.5,angle = 90, size = 10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))+
          stat_compare_means(method = "anova", 
                     label.x.npc = "center",
                     label.y.npc = "top", vjust = 1.0, hjust = 0.5)

p

pdf("PanCan_TRIM45_2.pdf", height = 6, width = 8)
print(p)
dev.off()




