#############
# Libraries #
#############

library(tibble)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidytext)
library(Hmisc)
library(data.table)
library(tools)

########################
# Read expression file #
########################

setwd("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/pancanAtlasHub")
merged_1 <- as.data.frame(fread("EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"))

df <- subset(merged_1, merged_1$sample == "TRIM45")
rownames(df) <- df$sample
df$sample <- NULL
df <- as.data.frame(t(df))
df <- tibble::rownames_to_column(df, "id")
df$id <- gsub('\\.',"-",df$id)

ann <- as.data.frame(fread("TCGA_phenotype_denseDataOnlyDownload.tsv"))
#ann <- ann[1:3]
colnames(ann)[1] <- "id"

merged <- merge(df, ann, by= "id")
merged <- merged %>% filter(!grepl('-11', id))
merged <- merged[complete.cases(merged),]
merged$TRIM45 <- as.numeric(merged$TRIM45)
merged <- merged[order(merged$TRIM45, decreasing = T),]
#head(merged)

colnames(merged)[5] <- "cancer_type"

merged$cancer_type <- as.character(merged$cancer_type)

merged$cancer_type <- tools::toTitleCase(merged$cancer_type)

fac <- with(merged, reorder(cancer_type, TRIM45, median, decreasing = F))
merged$cancer_type <- factor(merged$cancer_type, levels = levels(fac))

p <- ggboxplot(merged, x="cancer_type", y="TRIM45",outlier.shape = NA,
          fill ="cancer_type",
          palette = c("#008837", "#1B9E77", "#D95F02", "#7570B3" ,"#4DAC26", "#66A61E", "#E6AB02","#7FC97F", 
                      "#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFD92F", "#E0F3F8", "#ABD9E9",
                      "#74ADD1", "#4575B4", "#8DA0CB", "#E78AC3", "#762A83", "#9970AB", "#C2A5CF", "#E7D4E8",
                      "#F7F7F7", "#D9F0D3", "#A6DBA0", "#5AAE61", "#E7298A","#1B7837", "#80CDC1", "#35978F",
                      "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
                      "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A51A3", "#FFFF99", "#B15928", "#0571B0", "#CC4C02"), 
          ylab = "TRIM45 Expression", xlab = "Cancer Code") +
          guides(fill=guide_legend(title="Primary Disease"))+
          scale_x_reordered()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face ="italic"), #Italic if it is a gene. 
        axis.text.x = element_text(size=10, angle = 90), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10))+
          stat_compare_means(method = "anova",
                     label.x.npc = "center",
                     label.y.npc = "top", vjust = 1.0, hjust = 0.5)

p

pdf("PanCan_TRIM45_5.pdf", height = 10, width = 15)
print(p)
dev.off()
