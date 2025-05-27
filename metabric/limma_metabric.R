#############
# Libraries #
#############

library(edgeR)
library(data.table)
library(tidyverse)
library(ggpubr)
library(Glimma)
library(biomaRt)
library(EnhancedVolcano)
library(ggfortify)
library(dplyr)
###################
# Expression file #
###################

setwd("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/Metabric")
#Upload METABRIC file 
expr <-fread("Ensembl_Metabric_17724_duplicates removed .csv")
expr <- expr[,-c(2,3)]

genes <- subset(expr, expr$Hugo_Symbol %in% "TRIM45")
rownames(genes) <- genes[,1]
genes <- genes[,-1]
genes <- as.data.frame(t(genes))
colnames(genes) <- "TRIM45"


expr <- column_to_rownames(expr, "Hugo_Symbol")

expr_1 <- as.data.frame(complete.cases(expr)) #7samples have NAs expression 
expr <- na.omit(expr) #remove the seven genes

#Read clinical file 
info <- read.csv2("brca_metabric_clinical_data.csv",sep = ";", as.is = T, check.names = F)
info_rec <- info[,-c(1,3:9,11)]
info_rec <- info_rec[,-c(5:30)]

info_rec$`Pam50 + Claudin-low subtype` <- gsub("claudin-low", "Basal", info_rec$`Pam50 + Claudin-low subtype`)

genes <- tibble::rownames_to_column(genes, "ID")
merged <- merge(genes, info_rec, by.x= "ID", by.y = "Patient ID")

#Determine t45 expression level and filter for only ER positive
merged <- subset(merged, merged$`ER status measured by IHC` == "Positve")
merged$TRIM45_expression <- ifelse(merged$TRIM45 >= median(merged$TRIM45), 'High', "Low")

#Take top and bottom 25% of t45 expression
top <- merged %>% filter(merged$TRIM45 > quantile(merged$TRIM45, 0.75))
bottom <- merged %>% filter(merged$TRIM45 < quantile(merged$TRIM45, 0.25))

combined <- rbind(top, bottom)
merged <- combined

#merged <- merged %>% mutate_if(is.character, ~na_if(.,"NC"))

#Make sure the ids are similiar and in the same order 
filtered_colnames <- intersect(names(expr), merged$ID)
df1_filtered <- expr[, filtered_colnames, drop = FALSE]
expr <- df1_filtered

rownames(merged) <- merged[,1]
merged$ID <- NULL
merged <- merged %>% arrange(TRIM45_expression)
expr<- expr[,rownames(merged)]
all(rownames(merged) == colnames(expr))


#######
# PCA #
#######
str(merged)
pca <- prcomp(t(expr), scale = TRUE)
pca_summary <- summary(pca)
pca_summary$importance[, 1:100]

merged <- rownames_to_column(merged, "id")
merged <- as.data.frame(merged)
expr <- as.data.frame(expr)

pca_df <- data.frame(pca$x[,1:2]) %>%
  rownames_to_column("id") %>%
  inner_join(select(merged, c(id, TRIM45_expression, `Pam50 + Claudin-low subtype`)),by = "id")


pca_plot <- ggplot(pca_df,aes(x = PC1,y = PC2, color = `Pam50 + Claudin-low subtype`)) +
  geom_point() + 
  theme_classic() 
pca_plot


ggsave("pca_subtype_metabric_top_bottom_25.pdf",plot = pca_plot)

#make group information 
group <- as.factor(merged$TRIM45_expression)
subtype <- as.factor(merged$`Pam50 + Claudin-low subtype`)

design <- model.matrix(~0+group+subtype)
design
colnames(design) <- c("High", "Low","HER2_Enriched", "Luminal_A", "Luminal_B","NC", "Normal_Like")

contr <- makeContrasts(High - Low, levels = design)
contr

#########
# Limma #
#########
fit <- lmFit(expr, design)
fit.contr <- contrasts.fit(fit, contr)
fit.contr <- eBayes(fit.contr)

volcanoplot(fit.contr)
top.table <- topTable(fit.contr, sort.by = "P", n = Inf, coef = 1, adjust.method = "BH")


write.csv2(top.table, "limma_metabric_er_positive_top_bottom_t45_level_BH.csv")

result <- decideTests(fit.contr)
summary(result)

top.table <- subset(top.table, rownames(top.table) != "TRIM45")

v <- EnhancedVolcano(top.table,
                     lab = rownames(top.table),
                     labSize = 6.0,
                     labCol = 'black',
                     labFace = 'bold',
                     title = "High TRIM45 vs Low TRIM45",
                     x = 'logFC',
                     y = 'adj.P.Val',
                     caption = paste(
                       "Total significant genes (Padj): ",
                       sum(top.table$adj.P.Val < 0.05), ", Up: ",
                       sum(top.table$adj.P.Val < 0.05 & top.table$logFC > 0), " Down: ",
                       sum(top.table$adj.P.Val < 0.05 & top.table$logFC < -0)),
                     colAlpha = 0.5, 
                     FCcutoff = 0.5)
v


pdf("2025_03_25_volcanoplot_limma_metabric_er_positive_top_vs_bottom_t45_level_BH.pdf", heigh = 10, width =10)
print(v)
dev.off()

png("limma/volcanoplot_limma_metabric_er_positive_top_bottom_t45_level_no_corr.png",res= 200, heigh = 2000, width =1800)
print(v)
dev.off()