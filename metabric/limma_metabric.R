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

###################
# Expression file #
###################

setwd("C:/Users/abe186/UiT Office 365/O365-Bioinformatikk TRIM27 - General/Metabric")
#Upload METABRIC file 
expr <- read.csv2("Ensembl_Metabric_17724_duplicates removed .csv",sep = ";", as.is = T, check.names = F)
expr <- expr[,-c(2,3)]

genes <- subset(expr, expr$Hugo_Symbol %in% "TRIM45")
rownames(genes) <- genes[,1]
genes <- genes[,-1]
genes <- as.data.frame(t(genes))
colnames(genes) <- "TRIM45"

rownames(expr) <- expr[,1]
expr$Hugo_Symbol <- NULL

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
merged$TRIM45_expression <- ifelse(merged$TRIM45 >= median(merged$TRIM45), 'High', "Low")
merged <- subset(merged, merged$`ER Status` == "Positive")

#Take top and bottom 25% of t45 expression
top <- merged %>% filter(merged$TRIM45 > quantile(merged$TRIM45, 0.75))
bottom <- merged %>% filter(merged$TRIM45 < quantile(merged$TRIM45, 0.25))

combined <- rbind(top, bottom)
merged <- combined

merged <- merged %>% 
  mutate_if(is.character, ~na_if(.,"NC"))

merged <- merged[-357,]

#Make sure the ids are similiar and in the same order 
filtered_colnames <- intersect(names(expr), merged$ID)
df1_filtered <- expr[, filtered_colnames, drop = FALSE]
expr <- df1_filtered

rownames(merged) <- merged[,1]
merged$ID <- NULL
merged <- merged %>% arrange(TRIM45_expression)
expr<- expr[,rownames(merged)]
all(rownames(merged) == colnames(expr))

#make group information 
group <- as.factor(merged$TRIM45_expression)
subtype <- as.factor(merged$`Pam50 + Claudin-low subtype`)

design <- model.matrix(~0+group+subtype)
design
colnames(design) <- c("High", "Low","HER2_Enriched", "Luminal_A", "Luminal_B", "Normal_Like")

contr <- makeContrasts(High - Low, levels = design)
contr


#########
# Limma #
#########
fit <- lmFit(expr, design)
fit.contr <- contrasts.fit(fit, contr)
fit.contr <- eBayes(fit.contr)

volcanoplot(fit.contr)


top_genes <- topTable(fit.contr, number = Inf, adjust.method = "none")
print(top_genes)

write.csv2(top_genes, "limma_metabric_er_positive_top_bottom_t45_level_no_corr.csv")

result <- decideTests(fit.contr)
summary(result)

top_genes <- subset(top_genes, rownames(top_genes) != "TRIM45")

v <- EnhancedVolcano(top_genes,
                     lab = rownames(top_genes),
                     labSize = 6.0,
                     labCol = 'black',
                     labFace = 'bold',
                     title = "High TRIM45 vs Low TRIM45",
                     x = 'logFC',
                     y = 'adj.P.Val',
                     colAlpha = 0.5)
v


pdf("volcanoplot_limma_metabric_er_positive_top_bottom_t45_level_no_corr.pdf", heigh = 10, width =10)
print(v)
dev.off()

png("volcanoplot_limma_metabric_er_positive_top_bottom_t45_level_no_corr.png",res= 200, heigh = 2000, width =1800)
print(v)
dev.off()



pca <- prcomp(
  t(expr), # transpose our data frame to obtain PC scores for samples, not genes
  scale = TRUE # we want the data scaled to have unit variance for each gene
)
pca_summary <- summary(pca)
pca_summary$importance[, 1:100]

merged <- rownames_to_column(merged, "id")
meerged <- as.matrix(merged)
expr <- as.matrix(expr)

pca_df <- data.frame(pca$x[,1:2]) %>%
  # Turn samples IDs stored as row names into a column
  tibble::rownames_to_column("id") %>%
  # Bring only the variables that we want from the metadata into this data frame
  # here we are going to join by `refinebio_accession_code` values
  dplyr::inner_join(
    dplyr::select(merged, id, TRIM45_expression, `Pam50 + Claudin-low subtype`),
    by = "id"
  )


pca_plot <- ggplot(pca_df,aes(x = PC1,y = PC2, color = TRIM45_expression)) +
  geom_point() + # Plot individual points to make a scatterplot
  theme_classic() # Format as a classic-looking plot with no gridlines

# Print out the plot here
pca_plot


ggsave("pca_t45_level_metabric_top_bottom_25.pdf",plot = pca_plot)

