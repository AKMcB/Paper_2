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
expr <- as.data.frame(fread("Raw data/TCGA-BRCA.htseq_counts.tsv"))
expr$Ensembl_ID <- gsub("\\..*","", expr$Ensembl_ID)
#explore the expr file 
#test <- as.data.frame(colnames(expr))
#1097 primary cancer samples 
#7 metastatic 
#113 normal samples 

rownames(expr) <- expr[,1]
expr <- expr[,-1]

#There are duplicate samples but with different letter A, B or C 
#Remove the duplicated 01B and 01C samples 
cols_to_remove <- c("TCGA-A7-A0DB-01C", "TCGA-A7-A13E-01B",
                    "TCGA-A7-A13D-01B", "TCGA-A7-A0DC-01B",
                    "TCGA-A7-A26J-01B", "TCGA-A7-A26E-01B") 
#1091 primary cancer samples left


expr <- expr[, !(names(expr) %in% cols_to_remove), drop = FALSE]

#Have to remove the letter to late match with clin file 
#did not edit 11A, 11B and 06 since I do not want to include them 
colnames(expr) <- sub("-01A.*", "", colnames(expr))
colnames(expr) <- sub("-01B.*", "", colnames(expr))

#explore the expr file 
#test_edit <- as.data.frame(colnames(expr))
#After edited 1091 primary cancer saamples with -01 
# 7 metastatic samples with -06
#did not edit 11A and 11B since I do not want to include them 

rm(cols_to_remove)


#Xena browser uploads the data after log2(x+1) trasformation. 
#So converting it to linear scale to get the raw count data
expr <- 2^expr


expr <- subset(expr, !rownames(expr) %in% c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique"))
############################
# Create group information # 
############################

t45_1 <- read.csv2("limma/clinical_t45_ESR1_tcga_brca_receptor_status.csv", header = TRUE)
t45_1 <- t45_1[,-1]

t45_1$TRIM45_expression <- ifelse(t45_1$TRIM45 >= median(t45_1$TRIM45), 'High', "Low")
#Take the top 25% highest expression vs bottom 25% 
top <- t45_1 %>% filter(t45_1$TRIM45 > quantile(t45_1$TRIM45, 0.75))
bottom <- t45_1 %>% filter(t45_1$TRIM45 < quantile(t45_1$TRIM45, 0.25))

combined <- rbind(top, bottom)
t45_1 <- combined
#get same patients in expr as the clinical file
filtered_colnames <- intersect(names(expr), t45_1$id)
df1_filtered <- expr[, filtered_colnames, drop = FALSE]
expr <- df1_filtered


#Make sure the group info is in the same order as the expr file
t45_1 <- column_to_rownames(t45_1, "id")
t45_1 <- t45_1 %>% arrange(TRIM45_expression)
expr<- expr[,rownames(t45_1)]
all(rownames(t45_1) == colnames(expr))


group <- as.factor(t45_1$TRIM45_expression)
subtype <- as.factor(t45_1$BRCA_Subtype_PAM50)

####################
# gene annotations #
####################
raw_data <- expr
raw_data <- rownames_to_column(raw_data, "gene")
ensembl <- raw_data$gene


mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

G_list <- getBM(mart= mart, filters= "ensembl_gene_id", values= ensembl, attributes= c("ensembl_gene_id",
                                                                                      "hgnc_symbol", "entrezgene_id")
)

dup <-G_list[duplicated(G_list$hgnc_symbol)|duplicated(G_list$hgnc_symbol, fromLast=TRUE),]
not_empty <- dup[!(dup$hgnc_symbol == ""), ]
#of the 19537 duplicated samples, 16 582 have no gene symbols, 2955 have gene symbols 
#The genes have the same ensembl id, but different entrez_id 

G_list$entrezgene_id <- NULL

#Remove rows that does not have gene symbol 
G_list <- G_list %>% dplyr::filter(!(hgnc_symbol=="")) #16582 removed
dup <-G_list[duplicated(G_list$hgnc_symbol)|duplicated(G_list$hgnc_symbol, fromLast=TRUE),]
 

#order the genes in G-list_1 in the same order as in expr file 
expr <- rownames_to_column(expr, "id")
expr <- merge(expr, G_list, by.x = "id", by.y = "ensembl_gene_id")
expr <- expr %>% relocate("hgnc_symbol", .after = "id")

expr <- distinct(expr, id, .keep_all = TRUE)
#There are 32 cases where the ensembl ID is different but points to the same gene
#In most cases they have a count of 1 or similiar count
dup <-expr[duplicated(expr$hgnc_symbol)|duplicated(expr$hgnc_symbol, fromLast=TRUE),]

expr <- distinct(expr, hgnc_symbol, .keep_all = TRUE)
dup <-expr[duplicated(expr$hgnc_symbol)|duplicated(expr$hgnc_symbol, fromLast=TRUE),]
rm(dup, not_empty, raw_data)


#########################
# Create DGEList object #
#########################

d0 <- DGEList(expr)
d0

#Filtering
keep.exprs <- filterByExpr(d0, group= group) #By default it filters out genes that has total counts below 10 in the minimum number of samples
d0 <- d0[keep.exprs,, keep.lib.sizes=FALSE]
dim(d0)


#Normalizing the actual data now and visualizing
d0 <- calcNormFactors(d0, method = "TMM")
lcpm <- cpm(d0, log=TRUE)

#boxplot(lcpm,col=group)
library(RColorBrewer)
Shape <- c(0, 1, 2, 3, 4)
pam50 <- factor(c("LumA", "LumB", "Normal", "Her2", "Basal"))

col_sub <- brewer.pal(nlevels(pam50), "Set1")
names(col_sub) <- levels(pam50)

col_sub <- as.character(col_sub[pam50])
plotMDS(lcpm, col=col_sub, pch=Shape)
legend("bottomleft", col=c("#4DAF4A","#984EA3","#FF7F00","#377EB8","#E41A1C"),
       pch=c(0, 1, 2, 3, 4), legend=c("Luminal A", "Luminal B", "Normal-like", "HER2-enriched", "Basal-like"))
#luminal A = #4DAF4A
#Luminal B = #984EA3 
#Normal = #FF7F00
#Her2 = #377EB8
#Basal = #E41A1C 
#Levels: Basal Her2 LumA LumB Normal
#Levels: #E41A1C #377EB8 #4DAF4A #984EA3 #FF7F00


Shape_group <- c(0,1)
col_group = factor(c("High", "Low"))
myColors <- c("red", "blue")
levels(col_group) <- myColors
col_group <- as.character(col_group)
plotMDS(lcpm, col=col_group, pch=Shape_group)
legend("bottomleft", col=c("red","blue"),pch = c(0,1), legend=c("High TRIM45","Low TRIM45"))

#######
# PCA #
#######

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



design <- model.matrix(~0+group+subtype)

colnames(design) <- gsub("group","", colnames(design))
design


v <- voom(d0, design, plot=TRUE)
vfit <- lmFit(v, design)


contr <- makeContrasts(High - Low, levels = colnames(coef(vfit)))
contr


vfit <- contrasts.fit(vfit, contr)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
summary(decideTests(efit))

top.table <- topTable(efit, sort.by = "P", n = Inf, coef = 1, adjust.method = "BH")
write.csv2(top.table, "2025_03_25_limma_trim45_top_bottom_tcga_brca_er_positive_bh.csv")

top.table <- subset(top.table, top.table$hgnc_symbol != "TRIM45")

v <- EnhancedVolcano(top.table,
                lab = top.table$hgnc_symbol,
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
                colAlpha = 0.5)
v

pdf("2025_03_25_volcanoplot_limma_trim45_top_bottom_tcga_brca_er_positive_bh.pdf", heigh = 10, width =10)
print(v)
dev.off()

png("limma/limma_tcga_top_bottom_25/volcanoplot_limma_trim45_top_bottom_tcga_brca_er_positive_no_corr.png",res= 200, heigh = 2000, width =1800)
print(v)
dev.off()
